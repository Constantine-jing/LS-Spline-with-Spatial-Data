# ============================================================
# gibbs_interaction.R
# p=2 Proof-of-Concept: Collapsed Gibbs Sampler with Tensor-Product
# LS Interaction Surface
#
# Model:
#   y(s_r) = mu + f_1(X_1r) + f_2(X_2r) + f_{12}(X_1r, X_2r) + b(s_r) + eps_r
#
#   f_j   : LS natural cubic spline with RW2 prior (existing)
#   f_{12}: LS tensor-product interaction, W_uv * beta_uv, RW2 prior (NEW)
#   b     : Matern GP spatial effect, marginalized out (collapsed Gibbs)
#   eps   : iid N(0, tau2)
#
# Priors:
#   mu          ~ N(0, kappa2)
#   beta_j      ~ N(0, tau2_s_j * K_j^-)     j = 1, 2
#   beta_12     ~ N(0, tau2_s_12 * K_12^-)   (interaction, RW2 2D)
#   tau2_s_j    ~ IG(a_s, b_s)   (conjugate)
#   tau2_s_12   ~ IG(a_s, b_s)   (conjugate, separate smoothing variance)
#   sigma2      ~ IG(a_sigma, b_sigma)  via MH
#   tau2        ~ IG(a_tau, b_tau)      via MH
#   rho         ~ logN(log_rho_mu, log_rho_sd^2)  via MH
#
# Gibbs steps (identical structure to gibbs_stage_c_full.R):
#   1) eta = (mu, beta_1, beta_2, beta_12) | y, sigma2, tau2, rho, {tau2_s}  ~ MVN
#   2) sigma2 | rest   MH
#   3) tau2   | rest   MH
#   4) rho    | rest   MH
#   5) tau2_s_j | beta_j  ~ IG  (for j=1,2 and j=12)
#   6) b_hat (derived posterior mean, stored)
#
# Depends on:
#   ls_basis.R          (ls_build_one_full, ls_additive_build)
#   ls_interaction.R    (ls_build_interaction, khatri_rao_rowwise_R)
#   spatial_utils.R     (matern_cor)
#   ls_interaction_core.cpp  (optional, loaded via Rcpp::sourceCpp)
#
# Usage:
#   source("ls_basis.R"); source("ls_interaction.R"); source("spatial_utils.R")
#   source("gibbs_interaction.R")
#   # Optionally compile C++:
#   # Rcpp::sourceCpp("ls_interaction_core.cpp")
#   result <- run_interaction_poc(n=500, M=10, seed=1)
# ============================================================

# Try to load C++ functions (optional; R fallbacks used if unavailable)
.has_cpp <- tryCatch({
  requireNamespace("Rcpp", quietly = TRUE) &&
  exists("khatri_rao_cpp")
}, error = function(e) FALSE)

if (.has_cpp) message("Using C++ hot loops.") else
  message("Using pure R fallbacks (compile ls_interaction_core.cpp for speedup).")


# ============================================================
# build_interaction_prior_precision()
#
# Builds the FULL prior precision Q0 for eta = (mu, beta_1, beta_2, beta_12).
# Extends build_block_prior_precision() from gibbs_stage_c_full.R.
#
# col_map_full: list with
#   $main      : list of 0-indexed column ranges for main effect blocks
#   $interaction: list of 0-indexed column ranges for interaction blocks
# K_int_list: list of identified 2D RW2 penalty matrices (one per interaction pair)
# tau2_s: named/ordered vector (tau2_s_1, tau2_s_2, tau2_s_12)
# ============================================================
build_interaction_prior_precision <- function(col_map_full, K_main_list,
                                              K_int_list, tau2_s_main,
                                              tau2_s_int, kappa2 = 1e6,
                                              eps_ridge = 1e-6) {
  p_total <- 1L + sum(sapply(K_main_list, nrow)) + sum(sapply(K_int_list, nrow))
  Q0 <- matrix(0, p_total, p_total)

  # Intercept
  Q0[1, 1] <- 1 / kappa2

  # Main effect blocks
  for (j in seq_along(K_main_list)) {
    idx <- 1 + col_map_full$main[[j]]   # 0-indexed -> 1-indexed
    Kj  <- K_main_list[[j]]
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s_main[j]) * Kj + eps_ridge * diag(length(idx))
  }

  # Interaction blocks
  for (k in seq_along(K_int_list)) {
    idx <- 1 + col_map_full$interaction[[k]]
    Kk  <- K_int_list[[k]]
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s_int[k]) * Kk + eps_ridge * diag(length(idx))
  }

  Q0
}


# ============================================================
# gibbs_interaction_sampler()
#
# Main sampler. Arguments mirror gibbs_full_sampler() from
# gibbs_stage_c_full.R, with additions for interaction terms.
#
# New arguments:
#   W_int_list  : list of n x (M_u-1)*(M_v-1) interaction design matrices
#   K_int_list  : list of corresponding identified 2D RW2 penalty matrices
#   col_map_int : list of 0-indexed column ranges for each interaction block in H
# ============================================================
gibbs_interaction_sampler <- function(
    y, H, D, nu = 1.5,
    col_map_main,    # list: 0-indexed col ranges for main effects in H (after intercept)
    K_main_list,     # list: identified 1D RW2 penalty for each main effect
    col_map_int,     # list: 0-indexed col ranges for each interaction block in H
    K_int_list,      # list: identified 2D RW2 penalty for each interaction
    n_iter  = 5000,
    n_burn  = 1000,
    n_thin  = 1,
    kappa2  = 1e6,
    a_sigma = 2, b_sigma = 1,
    a_tau   = 2, b_tau   = 0.3,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = -1.6, log_rho_sd = 1.0,
    mh_sd_log_sigma2 = 0.3,
    mh_sd_log_tau2   = 0.3,
    mh_sd_log_rho    = 0.2,
    eps_ridge = 1e-6,
    init      = NULL,
    jitter    = 1e-8,
    verbose   = TRUE
) {
  y <- as.numeric(y)
  H <- as.matrix(H)
  D <- as.matrix(D)
  n <- length(y)
  p <- ncol(H)

  n_main <- length(col_map_main)
  n_int  <- length(col_map_int)

  # Pre-build col_map_full for Q0 constructor
  col_map_full <- list(main = col_map_main, interaction = col_map_int)

  # Helper: Matern correlation
  compute_R <- function(rho) matern_cor(D, rho = rho, nu = nu)

  # Log marginal likelihood (quadratic part only; log-det separated)
  log_marg_lik <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- tryCatch(chol(Sigma + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L)) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    alpha  <- forwardsolve(t(L), resid)
    -0.5 * (logdet + sum(alpha^2))
  }

  log_ig_prior <- function(x, a, b) -(a + 1) * log(x) - b / x

  log_lognormal_prior <- function(rho, mu, sd)
    dnorm(log(rho), mean = mu, sd = sd, log = TRUE) - log(rho)

  compute_b_postmean <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    Sigma_inv <- chol2inv(chol(Sigma + diag(jitter, n)))
    as.vector(sigma2 * R %*% Sigma_inv %*% resid)
  }

  # --- Initialization ---
  if (is.null(init)) {
    sigma2      <- 1.0
    tau2        <- 0.5
    rho         <- 0.2
    eta         <- rep(0, p)
    tau2_s_main <- rep(1.0, n_main)
    tau2_s_int  <- rep(1.0, n_int)
  } else {
    sigma2      <- init$sigma2
    tau2        <- init$tau2
    rho         <- init$rho
    eta         <- init$eta
    tau2_s_main <- init$tau2_s_main
    tau2_s_int  <- init$tau2_s_int
  }
  R <- compute_R(rho)

  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples      <- matrix(NA, n_keep, p)
  b_samples        <- matrix(NA, n_keep, n)
  sigma2_samples   <- numeric(n_keep)
  tau2_samples     <- numeric(n_keep)
  rho_samples      <- numeric(n_keep)
  tau2_s_main_samp <- matrix(NA, n_keep, n_main)
  tau2_s_int_samp  <- matrix(NA, n_keep, n_int)

  accept <- c(sigma2 = 0, tau2 = 0, rho = 0)
  keep_idx <- 0

  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in seq_len(n_iter)) {

    # ----------------------------------------------------------
    # Step 1: Draw eta | y, sigma2, tau2, rho, {tau2_s}
    #
    # eta = (mu, beta_1, beta_2, beta_12) ~ MVN
    # Q_eta = H^T Sigma^{-1} H + Q0
    # m_eta = Q_eta^{-1} H^T Sigma^{-1} y
    # ----------------------------------------------------------
    Q0 <- build_interaction_prior_precision(
      col_map_full, K_main_list, K_int_list,
      tau2_s_main, tau2_s_int, kappa2, eps_ridge
    )

    Sigma    <- sigma2 * R + tau2 * diag(n)
    L_Sigma  <- chol(Sigma + diag(jitter, n))
    y_w      <- forwardsolve(t(L_Sigma), y)
    H_w      <- forwardsolve(t(L_Sigma), H)

    Q_eta <- crossprod(H_w) + Q0
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% crossprod(H_w, y_w))

    z   <- rnorm(p)
    eta <- m_eta + as.vector(backsolve(U_eta, z))

    resid <- as.vector(y - H %*% eta)

    # ----------------------------------------------------------
    # Steps 2-4: MH for sigma2, tau2, rho  (identical to gibbs_stage_c_full.R)
    # ----------------------------------------------------------

    # sigma2
    sigma2_prop <- exp(log(sigma2) + rnorm(1, 0, mh_sd_log_sigma2))
    lp_c <- log_marg_lik(resid, sigma2,      tau2, R) + log_ig_prior(sigma2,      a_sigma, b_sigma)
    lp_p <- log_marg_lik(resid, sigma2_prop, tau2, R) + log_ig_prior(sigma2_prop, a_sigma, b_sigma)
    if (log(runif(1)) < (lp_p - lp_c)) { sigma2 <- sigma2_prop; accept["sigma2"] <- accept["sigma2"] + 1 }

    # tau2
    tau2_prop <- exp(log(tau2) + rnorm(1, 0, mh_sd_log_tau2))
    lp_c <- log_marg_lik(resid, sigma2, tau2,      R) + log_ig_prior(tau2,      a_tau, b_tau)
    lp_p <- log_marg_lik(resid, sigma2, tau2_prop, R) + log_ig_prior(tau2_prop, a_tau, b_tau)
    if (log(runif(1)) < (lp_p - lp_c)) { tau2 <- tau2_prop; accept["tau2"] <- accept["tau2"] + 1 }

    # rho
    rho_prop <- exp(log(rho) + rnorm(1, 0, mh_sd_log_rho))
    R_prop   <- compute_R(rho_prop)
    lp_c <- log_marg_lik(resid, sigma2, tau2, R)      + log_lognormal_prior(rho,      log_rho_mu, log_rho_sd)
    lp_p <- log_marg_lik(resid, sigma2, tau2, R_prop) + log_lognormal_prior(rho_prop, log_rho_mu, log_rho_sd)
    if (log(runif(1)) < (lp_p - lp_c)) { rho <- rho_prop; R <- R_prop; accept["rho"] <- accept["rho"] + 1 }

    # ----------------------------------------------------------
    # Step 5a: Draw tau2_s_j | beta_j  (main effects, conjugate IG)
    # ----------------------------------------------------------
    for (j in seq_len(n_main)) {
      idx    <- 1 + col_map_main[[j]]
      beta_j <- eta[idx]
      K_j    <- K_main_list[[j]]
      qf     <- as.numeric(t(beta_j) %*% K_j %*% beta_j)
      rank_j <- nrow(K_j) - 2
      if (rank_j < 1) rank_j <- 1
      tau2_s_main[j] <- 1 / rgamma(1, shape = a_smooth + rank_j / 2,
                                       rate  = b_smooth + qf / 2)
    }

    # ----------------------------------------------------------
    # Step 5b: Draw tau2_s_12 | beta_12  (interaction, conjugate IG)
    # ----------------------------------------------------------
    for (k in seq_len(n_int)) {
      idx    <- 1 + col_map_int[[k]]
      beta_k <- eta[idx]
      K_k    <- K_int_list[[k]]
      qf     <- as.numeric(t(beta_k) %*% K_k %*% beta_k)
      # Rank of 2D RW2 penalty: (M_u-1)*(M_v-1) - 1
      rank_k <- nrow(K_k) - 1
      if (rank_k < 1) rank_k <- 1
      tau2_s_int[k] <- 1 / rgamma(1, shape = a_smooth + rank_k / 2,
                                      rate  = b_smooth + qf / 2)
    }

    # ----------------------------------------------------------
    # Step 6: Posterior mean of b (derived)
    # ----------------------------------------------------------
    b_mean <- compute_b_postmean(resid, sigma2, tau2, R)

    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]          <- eta
      b_samples[keep_idx, ]            <- b_mean
      sigma2_samples[keep_idx]         <- sigma2
      tau2_samples[keep_idx]           <- tau2
      rho_samples[keep_idx]            <- rho
      tau2_s_main_samp[keep_idx, ]     <- tau2_s_main
      tau2_s_int_samp[keep_idx, ]      <- tau2_s_int
    }

    if (verbose && (iter %% 500 == 0)) {
      cat(sprintf("  iter %d/%d  rho=%.3f  sigma2=%.3f  tau2=%.3f  tau2_s_int=%.4f\n",
                  iter, n_iter, rho, sigma2, tau2, tau2_s_int[1]))
    }
  }

  cat(sprintf("  Accept rates: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              accept["sigma2"] / n_iter, accept["tau2"] / n_iter, accept["rho"] / n_iter))

  list(
    eta_samples      = eta_samples,
    b_samples        = b_samples,
    sigma2_samples   = sigma2_samples,
    tau2_samples     = tau2_samples,
    rho_samples      = rho_samples,
    tau2_s_main_samp = tau2_s_main_samp,
    tau2_s_int_samp  = tau2_s_int_samp,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
    accept_rate = accept / n_iter
  )
}


# ============================================================
# run_interaction_poc()
#
# Self-contained proof-of-concept simulation with p=2 covariates
# and a KNOWN interaction surface. Runs the sampler and produces
# summary diagnostics.
#
# True model:
#   f_1(x) = sin(2*pi*x)
#   f_2(x) = (2*x - 1)^2
#   f_12(x1, x2) = 1.5 * sin(pi*x1) * cos(pi*x2)   [interaction surface]
#   b(s): Matern(rho=0.3, nu=1.5), sigma2=0.5
#   eps ~ N(0, tau2=0.25)
#
# Locations: n uniformly scattered on [0,1]^2
# ============================================================
run_interaction_poc <- function(n = 500, M = 10, n_iter = 3000, n_burn = 500,
                                 seed = 42, verbose = TRUE) {
  set.seed(seed)
  cat("=== Interaction POC: p=2, n=", n, ", M=", M, "===\n")

  # ---- 1. Simulate data ----
  # Locations
  locs <- matrix(runif(2 * n), n, 2)
  D    <- as.matrix(dist(locs))

  # Covariates (independent of location for POC)
  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)

  # True functions
  true_f1  <- sin(2 * pi * X1)
  true_f2  <- (2 * X2 - 1)^2
  true_f12 <- 1.5 * sin(pi * X1) * cos(pi * X2)

  # Spatial effect (Matern GP)
  true_sigma2 <- 0.5
  true_rho    <- 0.3
  true_tau2   <- 0.25
  true_nu     <- 1.5

  R_true <- matern_cor(D, rho = true_rho, nu = true_nu)
  L_R    <- chol(R_true + diag(1e-8, n))
  b_true <- as.vector(t(L_R) %*% rnorm(n)) * sqrt(true_sigma2)
  eps    <- rnorm(n, 0, sqrt(true_tau2))

  y <- true_f1 + true_f2 + true_f12 + b_true + eps
  cat(sprintf("  y: mean=%.3f, sd=%.3f\n", mean(y), sd(y)))

  # ---- 2. Build LS bases ----
  obj1 <- ls_build_one_full(X1, M = M)   # main effect for X1
  obj2 <- ls_build_one_full(X2, M = M)   # main effect for X2

  W1 <- obj1$W   # n x (M-1)
  W2 <- obj2$W   # n x (M-1)

  # 1D RW2 penalties (on identified beta)
  T1 <- obj1$T; T2 <- obj2$T
  K1_raw  <- build_rw2_penalty_1d(M)
  K2_raw  <- build_rw2_penalty_1d(M)
  K1_beta <- t(T1) %*% K1_raw %*% T1   # (M-1) x (M-1)
  K2_beta <- t(T2) %*% K2_raw %*% T2

  # ---- 3. Build interaction basis ----
  int12 <- ls_build_interaction(obj1, obj2)
  W12   <- int12$W_uv    # n x (M-1)^2
  K12   <- int12$K_uv    # (M-1)^2 x (M-1)^2

  cat(sprintf("  Bases: W1=%dx%d, W2=%dx%d, W12=%dx%d\n",
              nrow(W1), ncol(W1), nrow(W2), ncol(W2), nrow(W12), ncol(W12)))

  # ---- 4. Assemble H and col_map ----
  # H = [1 | W1 | W2 | W12]
  H <- cbind(1, W1, W2, W12)
  d1 <- ncol(W1); d2 <- ncol(W2); d12 <- ncol(W12)

  # 1-indexed column ranges in H (H = [1 | W1 | W2 | W12], intercept at col 1)
  col_map_main <- list(
    seq(2,        d1 + 1),                    # W1 block
    seq(d1 + 2,   d1 + d2 + 1)                # W2 block
  )
  col_map_int <- list(
    seq(d1 + d2 + 2, d1 + d2 + d12 + 1)       # W12 block
  )
  
  # Convert to 0-indexed (sampler uses 1 + idx internally)
  col_map_main <- lapply(col_map_main, function(x) x - 1L)
  col_map_int  <- lapply(col_map_int,  function(x) x - 1L)

  K_main_list <- list(K1_beta, K2_beta)
  K_int_list  <- list(K12)

  cat(sprintf("  H: %d x %d\n", nrow(H), ncol(H)))

  # ---- 5. Run sampler ----
  t_start <- proc.time()
  gs <- gibbs_interaction_sampler(
    y = y, H = H, D = D, nu = true_nu,
    col_map_main = col_map_main,
    K_main_list  = K_main_list,
    col_map_int  = col_map_int,
    K_int_list   = K_int_list,
    n_iter = n_iter, n_burn = n_burn,
    verbose = verbose
  )
  t_elapsed <- (proc.time() - t_start)["elapsed"]
  cat(sprintf("  Elapsed: %.1f sec\n", t_elapsed))

  # ---- 6. Posterior summaries ----
  eta_mean <- colMeans(gs$eta_samples)

  # Recover fitted curves
  f1_hat  <- W1  %*% eta_mean[2:(d1 + 1)]
  f2_hat  <- W2  %*% eta_mean[(d1 + 2):(d1 + d2 + 1)]
  f12_hat <- W12 %*% eta_mean[(d1 + d2 + 2):(d1 + d2 + d12 + 1)]

  rmse_f1  <- sqrt(mean((f1_hat  - true_f1)^2))
  rmse_f2  <- sqrt(mean((f2_hat  - true_f2)^2))
  rmse_f12 <- sqrt(mean((f12_hat - true_f12)^2))

  cat(sprintf("  RMSE f1=%.4f  f2=%.4f  f12=%.4f\n", rmse_f1, rmse_f2, rmse_f12))
  cat(sprintf("  sigma2: mean=%.3f (true=%.2f)\n",
              mean(gs$sigma2_samples), true_sigma2))
  cat(sprintf("  tau2:   mean=%.3f (true=%.2f)\n",
              mean(gs$tau2_samples), true_tau2))
  cat(sprintf("  rho:    mean=%.3f (true=%.2f)\n",
              mean(gs$rho_samples), true_rho))
  cat(sprintf("  tau2_s_int (interaction smoothing): mean=%.4f\n",
              mean(gs$tau2_s_int_samp)))

  invisible(list(
    gs = gs, H = H, y = y, X1 = X1, X2 = X2,
    true_f1 = true_f1, true_f2 = true_f2, true_f12 = true_f12,
    f1_hat = f1_hat, f2_hat = f2_hat, f12_hat = f12_hat,
    rmse = c(f1 = rmse_f1, f2 = rmse_f2, f12 = rmse_f12),
    obj1 = obj1, obj2 = obj2, int12 = int12,
    col_map_main = col_map_main, col_map_int = col_map_int
  ))
}


# ============================================================
# plot_interaction_poc()
#
# Diagnostics and surface plots for the POC result.
# Saved to PDF (no X11 needed on Hellbender).
# ============================================================
plot_interaction_poc <- function(poc, outfile = "interaction_poc_plots.pdf") {
  pdf(outfile, width = 10, height = 8)

  gs <- poc$gs

  # Trace plots
  par(mfrow = c(3, 2))
  plot(gs$sigma2_samples, type = "l", main = "sigma2 trace", ylab = "sigma2")
  abline(h = 0.5, col = "red", lty = 2)
  plot(gs$tau2_samples,   type = "l", main = "tau2 trace",   ylab = "tau2")
  abline(h = 0.25, col = "red", lty = 2)
  plot(gs$rho_samples,    type = "l", main = "rho trace",    ylab = "rho")
  abline(h = 0.3, col = "red", lty = 2)
  plot(gs$tau2_s_int_samp[, 1], type = "l",
       main = "tau2_s_12 trace (interaction smoothing)", ylab = "tau2_s_12")
  plot(gs$tau2_s_main_samp[, 1], type = "l",
       main = "tau2_s_1 trace", ylab = "tau2_s_1")
  plot(gs$tau2_s_main_samp[, 2], type = "l",
       main = "tau2_s_2 trace", ylab = "tau2_s_2")

  # Main effect fits
  par(mfrow = c(1, 2))
  ord1 <- order(poc$X1)
  plot(poc$X1[ord1], poc$true_f1[ord1], type = "l", lwd = 2, col = "red",
       main = "f1(X1): true vs estimated", xlab = "X1", ylab = "f1")
  lines(poc$X1[ord1], poc$f1_hat[ord1], col = "steelblue", lwd = 2, lty = 2)
  legend("topright", c("True", "Posterior mean"), col = c("red","steelblue"),
         lwd = 2, lty = 1:2)

  ord2 <- order(poc$X2)
  plot(poc$X2[ord2], poc$true_f2[ord2], type = "l", lwd = 2, col = "red",
       main = "f2(X2): true vs estimated", xlab = "X2", ylab = "f2")
  lines(poc$X2[ord2], poc$f2_hat[ord2], col = "darkorange", lwd = 2, lty = 2)
  legend("topright", c("True", "Posterior mean"), col = c("red","darkorange"),
         lwd = 2, lty = 1:2)

  # Interaction surface: scatter plot of true vs estimated (coloured by value)
  par(mfrow = c(1, 2))
  zlim <- range(c(poc$true_f12, poc$f12_hat))
  col_ramp <- colorRampPalette(c("blue","white","red"))(100)
  col_true <- col_ramp[cut(poc$true_f12, 100, labels = FALSE)]
  col_hat  <- col_ramp[cut(poc$f12_hat,  100, labels = FALSE)]
  plot(poc$X1, poc$X2, col = col_true, pch = 16, cex = 0.5,
       main = "True f_{12}(X1,X2)", xlab = "X1", ylab = "X2")
  plot(poc$X1, poc$X2, col = col_hat,  pch = 16, cex = 0.5,
       main = "Estimated f_{12}(X1,X2)", xlab = "X1", ylab = "X2")

  dev.off()
  cat(sprintf("  Plots saved to %s\n", outfile))
}
