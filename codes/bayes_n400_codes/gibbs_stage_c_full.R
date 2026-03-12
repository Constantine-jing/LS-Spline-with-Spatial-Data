# ============================================================
# gibbs_stage_c_full.R
# Full Bayesian model: collapsed Gibbs + RW2 smoothing + rho sampling
#
# Everything from gibbs_stage_c.R PLUS:
#   - rho (Matern range) sampled via MH on log(rho)
#   - When rho changes, R must be recomputed
#
# Model:
#   y = H eta + b + eps
#   b ~ N(0, sigma2 * R(rho, nu))     [marginalized out]
#   eps ~ N(0, tau2 * I)
#
# Priors:
#   mu ~ N(0, kappa2)                  [vague]
#   beta_j ~ N(0, tau2_s_j * K_j^-)   [RW2 smoothing]
#   tau2_s_j ~ IG(a_s, b_s)
#   sigma2 ~ IG(a_sigma, b_sigma)
#   tau2   ~ IG(a_tau, b_tau)
#   rho    ~ logN(log_rho_mu, log_rho_sd)  [log-normal prior]
#   nu     : FIXED (standard practice)
#
# Gibbs steps:
#   1) eta      | y, sigma2, tau2, rho, {tau2_s_j}  ~ MVN
#   2) sigma2   | y, eta, tau2, rho                  via MH
#   3) tau2     | y, eta, sigma2, rho                via MH
#   4) rho      | y, eta, sigma2, tau2               via MH  [NEW]
#   5) tau2_s_j | beta_j                             ~ IG
#   6) b        (derived posterior mean)
#
# Depends on: spatial_utils.R, ls_basis.R
# ============================================================


# --- Re-use penalty builders from gibbs_stage_c.R ---
# (source that file, or copy these if standalone)

build_rw2_penalty <- function(d) {
  if (d < 3) return(diag(d) * 0)
  D2 <- matrix(0, d - 2, d)
  for (i in 1:(d - 2)) {
    D2[i, i]     <-  1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <-  1
  }
  crossprod(D2)
}

build_block_prior_precision <- function(col_map, tau2_s, kappa2 = 1e6,
                                         eps_ridge = 1e-6) {
  p_covs <- length(col_map)
  p_total <- 1 + max(unlist(col_map))
  Q0 <- matrix(0, p_total, p_total)
  Q0[1, 1] <- 1 / kappa2
  for (j in 1:p_covs) {
    idx <- 1 + col_map[[j]]
    d_j <- length(idx)
    K_j <- build_rw2_penalty(d_j)
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s[j]) * K_j + eps_ridge * diag(d_j)
  }
  Q0
}


# ------------------------------------------------------------
# gibbs_full_sampler()
#
# Full Bayesian sampler with rho sampling.
# ------------------------------------------------------------
gibbs_full_sampler <- function(y, H, D, nu = 1.5, col_map,
                                n_iter  = 6000,
                                n_burn  = 1500,
                                n_thin  = 1,
                                kappa2  = 1e6,
                                # IG priors on sigma2, tau2
                                a_sigma = 2, b_sigma = 1,
                                a_tau   = 2, b_tau   = 0.3,
                                # IG priors on smoothing variances
                                a_smooth = 1, b_smooth = 0.005,
                                # Log-normal prior on rho:
                                # log(rho) ~ N(log_rho_mu, log_rho_sd^2)
                                log_rho_mu = -1.6,  # ~ exp(-1.6) = 0.20
                                log_rho_sd = 1.0,   # fairly diffuse
                                # MH tuning
                                mh_sd_log_sigma2 = 0.3,
                                mh_sd_log_tau2   = 0.3,
                                mh_sd_log_rho    = 0.2,
                                eps_ridge = 1e-6,
                                init    = NULL,
                                jitter  = 1e-8,
                                verbose = TRUE) {
  
  y <- as.numeric(y)
  H <- as.matrix(H)
  D <- as.matrix(D)  # pairwise distance matrix (fixed)
  n <- length(y)
  p <- ncol(H)
  p_covs <- length(col_map)
  
  # Pre-compute K_j matrices
  K_list <- vector("list", p_covs)
  d_list <- integer(p_covs)
  for (j in 1:p_covs) {
    d_j <- length(col_map[[j]])
    K_list[[j]] <- build_rw2_penalty(d_j)
    d_list[j] <- d_j
  }
  
  # --- Helpers ---
  
  # Compute R from rho (called whenever rho changes)
  compute_R <- function(rho) {
    matern_cor(D, rho = rho, nu = nu)
  }
  
  # Log marginal likelihood: y - H*eta ~ N(0, sigma2*R + tau2*I)
  log_marg_lik <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- tryCatch(chol(Sigma + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L)) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    alpha <- forwardsolve(t(L), resid)
    -0.5 * (logdet + sum(alpha^2))
  }
  
  log_ig_prior <- function(x, a, b) {
    -(a + 1) * log(x) - b / x
  }
  
  # Log-normal prior on rho: log(rho) ~ N(mu, sd^2)
  log_lognormal_prior <- function(rho, mu, sd) {
    dnorm(log(rho), mean = mu, sd = sd, log = TRUE) - log(rho)
  }
  
  compute_b_postmean <- function(resid, sigma2, tau2, R) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    Sigma_inv <- chol2inv(chol(Sigma + diag(jitter, n)))
    as.vector(sigma2 * R %*% Sigma_inv %*% resid)
  }
  
  # --- Initialize ---
  if (is.null(init)) {
    sigma2 <- 1.0
    tau2   <- 0.5
    rho    <- 0.2
    eta    <- rep(0, p)
    tau2_s <- rep(1.0, p_covs)
  } else {
    sigma2 <- init$sigma2
    tau2   <- init$tau2
    rho    <- if (!is.null(init$rho)) init$rho else 0.2
    eta    <- init$eta
    tau2_s <- if (!is.null(init$tau2_s)) init$tau2_s else rep(1.0, p_covs)
  }
  
  # Compute initial R
  R <- compute_R(rho)
  
  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples    <- matrix(NA, n_keep, p)
  b_samples      <- matrix(NA, n_keep, n)
  sigma2_samples <- numeric(n_keep)
  tau2_samples   <- numeric(n_keep)
  rho_samples    <- numeric(n_keep)
  tau2_s_samples <- matrix(NA, n_keep, p_covs)
  
  accept_sigma2 <- 0
  accept_tau2   <- 0
  accept_rho    <- 0
  keep_idx      <- 0
  
  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in 1:n_iter) {
    
    # ----------------------------------------------------------
    # Step 1: Draw eta | y, sigma2, tau2, rho, {tau2_s_j}
    # ----------------------------------------------------------
    Q0 <- build_block_prior_precision(col_map, tau2_s, kappa2, eps_ridge)
    
    Sigma <- sigma2 * R + tau2 * diag(n)
    L_Sigma <- chol(Sigma + diag(jitter, n))
    y_w <- forwardsolve(t(L_Sigma), y)
    H_w <- forwardsolve(t(L_Sigma), H)
    
    Q_eta <- crossprod(H_w) + Q0
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% crossprod(H_w, y_w))
    
    z <- rnorm(p)
    eta <- m_eta + as.vector(backsolve(U_eta, z))
    
    resid <- as.vector(y - H %*% eta)
    
    # ----------------------------------------------------------
    # Step 2: MH for sigma2 | y, eta, tau2, rho
    # ----------------------------------------------------------
    log_sigma2_prop <- log(sigma2) + rnorm(1, 0, mh_sd_log_sigma2)
    sigma2_prop <- exp(log_sigma2_prop)
    
    lp_curr <- log_marg_lik(resid, sigma2, tau2, R) +
               log_ig_prior(sigma2, a_sigma, b_sigma)
    lp_prop <- log_marg_lik(resid, sigma2_prop, tau2, R) +
               log_ig_prior(sigma2_prop, a_sigma, b_sigma)
    
    if (log(runif(1)) < (lp_prop - lp_curr)) {
      sigma2 <- sigma2_prop
      accept_sigma2 <- accept_sigma2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 3: MH for tau2 | y, eta, sigma2, rho
    # ----------------------------------------------------------
    log_tau2_prop <- log(tau2) + rnorm(1, 0, mh_sd_log_tau2)
    tau2_prop <- exp(log_tau2_prop)
    
    lp_curr <- log_marg_lik(resid, sigma2, tau2, R) +
               log_ig_prior(tau2, a_tau, b_tau)
    lp_prop <- log_marg_lik(resid, sigma2, tau2_prop, R) +
               log_ig_prior(tau2_prop, a_tau, b_tau)
    
    if (log(runif(1)) < (lp_prop - lp_curr)) {
      tau2 <- tau2_prop
      accept_tau2 <- accept_tau2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 4: MH for rho | y, eta, sigma2, tau2       [NEW]
    #
    # Propose on log scale: log(rho*) = log(rho) + N(0, mh_sd^2)
    # When rho changes, must recompute R.
    # Prior: log(rho) ~ N(log_rho_mu, log_rho_sd^2)
    # ----------------------------------------------------------
    log_rho_prop <- log(rho) + rnorm(1, 0, mh_sd_log_rho)
    rho_prop <- exp(log_rho_prop)
    
    # Recompute R at proposed rho
    R_prop <- compute_R(rho_prop)
    
    lp_curr_rho <- log_marg_lik(resid, sigma2, tau2, R) +
                   log_lognormal_prior(rho, log_rho_mu, log_rho_sd)
    lp_prop_rho <- log_marg_lik(resid, sigma2, tau2, R_prop) +
                   log_lognormal_prior(rho_prop, log_rho_mu, log_rho_sd)
    
    if (log(runif(1)) < (lp_prop_rho - lp_curr_rho)) {
      rho <- rho_prop
      R   <- R_prop    # update R for subsequent steps
      accept_rho <- accept_rho + 1
    }
    
    # ----------------------------------------------------------
    # Step 5: Draw tau2_s_j | beta_j  (conjugate IG)
    # ----------------------------------------------------------
    for (j in 1:p_covs) {
      idx <- 1 + col_map[[j]]
      beta_j <- eta[idx]
      K_j <- K_list[[j]]
      
      qform <- as.numeric(t(beta_j) %*% K_j %*% beta_j)
      rank_Kj <- d_list[j] - 2
      if (rank_Kj < 1) rank_Kj <- 1
      
      a_post <- a_smooth + rank_Kj / 2
      b_post <- b_smooth + qform / 2
      
      tau2_s[j] <- 1 / rgamma(1, shape = a_post, rate = b_post)
    }
    
    # ----------------------------------------------------------
    # Step 6: Compute b posterior mean (derived)
    # ----------------------------------------------------------
    b_mean <- compute_b_postmean(resid, sigma2, tau2, R)
    
    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]    <- eta
      b_samples[keep_idx, ]      <- b_mean
      sigma2_samples[keep_idx]   <- sigma2
      tau2_samples[keep_idx]     <- tau2
      rho_samples[keep_idx]      <- rho
      tau2_s_samples[keep_idx, ] <- tau2_s
    }
    
    if (verbose && (iter %% 1000 == 0)) {
      cat(sprintf("  iter %d/%d  rho=%.4f  sigma2=%.4f  tau2=%.4f  tau2_s=%s\n",
                  iter, n_iter, rho, sigma2, tau2,
                  paste(sprintf("%.3f", tau2_s), collapse=",")))
    }
  }
  
  cat(sprintf("  MH acceptance: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              accept_sigma2 / n_iter, accept_tau2 / n_iter,
              accept_rho / n_iter))
  
  list(
    eta_samples    = eta_samples,
    b_samples      = b_samples,
    sigma2_samples = sigma2_samples,
    tau2_samples   = tau2_samples,
    rho_samples    = rho_samples,
    tau2_s_samples = tau2_s_samples,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
    accept_rate = c(sigma2 = accept_sigma2 / n_iter,
                    tau2   = accept_tau2 / n_iter,
                    rho    = accept_rho / n_iter)
  )
}


# ------------------------------------------------------------
# plot_gibbs_trace_full()
#
# Trace + density for sigma2, tau2, AND rho.
# ------------------------------------------------------------
plot_gibbs_trace_full <- function(gs, true_sigma2 = NULL, true_tau2 = NULL,
                                   true_rho = NULL) {
  
  par(mfrow = c(3, 2))
  
  # sigma2
  plot(gs$sigma2_samples, type = "l", col = "steelblue",
       main = expression(paste("Trace: ", sigma^2)),
       ylab = expression(sigma^2), xlab = "iteration")
  if (!is.null(true_sigma2)) abline(h = true_sigma2, col = "red", lty = 2, lwd = 2)
  
  plot(density(gs$sigma2_samples), col = "steelblue", lwd = 2,
       main = expression(paste("Density: ", sigma^2)))
  if (!is.null(true_sigma2)) abline(v = true_sigma2, col = "red", lty = 2, lwd = 2)
  
  # tau2
  plot(gs$tau2_samples, type = "l", col = "darkorange",
       main = expression(paste("Trace: ", tau^2)),
       ylab = expression(tau^2), xlab = "iteration")
  if (!is.null(true_tau2)) abline(h = true_tau2, col = "red", lty = 2, lwd = 2)
  
  plot(density(gs$tau2_samples), col = "darkorange", lwd = 2,
       main = expression(paste("Density: ", tau^2)))
  if (!is.null(true_tau2)) abline(v = true_tau2, col = "red", lty = 2, lwd = 2)
  
  # rho
  plot(gs$rho_samples, type = "l", col = "forestgreen",
       main = expression(paste("Trace: ", rho)),
       ylab = expression(rho), xlab = "iteration")
  if (!is.null(true_rho)) abline(h = true_rho, col = "red", lty = 2, lwd = 2)
  
  plot(density(gs$rho_samples), col = "forestgreen", lwd = 2,
       main = expression(paste("Density: ", rho)))
  if (!is.null(true_rho)) abline(v = true_rho, col = "red", lty = 2, lwd = 2)
}
