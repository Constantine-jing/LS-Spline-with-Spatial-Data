# ============================================================
# run_interaction_null_test.R
#
# NULL INTERACTION TEST.
#
# Goal: verify that when the TRUE interaction surface is f_12 == 0,
# the sampler:
#   (a) recovers a posterior mean f_12 surface near zero,
#   (b) produces 95% credible bands for f_12 that cover 0 across the grid,
#   (c) shrinks tau2_s_int strongly toward 0 (or to the prior),
#   (d) does NOT degrade main-effect recovery vs. a sanity main-effects fit.
#
# This is the highest-stakes correctness check for the interaction
# extension: a positive false-detection of interaction would invalidate
# all downstream applications.
#
# DGP (modified from run_interaction_poc):
#   f_1(x)  = sin(2*pi*x)
#   f_2(x)  = (2*x - 1)^2
#   f_12    = 0                          [NULL]
#   b(s) ~ Matern(rho=0.3, nu=1.5), sigma2=0.5
#   eps  ~ N(0, tau2=0.25)
#
# Output:
#   - Console: RMSE/coverage/credible-band summary
#   - PDF:     trace plots, f_12 posterior surface, pointwise bands,
#              tau2_s_int trace + posterior density vs. prior
# ============================================================

# Required sources (assumed available in working directory)
# source("spatial_utils.R")
# source("ls_basis.R")
# source("ls_interaction.R")
# source("gibbs_interaction.R")

run_interaction_null_test <- function(n        = 300,
                                       M        = 8,
                                       n_iter   = 3000,
                                       n_burn   = 1000,
                                       seed     = 42,
                                       grid_len = 25,
                                       outfile  = "interaction_null_test.pdf",
                                       verbose  = TRUE) {

  set.seed(seed)
  cat("=== Interaction NULL test: f_12 == 0,  n =", n, " M =", M, "===\n")

  # ---- 1. Simulate data with NULL interaction ----
  locs <- matrix(runif(2 * n), n, 2)
  D    <- as.matrix(dist(locs))

  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)

  true_f1  <- sin(2 * pi * X1)
  true_f2  <- (2 * X2 - 1)^2
  true_f12 <- rep(0, n)               # <-- NULL

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

  # ---- 2. Build LS bases (identical to POC) ----
  obj1 <- ls_build_one_full(X1, M = M)
  obj2 <- ls_build_one_full(X2, M = M)

  W1 <- obj1$W
  W2 <- obj2$W
  T1 <- obj1$T; T2 <- obj2$T

  K1_raw  <- build_rw2_penalty_1d(M)
  K2_raw  <- build_rw2_penalty_1d(M)
  K1_beta <- t(T1) %*% K1_raw %*% T1
  K2_beta <- t(T2) %*% K2_raw %*% T2

  int12 <- ls_build_interaction(obj1, obj2)
  W12   <- int12$W_uv
  K12   <- int12$K_uv

  cat(sprintf("  Bases: W1=%dx%d, W2=%dx%d, W12=%dx%d\n",
              nrow(W1), ncol(W1), nrow(W2), ncol(W2), nrow(W12), ncol(W12)))

  # ---- 3. Assemble H and col_map (uses fixed indexing) ----
  H  <- cbind(1, W1, W2, W12)
  d1 <- ncol(W1); d2 <- ncol(W2); d12 <- ncol(W12)

  # 1-indexed ranges in H
  col_map_main <- list(
    seq(2,        d1 + 1),
    seq(d1 + 2,   d1 + d2 + 1)
  )
  col_map_int <- list(
    seq(d1 + d2 + 2, d1 + d2 + d12 + 1)
  )
  # Convert to 0-indexed (sampler convention)
  col_map_main <- lapply(col_map_main, function(x) x - 1L)
  col_map_int  <- lapply(col_map_int,  function(x) x - 1L)

  K_main_list <- list(K1_beta, K2_beta)
  K_int_list  <- list(K12)

  cat(sprintf("  H: %d x %d\n", nrow(H), ncol(H)))

  # ---- 4. Run sampler ----
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

  # ---- 5. Posterior summaries: in-sample f_12 ----
  # Indices in eta of each block
  idx_mu   <- 1
  idx_b1   <- 2:(d1 + 1)
  idx_b2   <- (d1 + 2):(d1 + d2 + 1)
  idx_b12  <- (d1 + d2 + 2):(d1 + d2 + d12 + 1)

  beta12_post <- gs$eta_samples[, idx_b12, drop = FALSE]   # n_keep x d12
  # In-sample f_12 surface samples: n_keep x n
  f12_samp_in <- beta12_post %*% t(W12)
  f12_mean_in <- colMeans(f12_samp_in)
  f12_lo_in   <- apply(f12_samp_in, 2, quantile, 0.025)
  f12_hi_in   <- apply(f12_samp_in, 2, quantile, 0.975)

  rmse_f12_in   <- sqrt(mean(f12_mean_in^2))                 # truth = 0
  cov_f12_in    <- mean(f12_lo_in <= 0 & 0 <= f12_hi_in)     # should be ~0.95
  band_width_in <- mean(f12_hi_in - f12_lo_in)

  # ---- 6. Posterior summaries: f_12 on a regular grid ----
  g  <- seq(0, 1, length.out = grid_len)
  gg <- expand.grid(X1 = g, X2 = g)
  W12_g <- ls_interaction_design_new(gg$X1, gg$X2, int12$recipe, clip = TRUE)

  f12_samp_grid <- beta12_post %*% t(W12_g)                  # n_keep x (grid_len^2)
  f12_mean_grid <- colMeans(f12_samp_grid)
  f12_lo_grid   <- apply(f12_samp_grid, 2, quantile, 0.025)
  f12_hi_grid   <- apply(f12_samp_grid, 2, quantile, 0.975)

  rmse_f12_grid   <- sqrt(mean(f12_mean_grid^2))
  cov_f12_grid    <- mean(f12_lo_grid <= 0 & 0 <= f12_hi_grid)
  band_width_grid <- mean(f12_hi_grid - f12_lo_grid)

  # ---- 7. Main-effect recovery (sanity: should be unaffected) ----
  beta1_post <- gs$eta_samples[, idx_b1, drop = FALSE]
  beta2_post <- gs$eta_samples[, idx_b2, drop = FALSE]
  f1_hat <- W1 %*% colMeans(beta1_post)
  f2_hat <- W2 %*% colMeans(beta2_post)
  rmse_f1 <- sqrt(mean((f1_hat - true_f1)^2))
  rmse_f2 <- sqrt(mean((f2_hat - true_f2)^2))

  # ---- 8. Variance components ----
  sigma2_mean <- mean(gs$sigma2_samples)
  tau2_mean   <- mean(gs$tau2_samples)
  rho_mean    <- mean(gs$rho_samples)
  tau2_s_int_mean   <- mean(gs$tau2_s_int_samp)
  tau2_s_int_median <- median(gs$tau2_s_int_samp)
  tau2_s_int_q025   <- quantile(gs$tau2_s_int_samp, 0.025)
  tau2_s_int_q975   <- quantile(gs$tau2_s_int_samp, 0.975)

  # ---- 9. Print summary ----
  cat("\n--- NULL test summary ---\n")
  cat(sprintf("  In-sample f_12 (truth=0):\n"))
  cat(sprintf("    RMSE        = %.4f   [should be small]\n", rmse_f12_in))
  cat(sprintf("    band cov 0  = %.3f   [target ~0.95]\n",   cov_f12_in))
  cat(sprintf("    avg band w  = %.3f\n",                     band_width_in))
  cat(sprintf("  Grid f_12 (truth=0, %dx%d):\n", grid_len, grid_len))
  cat(sprintf("    RMSE        = %.4f\n",                     rmse_f12_grid))
  cat(sprintf("    band cov 0  = %.3f   [target ~0.95]\n",   cov_f12_grid))
  cat(sprintf("    avg band w  = %.3f\n",                     band_width_grid))
  cat(sprintf("  Main effects (should be unaffected by null int):\n"))
  cat(sprintf("    RMSE f1     = %.4f\n", rmse_f1))
  cat(sprintf("    RMSE f2     = %.4f\n", rmse_f2))
  cat(sprintf("  Variance components:\n"))
  cat(sprintf("    sigma2 mean = %.3f (true=%.2f)\n", sigma2_mean, true_sigma2))
  cat(sprintf("    tau2   mean = %.3f (true=%.2f)\n", tau2_mean,   true_tau2))
  cat(sprintf("    rho    mean = %.3f (true=%.2f)\n", rho_mean,    true_rho))
  cat(sprintf("    tau2_s_int  : mean=%.4f median=%.4f 95%%CI=(%.4f, %.4f)\n",
              tau2_s_int_mean, tau2_s_int_median,
              tau2_s_int_q025, tau2_s_int_q975))
  cat(sprintf("  MH accept rates: sigma2=%.3f tau2=%.3f rho=%.3f\n",
              gs$accept_sigma2, gs$accept_tau2, gs$accept_rho))

  # ---- 10. Diagnostic verdict ----
  pass_cov_in   <- cov_f12_in   >= 0.90
  pass_cov_grid <- cov_f12_grid >= 0.90
  pass_rmse     <- rmse_f12_grid < 0.30   # generous; tighten if you want
  cat("\n--- Verdict ---\n")
  cat(sprintf("  Coverage in-sample >= 0.90 : %s\n", ifelse(pass_cov_in,   "PASS", "FAIL")))
  cat(sprintf("  Coverage grid      >= 0.90 : %s\n", ifelse(pass_cov_grid, "PASS", "FAIL")))
  cat(sprintf("  RMSE grid          <  0.30 : %s\n", ifelse(pass_rmse,     "PASS", "FAIL")))

  # ---- 11. PDF diagnostics ----
  pdf(outfile, width = 10, height = 8)

  # Traces of variance components
  par(mfrow = c(2, 2))
  plot(gs$sigma2_samples, type = "l", main = "sigma2 trace", ylab = "sigma2")
  abline(h = true_sigma2, col = "red", lty = 2)
  plot(gs$tau2_samples,   type = "l", main = "tau2 trace",   ylab = "tau2")
  abline(h = true_tau2, col = "red", lty = 2)
  plot(gs$rho_samples,    type = "l", main = "rho trace",    ylab = "rho")
  abline(h = true_rho, col = "red", lty = 2)
  plot(gs$tau2_s_int_samp[, 1], type = "l",
       main = "tau2_s_int trace (NULL)", ylab = "tau2_s_int")

  # Posterior density of tau2_s_int with prior overlay
  par(mfrow = c(1, 1))
  ts <- gs$tau2_s_int_samp[, 1]
  d_post <- density(ts)
  # Inverse-Gamma(a=1, b=0.005) prior density (default in sampler)
  a_pr <- 1; b_pr <- 0.005
  xx <- seq(min(ts), max(ts), length.out = 200)
  ig_dens <- (b_pr^a_pr / gamma(a_pr)) * xx^(-a_pr - 1) * exp(-b_pr / xx)
  plot(d_post, main = "tau2_s_int posterior (NULL) vs IG prior",
       xlab = "tau2_s_int", lwd = 2)
  lines(xx, ig_dens, col = "red", lty = 2, lwd = 2)
  legend("topright", c("posterior", "IG(1, 0.005) prior"),
         col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n")

  # In-sample pointwise bands for f_12 (sorted by mean for visual clarity)
  par(mfrow = c(1, 1))
  ord <- order(f12_mean_in)
  plot(seq_along(ord), f12_mean_in[ord], type = "l",
       ylim = range(c(f12_lo_in, f12_hi_in)),
       main = sprintf("In-sample f_12 posterior (truth=0); coverage of 0 = %.2f",
                      cov_f12_in),
       xlab = "obs (sorted by post mean)", ylab = "f_12")
  lines(seq_along(ord), f12_lo_in[ord], lty = 2)
  lines(seq_along(ord), f12_hi_in[ord], lty = 2)
  abline(h = 0, col = "red")

  # Grid surface heatmap of posterior mean
  z_mean <- matrix(f12_mean_grid, grid_len, grid_len)
  image(g, g, z_mean,
        main = sprintf("Posterior mean of f_12 on grid (truth=0)\nrange=[%.3f, %.3f]",
                       min(z_mean), max(z_mean)),
        xlab = "X1", ylab = "X2",
        col = hcl.colors(50, "Blue-Red"))
  contour(g, g, z_mean, add = TRUE)

  # Grid map of band coverage of zero
  cov_grid_mat <- matrix(as.numeric(f12_lo_grid <= 0 & 0 <= f12_hi_grid),
                         grid_len, grid_len)
  image(g, g, cov_grid_mat,
        main = sprintf("Where 95%% CI covers 0 (1=yes); overall = %.2f",
                       cov_f12_grid),
        xlab = "X1", ylab = "X2",
        col = c("tomato", "lightgreen"), zlim = c(0, 1))

  dev.off()
  cat(sprintf("\n  Plots written to: %s\n", outfile))

  invisible(list(
    gs = gs, y = y, X1 = X1, X2 = X2,
    f12_mean_in = f12_mean_in, f12_lo_in = f12_lo_in, f12_hi_in = f12_hi_in,
    f12_mean_grid = f12_mean_grid,
    f12_lo_grid = f12_lo_grid, f12_hi_grid = f12_hi_grid,
    rmse_f12_in = rmse_f12_in, rmse_f12_grid = rmse_f12_grid,
    cov_f12_in = cov_f12_in,   cov_f12_grid = cov_f12_grid,
    rmse_f1 = rmse_f1, rmse_f2 = rmse_f2,
    tau2_s_int_summary = c(mean = tau2_s_int_mean,
                           median = tau2_s_int_median,
                           q025 = tau2_s_int_q025,
                           q975 = tau2_s_int_q975),
    elapsed_sec = as.numeric(t_elapsed)
  ))
}


# ============================================================
# Optional driver: run on Hellbender via SBATCH
# ============================================================
if (!interactive()) {
  source("spatial_utils.R")
  source("ls_basis.R")
  source("ls_interaction.R")
  source("gibbs_interaction.R")

  res <- run_interaction_null_test(
    n      = 300,
    M      = 8,
    n_iter = 3000,
    n_burn = 1000,
    seed   = 42,
    outfile = "interaction_null_test.pdf"
  )
  saveRDS(res, file = "interaction_null_test.rds")
}
