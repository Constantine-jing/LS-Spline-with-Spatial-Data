# ============================================================
# run_sim_poisson_2.R
# Sim_Poisson_2: Low-count vs High-count comparison
#
# Same Sim4 DGP structure (4 real + 2 garbage, spatial GP),
# but we run TWO scenarios to test count-regime sensitivity:
#
#   (A) LOW counts:  intercept mu = 0.8  => mean(y) ~ 3-5
#   (B) HIGH counts: intercept mu = 2.5  => mean(y) ~ 15-30
#
# The nonlinear function scale is adjusted so that both
# scenarios have comparable signal-to-noise on the log scale.
#
# Key question: does the IWLS proposal degrade when counts
# are sparse?  We expect:
#   - Lower b acceptance in the low-count regime
#   - Wider credible bands / higher RMSE for curves
#   - Possibly worse rho/sigma2 recovery
#
# Run: Rscript run_sim_poisson_2.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("gibbs_stage_c_poisson.R")
library(mvtnorm)

# ============================================================
# Settings (shared)
# ============================================================
M_C    <- 20
n      <- 500
n_iter <- 10000
n_burn <- 2500
n_thin <- 1
seed   <- 42
nu     <- 1.5
p      <- 6

x_grid <- seq(0, 1, length.out = 101)

sigma2_true <- 0.3
rho_true    <- 0.2
rho_X       <- 0.10
nu_X        <- 1.0

var_names <- paste0("X", 1:p)

dir.create("poisson_out", showWarnings = FALSE)
stopifnot(ls_tests())


# ============================================================
# Scenario definitions
# ============================================================
scenarios <- list(
  low = list(
    name     = "LOW counts (mean~3-5)",
    tag      = "low",
    mu       = 0.8,     # baseline log(lambda) ~ 0.8 => lambda ~ 2.2
    scale_f  = 0.3,     # tighter function range on log scale
    sigma2   = 0.2      # slightly less spatial variance
  ),
  high = list(
    name     = "HIGH counts (mean~15-30)",
    tag      = "high",
    mu       = 2.5,     # baseline log(lambda) ~ 2.5 => lambda ~ 12
    scale_f  = 0.5,     # same as Sim_Poisson_1
    sigma2   = 0.3      # same as Sim_Poisson_1
  )
)


# ============================================================
# DGP function (parametrized by scenario)
# ============================================================
simulate_poisson_scenario <- function(n, p, scenario, seed = 42) {
  mu      <- scenario$mu
  scale_f <- scenario$scale_f
  sigma2  <- scenario$sigma2

  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D      <- pairdist(coords)

  # Spatially correlated covariates
  R_X     <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(1e-8, n)
  L_X     <- chol(Sigma_X)
  X <- matrix(NA, n, p)
  for (j in 1:p) {
    z <- rnorm(n)
    raw <- as.vector(t(L_X) %*% z)
    X[, j] <- (raw - min(raw)) / (max(raw) - min(raw))
  }
  colnames(X) <- paste0("X", 1:p)

  # Additive predictor
  eta_true <- mu + scale_f * (
    2 * sin(pi * X[,1]) +
    1.5 * exp(X[,2] - 0.5) +
    0.7 * (X[,3]^2) +
    0.5 * sin(2 * pi * X[,4])
  )

  # Spatial random effect
  R_b <- matern_cor(D, rho = rho_true, nu = nu)
  L_b <- chol(R_b + diag(1e-8, n))
  b_true <- as.vector(sqrt(sigma2) * t(L_b) %*% rnorm(n))

  log_lambda <- eta_true + b_true
  y <- rpois(n, lambda = exp(log_lambda))

  # Truth functions for this scenario
  truth_f <- list(
    function(x) scale_f * 2 * sin(pi * x),
    function(x) scale_f * 1.5 * exp(x - 0.5),
    function(x) scale_f * 0.7 * (x^2),
    function(x) scale_f * 0.5 * sin(2 * pi * x),
    function(x) rep(0, length(x)),
    function(x) rep(0, length(x))
  )

  cat(sprintf("  [%s] DGP: mean(y)=%.1f  median(y)=%d  max(y)=%d  zeros=%d/%d\n",
              scenario$tag, mean(y), median(y), max(y), sum(y == 0), n))
  cat(sprintf("  [%s] log(lam) range: [%.2f, %.2f]  mean=%.2f\n",
              scenario$tag, min(log_lambda), max(log_lambda), mean(log_lambda)))

  list(y = y, X = X, coords = coords, D = D,
       eta_true = eta_true, b_true = b_true,
       log_lambda = log_lambda, truth_f = truth_f,
       scenario = scenario)
}


# ============================================================
# Run both scenarios
# ============================================================
all_results <- list()

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  cat(sprintf("\n========== %s ==========\n", sc$name))

  sim <- simulate_poisson_scenario(n, p, sc, seed = seed)

  # Build LS basis
  des <- ls_additive_build(sim$X, M_vec = M_C)
  H   <- cbind(1, des$W)
  col_map <- des$col_map

  # Run sampler
  t0 <- proc.time()

  gs <- gibbs_poisson_sampler(
    y       = sim$y,
    H       = H,
    D       = sim$D,
    nu      = nu,
    col_map = col_map,
    offset  = NULL,
    n_iter  = n_iter,
    n_burn  = n_burn,
    n_thin  = n_thin,
    a_sigma = 2, b_sigma = 0.5,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = -1.6,
    log_rho_sd = 1.0,
    iwls_scale_eta = 1.0,
    iwls_scale_b   = 1.0,
    mh_sd_log_rho  = 0.2,
    verbose = TRUE
  )

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  [%s] Elapsed: %.1f seconds\n", sc$tag, elapsed))

  # --- Summaries ---
  b_post_mean <- colMeans(gs$b_samples)
  b_cor <- cor(b_post_mean, sim$b_true)

  # Curve RMSE
  eta_post_mean <- colMeans(gs$eta_samples)
  rmse_vec <- numeric(p)
  coverage_vec <- numeric(p)

  for (j in 1:p) {
    idx_col <- 1 + col_map[[j]]
    W_grid_j <- des$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)

    # Posterior mean curve
    f_hat_j <- as.vector(W_grid_j %*% eta_post_mean[idx_col])
    f_hat_j <- f_hat_j - mean(f_hat_j)
    f_true_j <- sim$truth_f[[j]](x_grid)
    f_true_j <- f_true_j - mean(f_true_j)
    rmse_vec[j] <- sqrt(mean((f_hat_j - f_true_j)^2))

    # 95% band coverage
    f_samples <- gs$eta_samples[, idx_col, drop = FALSE] %*% t(W_grid_j)
    f_samples <- f_samples - rowMeans(f_samples)
    f_lo <- apply(f_samples, 2, quantile, 0.025)
    f_hi <- apply(f_samples, 2, quantile, 0.975)
    coverage_vec[j] <- mean(f_true_j >= f_lo & f_true_j <= f_hi)
  }

  result <- list(
    gs = gs, sim = sim, des = des, col_map = col_map,
    elapsed = elapsed,
    sigma2_post = mean(gs$sigma2_samples),
    rho_post = mean(gs$rho_samples),
    b_cor = b_cor,
    rmse = rmse_vec,
    coverage = coverage_vec,
    tau2_s_means = colMeans(gs$tau2_s_samples),
    scenario = sc
  )
  all_results[[sc_name]] <- result

  cat(sprintf("  [%s] sigma2: post=%.4f true=%.4f\n", sc$tag,
              result$sigma2_post, sc$sigma2))
  cat(sprintf("  [%s] rho:    post=%.4f true=%.4f\n", sc$tag,
              result$rho_post, rho_true))
  cat(sprintf("  [%s] b_cor:  %.4f\n", sc$tag, b_cor))
  cat(sprintf("  [%s] RMSE:   %s\n", sc$tag,
              paste(sprintf("%.5f", rmse_vec), collapse = "  ")))
  cat(sprintf("  [%s] Cover:  %s\n", sc$tag,
              paste(sprintf("%.3f", coverage_vec), collapse = "  ")))
}


# ============================================================
# Comparison summary table
# ============================================================
cat("\n\n#########################################\n")
cat("# Sim_Poisson_2: Low vs High Count Comparison\n")
cat("#########################################\n\n")

cat(sprintf("%-12s  %8s  %8s  %8s  %8s  %8s\n",
            "Scenario", "sigma2", "rho", "b_cor", "acc_eta", "acc_b"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (sc_name in names(all_results)) {
  r <- all_results[[sc_name]]
  cat(sprintf("%-12s  %8.4f  %8.4f  %8.4f  %8.3f  %8.3f\n",
              r$scenario$tag,
              r$sigma2_post, r$rho_post, r$b_cor,
              r$gs$accept_rate["eta"], r$gs$accept_rate["b"]))
}
cat(sprintf("%-12s  %8.4f  %8.4f\n", "TRUE(low)",
            scenarios$low$sigma2, rho_true))
cat(sprintf("%-12s  %8.4f  %8.4f\n", "TRUE(high)",
            scenarios$high$sigma2, rho_true))

cat("\n--- Curve RMSE: X1-X4 ---\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %8s\n", "Scenario", "X1", "X2", "X3", "X4"))
cat(paste(rep("-", 52), collapse = ""), "\n")
for (sc_name in names(all_results)) {
  r <- all_results[[sc_name]]
  cat(sprintf("%-12s  %8.5f  %8.5f  %8.5f  %8.5f\n",
              r$scenario$tag, r$rmse[1], r$rmse[2], r$rmse[3], r$rmse[4]))
}

cat("\n--- 95% Band Coverage: X1-X6 ---\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %8s  %8s  %8s\n",
            "Scenario", "X1", "X2", "X3", "X4", "X5", "X6"))
cat(paste(rep("-", 64), collapse = ""), "\n")
for (sc_name in names(all_results)) {
  r <- all_results[[sc_name]]
  cat(sprintf("%-12s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
              r$scenario$tag,
              r$coverage[1], r$coverage[2], r$coverage[3],
              r$coverage[4], r$coverage[5], r$coverage[6]))
}

cat("\n--- tau2_s shrinkage (garbage X5, X6) ---\n")
for (sc_name in names(all_results)) {
  r <- all_results[[sc_name]]
  cat(sprintf("  [%s] real(X1-X4): %s\n", r$scenario$tag,
              paste(sprintf("%.4f", r$tau2_s_means[1:4]), collapse = " ")))
  cat(sprintf("  [%s] garb(X5-X6): %s\n", r$scenario$tag,
              paste(sprintf("%.4f", r$tau2_s_means[5:6]), collapse = " ")))
}

cat("\n--- Timing ---\n")
for (sc_name in names(all_results)) {
  r <- all_results[[sc_name]]
  cat(sprintf("  [%s] %.0f seconds\n", r$scenario$tag, r$elapsed))
}


# ============================================================
# Plots: side-by-side marginal curves
# ============================================================
pdf("poisson_out/sim_poisson_2_comparison.pdf", width = 14, height = 10)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

for (sc_name in c("low", "high")) {
  r <- all_results[[sc_name]]
  for (j in 1:4) {
    idx_col <- 1 + r$col_map[[j]]
    W_grid_j <- r$des$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)

    f_samples <- r$gs$eta_samples[, idx_col, drop = FALSE] %*% t(W_grid_j)
    f_samples <- f_samples - rowMeans(f_samples)
    f_mean <- colMeans(f_samples)
    f_lo   <- apply(f_samples, 2, quantile, 0.025)
    f_hi   <- apply(f_samples, 2, quantile, 0.975)

    f_true <- r$sim$truth_f[[j]](x_grid)
    f_true <- f_true - mean(f_true)

    ylim <- range(c(f_lo, f_hi, f_true))

    plot(x_grid, f_true, type = "l", lty = 2, lwd = 2, col = "black",
         ylim = ylim, xlab = var_names[j],
         ylab = paste0("f(", var_names[j], ")"),
         main = sprintf("%s — %s (RMSE=%.4f)", toupper(sc_name),
                        var_names[j], r$rmse[j]))
    polygon(c(x_grid, rev(x_grid)),
            c(f_lo, rev(f_hi)),
            col = rgb(0.3, 0.5, 0.9, 0.3), border = NA)
    lines(x_grid, f_mean, col = "blue", lwd = 2)
  }
}

dev.off()
cat("\n  Saved: poisson_out/sim_poisson_2_comparison.pdf\n")


# ============================================================
# Plots: trace plots for both scenarios
# ============================================================
pdf("poisson_out/sim_poisson_2_traces.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

for (sc_name in c("low", "high")) {
  r <- all_results[[sc_name]]

  plot(r$gs$sigma2_samples, type = "l", col = "steelblue",
       main = sprintf("%s: sigma2 (true=%.2f)", toupper(sc_name),
                      r$scenario$sigma2),
       ylab = expression(sigma^2), xlab = "iteration")
  abline(h = r$scenario$sigma2, col = "red", lty = 2, lwd = 2)

  plot(r$gs$rho_samples, type = "l", col = "forestgreen",
       main = sprintf("%s: rho (true=%.2f)", toupper(sc_name), rho_true),
       ylab = expression(rho), xlab = "iteration")
  abline(h = rho_true, col = "red", lty = 2, lwd = 2)
}

dev.off()
cat("  Saved: poisson_out/sim_poisson_2_traces.pdf\n")


# ============================================================
# Plots: b recovery comparison
# ============================================================
pdf("poisson_out/sim_poisson_2_spatial_b.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

for (sc_name in c("low", "high")) {
  r <- all_results[[sc_name]]
  b_pm <- colMeans(r$gs$b_samples)
  plot(r$sim$b_true, b_pm,
       xlab = "True b", ylab = "Posterior mean b",
       main = sprintf("%s counts (cor=%.3f)", toupper(sc_name), r$b_cor),
       pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.4))
  abline(0, 1, col = "red", lty = 2, lwd = 2)
}

dev.off()
cat("  Saved: poisson_out/sim_poisson_2_spatial_b.pdf\n")

cat("\n=== Done ===\n")
