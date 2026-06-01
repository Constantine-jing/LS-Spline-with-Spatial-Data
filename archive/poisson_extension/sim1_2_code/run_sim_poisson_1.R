# ============================================================
# run_sim_poisson_1.R
# Sim_Poisson_1: Poisson extension of Sim4
#   - 4 real covariates (X1-X4) + 2 garbage (X5-X6)
#   - Poisson(lambda_i), log(lambda_i) = eta_i + b_i
#   - Spatial GP on b, Matern covariance
#   - n = 500 (Poisson is more expensive due to no collapse)
#   - M_C = 20 knots for Bayes
#
# True functions are SCALED DOWN from the Gaussian sims so that
# log(lambda) stays in a reasonable range.  With the Gaussian
# sims, eta ranged roughly from -2 to 4 (before b).  For
# Poisson, we want log(lambda) in [0.5, 4] => counts ~2 to 55.
# We use a global scale factor of 0.5 and an intercept of 1.5.
#
# Run: Rscript run_sim_poisson_1.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("gibbs_stage_c_poisson.R")
library(mvtnorm)

# ============================================================
# Settings
# ============================================================
M_C    <- 20
n      <- 500
n_iter <- 10000
n_burn <- 2500
n_thin <- 1
seed   <- 42
nu     <- 1.5
p      <- 6           # 4 real + 2 garbage

x_grid <- seq(0, 1, length.out = 101)

sigma2_true <- 0.3     # spatial variance (on log scale)
rho_true    <- 0.2     # Matern range
rho_X       <- 0.10    # covariate spatial correlation
nu_X        <- 1.0
mu_true     <- 1.5     # intercept => baseline log(lambda) ~ 1.5

# Scale factor: shrink the Gaussian truth functions so
# the total log(lambda) is in a reasonable count range
scale_f <- 0.5

var_names <- paste0("X", 1:p)

# Truth functions (scaled)
truth_f_list_poisson <- list(
  function(x) scale_f * 2 * sin(pi * x),           # max ~ 1.0
  function(x) scale_f * 1.5 * exp(x - 0.5),        # range ~ [0.45, 1.22]
  function(x) scale_f * 0.7 * (x^2),               # max ~ 0.35
  function(x) scale_f * 0.5 * sin(2 * pi * x),     # max ~ 0.25
  function(x) rep(0, length(x)),                    # garbage
  function(x) rep(0, length(x))                     # garbage
)

dir.create("poisson_out", showWarnings = FALSE)
stopifnot(ls_tests())


# ============================================================
# Simulation DGP
# ============================================================
eta_truth_poisson <- function(X, mu = mu_true) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + scale_f * (2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) +
                    0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4]))
}

simulate_poisson_1 <- function(n = 500, p = 6, mu = mu_true,
                                rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                                sigma2 = 0.3, rho = 0.2, nu = 1.5,
                                seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords  <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D       <- pairdist(coords)

  # Spatially correlated covariates
  R_X     <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  L_X     <- chol(Sigma_X)
  X <- matrix(NA, n, p)
  for (j in 1:p) {
    z <- rnorm(n)
    raw <- as.vector(t(L_X) %*% z)
    X[, j] <- (raw - min(raw)) / (max(raw) - min(raw))
  }
  colnames(X) <- paste0("X", 1:p)

  # Additive predictor (on log scale)
  eta_true <- eta_truth_poisson(X, mu)

  # Spatial random effect b ~ N(0, sigma2 * R)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  L_b <- chol(R_b + diag(jitter_b, n))
  b_true <- as.vector(sqrt(sigma2) * t(L_b) %*% rnorm(n))

  # Log-lambda and Poisson draw
  log_lambda <- eta_true + b_true
  y <- rpois(n, lambda = exp(log_lambda))

  cat(sprintf("  DGP summary: mean(y)=%.2f  median(y)=%.0f  max(y)=%d  zeros=%d/%d\n",
              mean(y), median(y), max(y), sum(y == 0), n))
  cat(sprintf("  log(lambda) range: [%.2f, %.2f],  mean=%.2f\n",
              min(log_lambda), max(log_lambda), mean(log_lambda)))

  list(y = y, X = X, coords = coords, D = D,
       eta_true = eta_true, b_true = b_true,
       log_lambda = log_lambda)
}


# ============================================================
# Run simulation
# ============================================================
cat("=== Sim_Poisson_1: Poisson spatial additive model ===\n")
cat(sprintf("  n=%d, M_C=%d, n_iter=%d\n", n, M_C, n_iter))

sim <- simulate_poisson_1(n = n, p = p, sigma2 = sigma2_true,
                           rho = rho_true, seed = seed)

# Build LS basis
des <- ls_additive_build(sim$X, M_vec = M_C)
H   <- cbind(1, des$W)
col_map <- des$col_map

cat("  H dimensions:", dim(H), "\n")

# ============================================================
# Run Poisson Gibbs sampler
# ============================================================
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
cat(sprintf("  Elapsed: %.1f seconds\n", elapsed))


# ============================================================
# Posterior summaries
# ============================================================
cat("\n--- Variance recovery ---\n")
cat(sprintf("  sigma2: posterior mean=%.4f, true=%.4f\n",
            mean(gs$sigma2_samples), sigma2_true))
cat(sprintf("  rho:    posterior mean=%.4f, true=%.4f\n",
            mean(gs$rho_samples), rho_true))

# Spatial random effect correlation
b_post_mean <- colMeans(gs$b_samples)
b_cor <- cor(b_post_mean, sim$b_true)
cat(sprintf("  b correlation (posterior mean vs true): %.4f\n", b_cor))

# tau2_s recovery
tau2_s_means <- colMeans(gs$tau2_s_samples)
cat("  tau2_s posterior means:", paste(sprintf("%.4f", tau2_s_means), collapse = ", "), "\n")


# ============================================================
# Marginal curve evaluation on x_grid
# ============================================================
cat("\n--- Curve RMSE (posterior mean vs truth) ---\n")
eta_post_mean <- colMeans(gs$eta_samples)

for (j in 1:p) {
  idx_col <- 1 + col_map[[j]]
  beta_j_post <- eta_post_mean[idx_col]

  # Evaluate on x_grid
  W_grid_j <- des$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)
  f_hat_j  <- as.vector(W_grid_j %*% beta_j_post)

  # Center (remove mean for comparison)
  f_hat_j <- f_hat_j - mean(f_hat_j)

  # Truth
  f_true_j <- truth_f_list_poisson[[j]](x_grid)
  f_true_j <- f_true_j - mean(f_true_j)

  rmse_j <- sqrt(mean((f_hat_j - f_true_j)^2))
  cat(sprintf("  %s: RMSE=%.5f\n", var_names[j], rmse_j))
}


# ============================================================
# Plots: trace plots
# ============================================================
pdf("poisson_out/sim_poisson_1_traces.pdf", width = 10, height = 6)
plot_gibbs_trace_poisson(gs, true_sigma2 = sigma2_true, true_rho = rho_true)
dev.off()
cat("  Saved: poisson_out/sim_poisson_1_traces.pdf\n")


# ============================================================
# Plots: marginal curves with 95% credible bands
# ============================================================
pdf("poisson_out/sim_poisson_1_marginals.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

for (j in 1:p) {
  idx_col <- 1 + col_map[[j]]

  # Posterior samples of f_j on grid
  W_grid_j <- des$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)
  f_samples <- gs$eta_samples[, idx_col, drop = FALSE] %*% t(W_grid_j)

  # Center each sample
  f_samples <- f_samples - rowMeans(f_samples)

  f_mean <- colMeans(f_samples)
  f_lo   <- apply(f_samples, 2, quantile, 0.025)
  f_hi   <- apply(f_samples, 2, quantile, 0.975)

  # Truth
  f_true <- truth_f_list_poisson[[j]](x_grid)
  f_true <- f_true - mean(f_true)

  ylim <- range(c(f_lo, f_hi, f_true))

  plot(x_grid, f_true, type = "l", lty = 2, lwd = 2, col = "black",
       ylim = ylim, xlab = var_names[j], ylab = paste0("f(", var_names[j], ")"),
       main = paste0(var_names[j], " (Poisson)"))
  polygon(c(x_grid, rev(x_grid)),
          c(f_lo, rev(f_hi)),
          col = rgb(0.3, 0.5, 0.9, 0.3), border = NA)
  lines(x_grid, f_mean, col = "blue", lwd = 2)
  legend("topright", legend = c("Truth", "Post mean", "95% CI"),
         col = c("black", "blue", rgb(0.3, 0.5, 0.9, 0.5)),
         lty = c(2, 1, NA), lwd = c(2, 2, NA),
         pch = c(NA, NA, 15), pt.cex = 2, cex = 0.7)
}

dev.off()
cat("  Saved: poisson_out/sim_poisson_1_marginals.pdf\n")


# ============================================================
# Plot: spatial random effect recovery
# ============================================================
pdf("poisson_out/sim_poisson_1_spatial_b.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# True vs estimated b
plot(sim$b_true, b_post_mean,
     xlab = "True b", ylab = "Posterior mean b",
     main = sprintf("Spatial random effect (cor=%.3f)", b_cor),
     pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.4))
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Spatial map of b
plot(sim$coords[,1], sim$coords[,2],
     col = hcl.colors(100, "RdYlBu")[
       cut(b_post_mean, breaks = 100, labels = FALSE)
     ],
     pch = 16, cex = 1.2,
     xlab = "x", ylab = "y",
     main = "Posterior mean b (spatial)")

dev.off()
cat("  Saved: poisson_out/sim_poisson_1_spatial_b.pdf\n")

cat("\n=== Done ===\n")
