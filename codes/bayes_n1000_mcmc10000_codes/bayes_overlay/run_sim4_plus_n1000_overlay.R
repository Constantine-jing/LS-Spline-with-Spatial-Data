# ============================================================
# run_sim4_plus_n1000.R
# Sim4_plus: 4 real covariates (X1-X4) + 6 garbage covariates (X5-X10)
#   n = 1000, M = 20 knots, p = 10
#   Full Bayes (rho sampled, RW2), 10000 MCMC iterations
#   vs REML (M=6)
#
# Stronger test than Sim4 (p=6): with 6 garbage covariates the RW2
# smoothing shrinkage must work harder to flatten X5-X10.
# Matches the spirit of Nandy et al. (2017) who test up to 35 covariates
# with only 4 real signals.
#
# MH proposals (tuned):
#   mh_sd_log_sigma2 = 0.5
#   mh_sd_log_tau2   = 0.3
#   mh_sd_log_rho    = 0.3
#
# Run: Rscript run_sim4_plus_n1000.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

# ============================================================
# Settings
# ============================================================
M_C    <- 20
M_reml <- 6
n      <- 1000
n_iter <- 10000
n_burn <- 2500        # 25% burn-in
n_thin <- 1
seed   <- 42
nu     <- 1.5
p      <- 10          # 4 real + 6 garbage

x_grid <- seq(0, 1, length.out = 101)

sigma2_true <- 0.8
rho_true    <- 0.2
tau2_true   <- 0.15
rho_X       <- 0.10
nu_X        <- 1.0

mh_sd_sigma2 <- 0.5
mh_sd_tau2   <- 0.3
mh_sd_rho    <- 0.3

var_names <- paste0("X", 1:p)

# Truth: X1-X4 real signals, X5-X10 identically zero
truth_f_list_10 <- c(
  default_truth_f_list(),                            # f1..f4
  replicate(6, function(x) rep(0, length(x)),
            simplify = FALSE)                         # f5..f10
)

dir.create("bayes_test", showWarnings = FALSE)
stopifnot(ls_tests())

# ============================================================
# Simulation: Sim4_plus
#   X1-X4: spatially correlated (Matern), rank-uniform-transformed
#   X5-X10: iid Uniform(0,1) -- garbage predictors
#   b:      Matern spatial random effect
#   eta:    depends on X1-X4 only; f5=...=f10=0
# ============================================================
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim4_plus <- function(n = 1000, p = 10, mu = 0,
                               rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                               sigma2 = 0.8, rho = 0.2, nu = 1.5,
                               tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords  <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D       <- pairdist(coords)

  # X1-X4: spatially correlated, rank-uniform-transformed
  R_X     <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })

  # X5-X10: iid Uniform(0,1) garbage
  X5_10 <- matrix(runif(n * 6), nrow = n, ncol = 6)

  X <- cbind(X14, X5_10)
  colnames(X) <- paste0("X", 1:p)

  eta <- eta_truth_additive(X[, 1:4, drop = FALSE], mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b   <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))
  eps <- rnorm(n, 0, sqrt(tau2))

  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = eta + b + eps, eta = eta, b = b, eps = eps, X)
}

# ============================================================
# Generate data
# ============================================================
cat(sprintf("\n=== Sim4_plus: n=%d, p=%d (4 real + 6 garbage) ===\n", n, p))
dat    <- simulate_sim4_plus(n = n, p = p, rho_X = rho_X, nu_X = nu_X,
                             sigma2 = sigma2_true, rho = rho_true, nu = nu,
                             tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat[, c("x", "y_coord")])
X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
y      <- dat$Y
b_true <- dat$b
D      <- pairdist(coords)
cat(sprintf("  Y range: [%.3f, %.3f]  mean=%.3f\n", min(y), max(y), mean(y)))

# ============================================================
# REML (M=6)
# ============================================================
cat("\n####################################################\n")
cat("# REML: M=6\n")
cat("####################################################\n")

t0_reml  <- proc.time()
reml_obj <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                           M_vec = rep(M_reml, p), nu = nu,
                           rho_init    = rho_true,
                           lambda_init = tau2_true / sigma2_true,
                           verbose = TRUE)
t_reml   <- (proc.time() - t0_reml)[3]
reml_obj$X_raw_for_marginal <- X_raw

cat(sprintf("  REML done in %.1fs\n", t_reml))
cat(sprintf("  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho))

reml_curves <- marginal_curves(reml_obj, x_grid = x_grid,
                               truth_f_list = truth_f_list_10, clip = TRUE)
reml_err    <- curve_error_table(reml_curves, var_names = var_names)
reml_rmse   <- reml_err$rmse_curve
cat("  REML curve RMSE:\n")
print(data.frame(var = var_names, rmse = round(reml_rmse, 5)))

# REML marginal plot (5 rows x 2 cols for 10 covariates)
# NOTE: Marginal plots deferred until after Bayes for overlay.

# ============================================================
# Full Bayes (M=20)
# ============================================================
cat(sprintf("\n####################################################\n"))
cat(sprintf("# Full Bayes: M=%d  n_iter=%d  n_burn=%d\n", M_C, n_iter, n_burn))
cat(sprintf("# MH sd: sigma2=%.1f  tau2=%.1f  rho=%.1f\n",
            mh_sd_sigma2, mh_sd_tau2, mh_sd_rho))
cat("####################################################\n")

des <- ls_additive_build(X_raw, M_vec = rep(M_C, p))
H   <- cbind(1, des$W)
cat(sprintf("  H: %d x %d\n", nrow(H), ncol(H)))

# REML initialisation
obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                       M_vec = rep(6, p), nu = nu,
                       rho_init    = rho_true,
                       lambda_init = tau2_true / max(sigma2_true, 0.01),
                       verbose = FALSE)
R_init     <- matern_cor(D, rho = obj6$fit$rho, nu = nu)
Sigma_init <- max(obj6$fit$sigma2, 0.01) * R_init + max(obj6$fit$tau2, 0.01) * diag(n)
L_init     <- chol(Sigma_init + diag(1e-8, n))
y_w        <- forwardsolve(t(L_init), y)
H_w        <- forwardsolve(t(L_init), H)
eta_init   <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H)),
                               crossprod(H_w, y_w)))
init <- list(eta    = eta_init,
             sigma2 = max(obj6$fit$sigma2, 0.01),
             tau2   = max(obj6$fit$tau2,   0.01),
             rho    = obj6$fit$rho,
             tau2_s = rep(1.0, p))
cat(sprintf("  Init: rho=%.4f  sigma2=%.4f  tau2=%.4f\n",
            init$rho, init$sigma2, init$tau2))

# Gibbs sampler
t0_bayes <- proc.time()
gs <- gibbs_full_sampler(
  y = y, H = H, D = D, nu = nu, col_map = des$col_map,
  n_iter = n_iter, n_burn = n_burn, n_thin = n_thin,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 0.3,
  a_smooth = 1, b_smooth = 0.005,
  log_rho_mu = log(0.2), log_rho_sd = 1.0,
  mh_sd_log_sigma2 = mh_sd_sigma2,
  mh_sd_log_tau2   = mh_sd_tau2,
  mh_sd_log_rho    = mh_sd_rho,
  init = init, verbose = TRUE
)
t_bayes <- (proc.time() - t0_bayes)[3]
cat(sprintf("  Full Bayes done in %.1fs\n", t_bayes))

gsm         <- gibbs_summary(gs)
rho_summary <- c(mean = mean(gs$rho_samples), sd = sd(gs$rho_samples),
                 quantile(gs$rho_samples, c(0.025, 0.5, 0.975)))
tau2s_means <- colMeans(gs$tau2_s_samples)
names(tau2s_means) <- var_names

cat(sprintf("  sigma2: mean=%.4f  95%%CI=[%.4f,%.4f]  true=%.2f\n",
            gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"], sigma2_true))
cat(sprintf("  tau2:   mean=%.4f  95%%CI=[%.4f,%.4f]  true=%.2f\n",
            gsm$tau2["mean"],   gsm$tau2["2.5%"],   gsm$tau2["97.5%"],   tau2_true))
cat(sprintf("  rho:    mean=%.4f  95%%CI=[%.4f,%.4f]  true=%.2f\n",
            rho_summary["mean"], rho_summary["2.5%"], rho_summary["97.5%"], rho_true))
cat(sprintf("  MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))

cat("  tau2_s (key: X5-X10 should be near 0):\n")
print(round(tau2s_means, 5))

b_post <- colMeans(gs$b_samples)
cat(sprintf("  b: cor=%.4f  RMSE=%.4f\n",
            cor(b_true, b_post), sqrt(mean((b_true - b_post)^2))))

obj_bayes <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)

cov_tab <- gibbs_curve_coverage(gs, obj_bayes, x_grid, truth_f_list_10,
                                X_raw = X_raw, var_names = var_names)
cat("  Curve metrics:\n"); print(cov_tab, digits = 4)

# ============================================================
# Overlay plot: REML + Bayes on the same panel per covariate
# ============================================================
bayes_bands <- vector("list", p)
for (j in 1:p) {
  tf <- truth_f_list_10[[j]]
  bayes_bands[[j]] <- gibbs_marginal_band(gs, obj_bayes, j, x_grid,
                                           X_raw = X_raw, truth_f = tf,
                                           level = 0.95)
}

all_y_reml  <- c(unlist(reml_curves$fhat_grid), unlist(reml_curves$ftrue_grid))
all_y_bayes <- unlist(lapply(bayes_bands, function(b)
                 c(b$lower, b$upper, b$f_hat, b$f_true)))
all_y <- c(all_y_reml, all_y_bayes)
all_y <- all_y[is.finite(all_y)]
shared_rng <- range(all_y)
a <- max(abs(shared_rng)); shared_rng <- c(-a, a)
span <- diff(shared_rng); if (span == 0) span <- 1
shared_ylim <- shared_rng + c(-1, 1) * 0.05 * span
cat(sprintf("  Shared ylim: [%.3f, %.3f]\n", shared_ylim[1], shared_ylim[2]))

pdf("bayes_test/sim4plus_n1000_overlay_marginals.pdf", width = 12, height = 15)
par(mfrow = c(5, 2), mar = c(4, 4, 3, 1))
for (j in 1:p) {
  b <- bayes_bands[[j]]
  reml_fhat <- reml_curves$fhat_grid[[j]]
  plot(b$x, b$f_hat, type = "n",
       main = paste0("Marginal: ", var_names[j]),
       xlab = "x", ylab = paste0("f_", j, "(x)"),
       ylim = shared_ylim)
  polygon(c(b$x, rev(b$x)), c(b$lower, rev(b$upper)),
          col = rgb(0.2, 0.4, 0.8, 0.15), border = NA)
  lines(x_grid, reml_fhat, lwd = 2, col = "darkorange")
  lines(b$x, b$f_hat, lwd = 2, col = "blue")
  if (!all(is.na(b$f_true)))
    lines(b$x, b$f_true, lty = 2, col = "black", lwd = 1.5)
  abline(h = 0, lty = 3, col = "gray50")
  r_reml  <- sqrt(mean((reml_fhat - b$f_true)^2))
  r_bayes <- sqrt(mean((b$f_hat   - b$f_true)^2))
  inside  <- b$f_true >= b$lower & b$f_true <= b$upper
  cov_pct <- mean(inside) * 100
  mtext(sprintf("RMSE: REML=%.4f  Bayes=%.4f  cov=%.0f%%",
                r_reml, r_bayes, cov_pct),
        side = 3, line = 0.2, cex = 0.7)
  legend("topleft",
         legend = c("truth", "REML (M=6)", "Bayes (M=20)", "95% band"),
         col    = c("black", "darkorange", "blue", rgb(0.2, 0.4, 0.8, 0.3)),
         lty    = c(2, 1, 1, NA),
         lwd    = c(1.5, 2, 2, 10),
         bty = "n", cex = 0.75)
}
dev.off()
cat("  Plot: bayes_test/sim4plus_n1000_overlay_marginals.pdf\n")

# Trace plot
pdf("bayes_test/sim4plus_n1000_bayes_trace.pdf", width = 10, height = 10)
plot_gibbs_trace_full(gs, true_sigma2 = sigma2_true,
                          true_tau2   = tau2_true,
                          true_rho    = rho_true)
dev.off()
cat("  Plot: bayes_test/sim4plus_n1000_bayes_trace.pdf\n")

# ============================================================
# Summary
# ============================================================
cat("\n\n####################################################\n")
cat("# REML vs Bayes: Sim4_plus n=1000\n")
cat("####################################################\n\n")

cat("--- Curve RMSE ---\n")
cat(sprintf("%-10s", ""))
for (v in var_names) cat(sprintf("  %7s", v)); cat("\n")
cat(sprintf("%-10s", "REML"))
for (r in reml_rmse) cat(sprintf("  %7.5f", r)); cat("\n")
cat(sprintf("%-10s", "Bayes"))
for (r in cov_tab$rmse_curve) cat(sprintf("  %7.5f", r)); cat("\n")

cat("\n--- tau2_s: real (X1-X4) vs garbage (X5-X10) ---\n")
cat("  Real:   "); cat(sprintf("%.5f  ", tau2s_means[1:4])); cat("\n")
cat("  Garbage:"); cat(sprintf("%.5f  ", tau2s_means[5:10])); cat("\n")

cat("\n--- Variance estimates ---\n")
cat(sprintf("  REML:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho))
cat(sprintf("  Bayes: sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            gsm$sigma2["mean"], gsm$tau2["mean"], rho_summary["mean"]))
cat(sprintf("  True:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            sigma2_true, tau2_true, rho_true))

cat(sprintf("\n--- Timing: REML=%.1fs  Bayes=%.1fs ---\n", t_reml, t_bayes))

# Save CSVs
rmse_row_reml  <- setNames(as.list(reml_rmse),  paste0("rmse_", var_names))
rmse_row_bayes <- setNames(as.list(cov_tab$rmse_curve), paste0("rmse_", var_names))
comp_table <- data.frame(
  method   = c("REML", "Bayes"),
  do.call(rbind, lapply(list(rmse_row_reml, rmse_row_bayes), as.data.frame)),
  sigma2   = c(reml_obj$fit$sigma2,  gsm$sigma2["mean"]),
  tau2     = c(reml_obj$fit$tau2,    gsm$tau2["mean"]),
  rho      = c(reml_obj$fit$rho,     rho_summary["mean"]),
  time_sec = c(t_reml, t_bayes),
  row.names = NULL
)
write.csv(comp_table, "bayes_test/sim4plus_n1000_reml_vs_bayes.csv", row.names = FALSE)
write.csv(cov_tab,   "bayes_test/sim4plus_n1000_bayes_coverage.csv",  row.names = FALSE)

cat("\nOutputs in bayes_test/:\n")
cat("  sim4plus_n1000_overlay_marginals.pdf\n")
cat("  sim4plus_n1000_bayes_trace.pdf\n")
cat("  sim4plus_n1000_reml_vs_bayes.csv\n")
cat("  sim4plus_n1000_bayes_coverage.csv\n")
