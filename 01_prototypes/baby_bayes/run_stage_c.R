# ============================================================
# run_stage_c.R
# Stage C: RW2 smoothing prior + more knots (M=20)
#
# Tests whether:
#   1) More knots + smoothing prior improves X1 coverage
#      (which was stuck at 11% due to bias with M=6)
#   2) Smoothing variances tau2_s_j adapt per covariate
#   3) sigma2 and tau2 still recovered
#
# Run: Rscript run_stage_c.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")       # for plotting/summary utilities
source("gibbs_bayes_v2.R")    # for collapsed sampler reference
source("gibbs_stage_c.R")
library(mvtnorm)

M_C    <- 20    # more knots for Stage C (was 6)
M_old  <- 6     # for Stage B comparison
seed   <- 42
nu     <- 1.5
x_grid <- seq(0, 1, length.out = 101)

truth_f_list <- default_truth_f_list()
var_names    <- paste0("X", 1:4)

dir.create("bayes_out", showWarnings = FALSE)

# ---- Simulation functions ----
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim2 <- function(n = 400, p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  X <- matrix(runif(n*p,0,1), n, p); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu = mu)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2*R + diag(jitter,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

stopifnot(ls_tests())

# ============================================================
# Sim2, n=400, M=20 knots + RW2
# ============================================================
cat("\n####################################################\n")
cat("# Stage C: RW2 smoothing + M=20 knots, Sim2 n=400\n")
cat("####################################################\n")

sigma2_true <- 0.8; rho_true <- 0.2; tau2_true <- 0.15; n <- 400

dat <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                     nu = nu, tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat[, c("x","y_coord")])
X_raw  <- as.matrix(dat[, paste0("X",1:4)])
y      <- dat$Y
b_true <- dat$b

# --- Build LS basis with M=20 knots ---
des_C <- ls_additive_build(X_raw, M_vec = rep(M_C, 4))
W_C   <- des_C$W
H_C   <- cbind(1, W_C)  # intercept + identified spline blocks

cat(sprintf("Design matrix H: %d x %d (M=%d, p=%d coefs per covariate)\n",
            nrow(H_C), ncol(H_C), M_C, M_C - 1))

# --- Also build M=6 version for REML init + comparison ---
obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                       M_vec = rep(M_old, 4), nu = nu,
                       rho_init = rho_true,
                       lambda_init = tau2_true / sigma2_true,
                       verbose = TRUE)
obj6$X_raw_for_marginal <- X_raw

# --- Build R (rho fixed at truth) ---
D <- pairdist(coords)
R <- matern_cor(D, rho = rho_true, nu = nu)

# --- Initialize ---
# Use OLS on the M=20 design as starting values for eta
# (crude but gets us in the right ballpark)
Sigma_init <- obj6$fit$sigma2 * R + obj6$fit$tau2 * diag(n)
L_init <- chol(Sigma_init + diag(1e-8, n))
y_w <- forwardsolve(t(L_init), y)
H_w <- forwardsolve(t(L_init), H_C)
eta_init <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H_C)),
                             crossprod(H_w, y_w)))

init_C <- list(
  eta    = eta_init,
  sigma2 = max(obj6$fit$sigma2, 0.01),
  tau2   = max(obj6$fit$tau2, 0.01),
  tau2_s = rep(1.0, 4)  # start smoothing variances at 1
)

# --- Run Stage C Gibbs ---
cat("\nRunning Stage C Gibbs (M=20, RW2 prior)...\n")
t0 <- proc.time()

gs_C <- gibbs_stage_c_sampler(
  y = y, H = H_C, R = R, col_map = des_C$col_map,
  n_iter = 6000, n_burn = 1500, n_thin = 1,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 0.3,
  a_smooth = 1, b_smooth = 0.005,
  mh_sd_log_sigma2 = 0.3,
  mh_sd_log_tau2   = 0.3,
  init = init_C, verbose = TRUE
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("Stage C done in %.1f seconds.\n", elapsed))

# --- Summaries ---
gsm_C <- gibbs_summary(gs_C)

cat("\n--- Variance posteriors (Stage C) ---\n")
cat(sprintf("  sigma2: mean=%.4f  sd=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm_C$sigma2["mean"], gsm_C$sigma2["sd"],
            gsm_C$sigma2["2.5%"], gsm_C$sigma2["97.5%"], sigma2_true))
cat(sprintf("  tau2:   mean=%.4f  sd=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm_C$tau2["mean"], gsm_C$tau2["sd"],
            gsm_C$tau2["2.5%"], gsm_C$tau2["97.5%"], tau2_true))

cat("\n--- Smoothing variances tau2_s (per covariate) ---\n")
for (j in 1:4) {
  cat(sprintf("  %s: mean=%.4f  median=%.4f  95%%CI=[%.4f, %.4f]\n",
              var_names[j],
              mean(gs_C$tau2_s_samples[, j]),
              median(gs_C$tau2_s_samples[, j]),
              quantile(gs_C$tau2_s_samples[, j], 0.025),
              quantile(gs_C$tau2_s_samples[, j], 0.975)))
}

# --- Build an obj-like structure for the M=20 basis (for band functions) ---
obj_C <- list(
  des   = des_C,
  X_fix = H_C,
  fit   = list(beta = colMeans(gs_C$eta_samples)),
  X_raw_for_marginal = X_raw
)

# --- Curve coverage (Stage C) ---
cov_C <- gibbs_curve_coverage(gs_C, obj_C, x_grid, truth_f_list,
                               X_raw = X_raw, var_names = var_names)

cat("\n--- Stage C curve metrics ---\n")
print(cov_C, digits = 4)

# --- Stage B (v2, M=6) comparison ---
# Re-run v2 quickly for fair comparison
gs_B <- gibbs_collapsed_sampler(
  y = y, H = obj6$X_fix, R = R,
  n_iter = 6000, n_burn = 1500, n_thin = 1,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 0.3,
  mh_sd_log_sigma2 = 0.3,
  mh_sd_log_tau2   = 0.3,
  init = list(eta = as.vector(obj6$fit$beta),
              sigma2 = max(obj6$fit$sigma2, 0.01),
              tau2 = max(obj6$fit$tau2, 0.01)),
  verbose = TRUE
)

cov_B <- gibbs_curve_coverage(gs_B, obj6, x_grid, truth_f_list,
                               X_raw = X_raw, var_names = var_names)

# --- Baby Bayes (M=6, fixed params) for reference ---
Sigma_true <- sigma2_true * R + tau2_true * diag(n)
bb <- baby_bayes_fit(y = y, H = obj6$X_fix, Sigma = Sigma_true)
cov_A <- bayes_curve_coverage(bb, obj6, x_grid, truth_f_list,
                               X_raw = X_raw, var_names = var_names)

# --- Coverage comparison table ---
cat("\n\n####################################################\n")
cat("# Coverage comparison: Stage A vs B vs C\n")
cat("####################################################\n\n")
compare <- data.frame(
  var       = var_names,
  rmse_A    = cov_A$rmse_curve,
  cov_A     = cov_A$band_coverage,
  rmse_B    = cov_B$rmse_curve,
  cov_B     = cov_B$band_coverage,
  rmse_C    = cov_C$rmse_curve,
  cov_C     = cov_C$band_coverage
)
print(compare, digits = 4)

write.csv(compare, "bayes_out/stage_abc_comparison.csv", row.names = FALSE)

# --- b recovery ---
b_post_C <- colMeans(gs_C$b_samples)
b_post_B <- colMeans(gs_B$b_samples)
cat(sprintf("\nb recovery:  Stage B cor=%.4f RMSE=%.4f  |  Stage C cor=%.4f RMSE=%.4f\n",
            cor(b_true, b_post_B), sqrt(mean((b_true - b_post_B)^2)),
            cor(b_true, b_post_C), sqrt(mean((b_true - b_post_C)^2))))

# --- Plots ---

# Stage C marginal curves
png("bayes_out/stage_c_marginals.png", width = 1200, height = 900)
par(mfrow = c(2, 2))
plot_gibbs_marginals(gs_C, obj_C, x_grid = x_grid,
                      truth_f_list = truth_f_list,
                      var_names = var_names, X_raw = X_raw)
dev.off()

# Stage C trace for sigma2, tau2
png("bayes_out/stage_c_trace.png", width = 1000, height = 700)
plot_gibbs_trace(gs_C, true_sigma2 = sigma2_true, true_tau2 = tau2_true)
dev.off()

# Smoothing variance traces
png("bayes_out/stage_c_tau2s_trace.png", width = 1200, height = 800)
par(mfrow = c(2, 2))
for (j in 1:4) {
  plot(gs_C$tau2_s_samples[, j], type = "l", col = "purple",
       main = paste0("tau2_s: ", var_names[j]),
       ylab = expression(tau[s]^2), xlab = "iteration")
}
dev.off()

cat("\nDone. Check bayes_out/ for:\n")
cat("  stage_abc_comparison.csv   (coverage table)\n")
cat("  stage_c_marginals.png      (M=20 credible bands)\n")
cat("  stage_c_trace.png          (variance traces)\n")
cat("  stage_c_tau2s_trace.png    (smoothing variance traces)\n")
