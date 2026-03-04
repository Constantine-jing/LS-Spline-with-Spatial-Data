# ============================================================
# run_gibbs_v2.R
# Stage B (v2): Collapsed Gibbs (marginalize b)
#
# Fixes the sigma2 explosion from v1 by:
#   - NOT sampling b explicitly in the Gibbs loop
#   - Working with marginal likelihood y ~ N(H eta, sigma2*R + tau2*I)
#   - Using MH for sigma2 and tau2
#   - Computing b as a derived posterior mean
#
# Run: Rscript run_gibbs_v2.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")       # for plotting/summary utilities
source("gibbs_bayes_v2.R")
library(mvtnorm)

M      <- 6
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

simulate_sim1 <- function(n = 400, p = 4, mu = 0, tau2 = 0.15, seed = 42) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  X <- matrix(runif(n*p,0,1), n, p); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu = mu)
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+eps,
             eta=eta, b=rep(0,n), eps=eps, X)
}

stopifnot(ls_tests())

# ============================================================
# Test 1: Sim2, n=400
# ============================================================
cat("\n####################################################\n")
cat("# Collapsed Gibbs (v2): Sim2, n=400\n")
cat("####################################################\n")

sigma2_true <- 0.8; rho_true <- 0.2; tau2_true <- 0.15

dat <- simulate_sim2(n = 400, sigma2 = sigma2_true, rho = rho_true,
                     nu = nu, tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat[, c("x","y_coord")])
X_raw  <- as.matrix(dat[, paste0("X",1:4)])
y      <- dat$Y
b_true <- dat$b

# REML fit
obj <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                      M_vec = rep(M,4), nu = nu,
                      rho_init = rho_true,
                      lambda_init = tau2_true / sigma2_true,
                      verbose = TRUE)
obj$X_raw_for_marginal <- X_raw

# Build R (rho fixed at truth)
D <- pairdist(coords)
R <- matern_cor(D, rho = rho_true, nu = nu)

# Initialize at REML
init <- list(
  eta    = as.vector(obj$fit$beta),
  sigma2 = max(obj$fit$sigma2, 0.01),
  tau2   = max(obj$fit$tau2, 0.01)
)

# Run collapsed Gibbs
cat("\nRunning collapsed Gibbs...\n")
t0 <- proc.time()

gs2 <- gibbs_collapsed_sampler(
  y = y, H = obj$X_fix, R = R,
  n_iter = 5000, n_burn = 1000, n_thin = 1,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,       # IG(2,1): mean=1, mode=0.33
  a_tau   = 2, b_tau   = 0.3,     # IG(2,0.3): mean=0.3, mode=0.1
  mh_sd_log_sigma2 = 0.3,
  mh_sd_log_tau2   = 0.3,
  init = init, verbose = TRUE
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("Done in %.1f seconds.\n", elapsed))

# Summaries
gsm2 <- gibbs_summary(gs2)

cat("\n--- Variance posteriors (v2) ---\n")
cat(sprintf("  sigma2: mean=%.4f  sd=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm2$sigma2["mean"], gsm2$sigma2["sd"],
            gsm2$sigma2["2.5%"], gsm2$sigma2["97.5%"], sigma2_true))
cat(sprintf("  tau2:   mean=%.4f  sd=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm2$tau2["mean"], gsm2$tau2["sd"],
            gsm2$tau2["2.5%"], gsm2$tau2["97.5%"], tau2_true))
cat(sprintf("  MH acceptance: sigma2=%.3f  tau2=%.3f\n",
            gs2$accept_rate["sigma2"], gs2$accept_rate["tau2"]))

# Curve coverage
cov_v2 <- gibbs_curve_coverage(gs2, obj, x_grid, truth_f_list,
                                X_raw = X_raw, var_names = var_names)

# Baby Bayes comparison
Sigma_true <- sigma2_true * R + tau2_true * diag(400)
bb <- baby_bayes_fit(y = y, H = obj$X_fix, Sigma = Sigma_true)
cov_baby <- bayes_curve_coverage(bb, obj, x_grid, truth_f_list,
                                  X_raw = X_raw, var_names = var_names)

cat("\n--- Coverage comparison: Baby vs Gibbs v2 ---\n")
print(data.frame(
  var   = var_names,
  baby  = cov_baby$band_coverage,
  gibbs = cov_v2$band_coverage
), digits = 4)

# b recovery
b_post_mean <- colMeans(gs2$b_samples)
cat(sprintf("\ncor(b_true, b_post_mean) = %.4f\n", cor(b_true, b_post_mean)))
cat(sprintf("RMSE(b) = %.4f\n", sqrt(mean((b_true - b_post_mean)^2))))

# Plots
png("bayes_out/gibbs_v2_marginals_sim2.png", width=1200, height=900)
par(mfrow=c(2,2))
plot_gibbs_marginals(gs2, obj, x_grid = x_grid,
                      truth_f_list = truth_f_list,
                      var_names = var_names, X_raw = X_raw)
dev.off()

png("bayes_out/gibbs_v2_trace_sim2.png", width=1000, height=700)
plot_gibbs_trace(gs2, true_sigma2 = sigma2_true, true_tau2 = tau2_true)
dev.off()

write.csv(data.frame(var = var_names,
                      rmse_baby = cov_baby$rmse_curve,
                      cov_baby  = cov_baby$band_coverage,
                      rmse_v2   = cov_v2$rmse_curve,
                      cov_v2    = cov_v2$band_coverage),
          "bayes_out/gibbs_v2_coverage_sim2.csv", row.names = FALSE)


# ============================================================
# Test 2: Sim1, n=400
# ============================================================
cat("\n####################################################\n")
cat("# Collapsed Gibbs (v2): Sim1, n=400\n")
cat("####################################################\n")

dat1 <- simulate_sim1(n = 400, tau2 = 0.15, seed = seed)
coords1 <- as.matrix(dat1[, c("x","y_coord")])
X_raw1  <- as.matrix(dat1[, paste0("X",1:4)])
y1      <- dat1$Y

rho_upper1 <- 1.0 * max(pairdist(coords1))
obj1 <- fit_ls_spatial(y = y1, X_raw = X_raw1, coords = coords1,
                       M_vec = rep(M,4), nu = nu,
                       rho_init = 0.2, lambda_init = 1.0,
                       rho_upper = rho_upper1, verbose = TRUE)
obj1$X_raw_for_marginal <- X_raw1

D1 <- pairdist(coords1)
R1 <- matern_cor(D1, rho = 0.2, nu = nu)

init1 <- list(
  eta    = as.vector(obj1$fit$beta),
  sigma2 = max(obj1$fit$sigma2, 0.001),
  tau2   = max(obj1$fit$tau2, 0.01)
)

cat("\nRunning collapsed Gibbs (Sim1)...\n")
gs2_sim1 <- gibbs_collapsed_sampler(
  y = y1, H = obj1$X_fix, R = R1,
  n_iter = 5000, n_burn = 1000, n_thin = 1,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 0.3,
  mh_sd_log_sigma2 = 0.3,
  mh_sd_log_tau2   = 0.3,
  init = init1, verbose = TRUE
)

gsm2_sim1 <- gibbs_summary(gs2_sim1)
cat("\n--- Sim1 variance posteriors ---\n")
cat(sprintf("  sigma2: mean=%.4f  95%%CI=[%.4f, %.4f]  (true=0)\n",
            gsm2_sim1$sigma2["mean"], gsm2_sim1$sigma2["2.5%"], gsm2_sim1$sigma2["97.5%"]))
cat(sprintf("  tau2:   mean=%.4f  95%%CI=[%.4f, %.4f]  (true=0.15)\n",
            gsm2_sim1$tau2["mean"], gsm2_sim1$tau2["2.5%"], gsm2_sim1$tau2["97.5%"]))

png("bayes_out/gibbs_v2_trace_sim1.png", width=1000, height=700)
plot_gibbs_trace(gs2_sim1, true_sigma2 = 0, true_tau2 = 0.15)
dev.off()

cat("\n--- ALL DONE ---\n")
cat("Check bayes_out/ for:\n")
cat("  gibbs_v2_marginals_sim2.png  (credible bands)\n")
cat("  gibbs_v2_trace_sim2.png      (trace + density for sigma2, tau2)\n")
cat("  gibbs_v2_trace_sim1.png      (Sim1 trace)\n")
cat("  gibbs_v2_coverage_sim2.csv   (coverage table)\n")
