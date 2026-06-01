# ============================================================
# run_sim4_n1000.R
# Sim4: 4 real covariates + 2 garbage covariates
#   n = 1000, M = 20 knots
#   Full Bayes (rho sampled, RW2), 10000 MCMC iterations
#   vs REML (M=6)
#   Outputs: marginal curves + trace plots + comparison CSV
#
# Run: Rscript run_sim4_n1000.R
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
M_C    <- 20          # knots for full Bayes
M_reml <- 6           # knots for REML
n      <- 1000        # sample size
n_iter <- 10000       # MCMC iterations
n_burn <- 2500        # burn-in (25%)
n_thin <- 1
seed   <- 42
nu     <- 1.5
p      <- 6           # 4 real + 2 garbage

x_grid <- seq(0, 1, length.out = 101)

# True parameter values
sigma2_true <- 0.8
rho_true    <- 0.2
tau2_true   <- 0.15
rho_X       <- 0.10
nu_X        <- 1.0

# Variable names
var_names <- paste0("X", 1:6)

# Truth functions: X1-X4 are real, X5-X6 are zero
truth_f_list_4 <- default_truth_f_list()   # f1..f4
truth_f_list_6 <- c(truth_f_list_4,
                    list(function(x) rep(0, length(x)),
                         function(x) rep(0, length(x))))

dir.create("bayes_out", showWarnings = FALSE)

stopifnot(ls_tests())

# ============================================================
# Simulation: Sim4
#   - X1-X4: spatially correlated (Matern), rank-uniform-transformed
#   - X5-X6: iid Uniform (garbage predictors)
#   - b:     spatial random effect (Matern)
#   - eta:   only uses X1-X4
# ============================================================
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim4 <- function(n = 1000, p = 6, mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D      <- pairdist(coords)

  # Spatially correlated X1-X4
  R_X    <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })

  # Garbage X5-X6
  X56 <- cbind(runif(n), runif(n))
  X   <- cbind(X14, X56)
  colnames(X) <- paste0("X", 1:p)

  # Mean surface
  eta <- eta_truth_additive(X[, 1:4, drop = FALSE], mu = mu)

  # Spatial random effect
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b   <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))

  # Noise
  eps <- rnorm(n, 0, sqrt(tau2))

  data.frame(x = coords[, 1], y_coord = coords[, 2],
             Y = eta + b + eps, eta = eta, b = b, eps = eps, X)
}

# ============================================================
# Generate data
# ============================================================
cat(sprintf("\n=== Simulating Sim4: n=%d, p=%d (4 real + 2 garbage) ===\n", n, p))
dat4   <- simulate_sim4(n = n, p = p, rho_X = rho_X, nu_X = nu_X,
                        sigma2 = sigma2_true, rho = rho_true, nu = nu,
                        tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat4[, c("x", "y_coord")])
X_raw  <- as.matrix(dat4[, paste0("X", 1:p)])
y      <- dat4$Y
b_true <- dat4$b
D      <- pairdist(coords)

cat(sprintf("  Y range: [%.3f, %.3f]  mean=%.3f\n", min(y), max(y), mean(y)))


# ============================================================
# REML (M=6)
# ============================================================
cat(sprintf("\n####################################################\n"))
cat("# REML: M=6\n")
cat("####################################################\n")

t0_reml <- proc.time()
reml_obj <- fit_ls_spatial(
  y       = y,
  X_raw   = X_raw,
  coords  = coords,
  M_vec   = rep(M_reml, p),
  nu      = nu,
  rho_init    = rho_true,
  lambda_init = tau2_true / sigma2_true,
  verbose = TRUE
)
t_reml <- (proc.time() - t0_reml)[3]
reml_obj$X_raw_for_marginal <- X_raw

cat(sprintf("  REML done in %.1fs\n", t_reml))
cat(sprintf("  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho))

# REML marginal curves
reml_curves <- marginal_curves(reml_obj, x_grid = x_grid,
                               truth_f_list = truth_f_list_6, clip = TRUE)
reml_err    <- curve_error_table(reml_curves, var_names = var_names)
reml_rmse   <- reml_err$rmse_curve
cat("  REML curve RMSE:\n")
print(data.frame(var = var_names, rmse = round(reml_rmse, 5)))

# REML marginal plot
pdf("bayes_out/sim4_n1000_reml_marginals.pdf", width = 12, height = 9)
par(mfrow = c(3, 2))
plot_marginals(reml_curves, var_names = var_names, show_rmse = TRUE,
               ylim_shared = FALSE)
dev.off()
cat("  Marginal plot: bayes_out/sim4_n1000_reml_marginals.pdf\n")


# ============================================================
# Full Bayes (M=20, n_iter=10000)
# ============================================================
cat(sprintf("\n####################################################\n"))
cat(sprintf("# Full Bayes: M=%d, n_iter=%d, n_burn=%d\n", M_C, n_iter, n_burn))
cat("####################################################\n")

des <- ls_additive_build(X_raw, M_vec = rep(M_C, p))
H   <- cbind(1, des$W)
cat(sprintf("  Design matrix H: %d x %d\n", nrow(H), ncol(H)))

# --- Initialize from REML (M=6) ---
obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                       M_vec   = rep(6, p), nu = nu,
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

init <- list(
  eta    = eta_init,
  sigma2 = max(obj6$fit$sigma2, 0.01),
  tau2   = max(obj6$fit$tau2, 0.01),
  rho    = obj6$fit$rho,
  tau2_s = rep(1.0, p)
)
cat(sprintf("  REML init: rho=%.4f  sigma2=%.4f  tau2=%.4f\n",
            init$rho, init$sigma2, init$tau2))

# --- Run Gibbs ---
t0_bayes <- proc.time()
gs <- gibbs_full_sampler(
  y        = y,
  H        = H,
  D        = D,
  nu       = nu,
  col_map  = des$col_map,
  n_iter   = n_iter,
  n_burn   = n_burn,
  n_thin   = n_thin,
  kappa2          = 1e6,
  a_sigma         = 2,    b_sigma  = 1,
  a_tau           = 2,    b_tau    = 0.3,
  a_smooth        = 1,    b_smooth = 0.005,
  log_rho_mu      = log(0.2),
  log_rho_sd      = 1.0,
  mh_sd_log_sigma2 = 0.3,
  mh_sd_log_tau2   = 0.3,
  mh_sd_log_rho    = 0.2,
  init    = init,
  verbose = TRUE
)
t_bayes <- (proc.time() - t0_bayes)[3]
cat(sprintf("  Full Bayes done in %.1fs\n", t_bayes))

# --- Posterior summaries ---
gsm <- gibbs_summary(gs)
rho_summary <- c(mean  = mean(gs$rho_samples),
                 sd    = sd(gs$rho_samples),
                 quantile(gs$rho_samples, c(0.025, 0.5, 0.975)))

cat(sprintf("  sigma2: mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"], sigma2_true))
cat(sprintf("  tau2:   mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            gsm$tau2["mean"], gsm$tau2["2.5%"], gsm$tau2["97.5%"], tau2_true))
cat(sprintf("  rho:    mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
            rho_summary["mean"], rho_summary["2.5%"], rho_summary["97.5%"], rho_true))

cat("  tau2_s (smoothing variances):\n")
tau2s_means <- colMeans(gs$tau2_s_samples)
names(tau2s_means) <- var_names
print(round(tau2s_means, 5))

cat(sprintf("  MH accept rates: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))

# --- b recovery ---
b_post <- colMeans(gs$b_samples)
cat(sprintf("  Spatial b: cor=%.4f  RMSE=%.4f\n",
            cor(b_true, b_post), sqrt(mean((b_true - b_post)^2))))

# --- Build obj for curve evaluation ---
obj_bayes <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)

# --- Curve coverage ---
cov_tab <- gibbs_curve_coverage(gs, obj_bayes, x_grid, truth_f_list_6,
                                X_raw = X_raw, var_names = var_names)
cat("  Curve metrics (Bayes):\n")
print(cov_tab, digits = 4)

# --- Marginal curves plot ---
pdf("bayes_out/sim4_n1000_bayes_marginals.pdf", width = 12, height = 9)
par(mfrow = c(3, 2))
plot_gibbs_marginals(gs, obj_bayes, x_grid = x_grid,
                     truth_f_list = truth_f_list_6,
                     var_names = var_names, X_raw = X_raw)
dev.off()
cat("  Marginal plot: bayes_out/sim4_n1000_bayes_marginals.pdf\n")

# --- Trace plots ---
pdf("bayes_out/sim4_n1000_bayes_trace.pdf", width = 10, height = 10)
plot_gibbs_trace_full(gs, true_sigma2 = sigma2_true,
                          true_tau2   = tau2_true,
                          true_rho    = rho_true)
dev.off()
cat("  Trace plot: bayes_out/sim4_n1000_bayes_trace.pdf\n")

# ============================================================
# REML vs Bayes: Comparison summary
# ============================================================
cat("\n\n####################################################\n")
cat("# REML vs Full Bayes: Sim4 n=1000 Summary\n")
cat("####################################################\n\n")

cat("--- Curve RMSE ---\n")
cat(sprintf("%-10s", ""))
for (v in var_names) cat(sprintf("  %8s", v))
cat("\n")
cat(sprintf("%-10s", "REML"))
for (r in reml_rmse) cat(sprintf("  %8.5f", r))
cat("\n")
cat(sprintf("%-10s", "Bayes"))
for (r in cov_tab$rmse_curve) cat(sprintf("  %8.5f", r))
cat("\n")

cat("\n--- Bayes Curve Coverage (95% CI) ---\n")
print(cov_tab, digits = 4)

cat("\n--- Variance Estimates ---\n")
cat(sprintf("  REML:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho))
cat(sprintf("  Bayes: sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            gsm$sigma2["mean"], gsm$tau2["mean"], rho_summary["mean"]))
cat(sprintf("  True:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            sigma2_true, tau2_true, rho_true))

cat(sprintf("\n--- Timing ---\n"))
cat(sprintf("  REML:       %.1fs\n", t_reml))
cat(sprintf("  Full Bayes: %.1fs\n", t_bayes))

# ============================================================
# Save comparison CSV
# ============================================================
comp_table <- data.frame(
  method   = c("REML", "Bayes"),
  rmse_X1  = c(reml_rmse[1], cov_tab$rmse_curve[1]),
  rmse_X2  = c(reml_rmse[2], cov_tab$rmse_curve[2]),
  rmse_X3  = c(reml_rmse[3], cov_tab$rmse_curve[3]),
  rmse_X4  = c(reml_rmse[4], cov_tab$rmse_curve[4]),
  rmse_X5  = c(reml_rmse[5], cov_tab$rmse_curve[5]),
  rmse_X6  = c(reml_rmse[6], cov_tab$rmse_curve[6]),
  sigma2   = c(reml_obj$fit$sigma2,  gsm$sigma2["mean"]),
  tau2     = c(reml_obj$fit$tau2,    gsm$tau2["mean"]),
  rho      = c(reml_obj$fit$rho,     rho_summary["mean"]),
  time_sec = c(t_reml, t_bayes)
)
write.csv(comp_table, "bayes_out/sim4_n1000_reml_vs_bayes.csv", row.names = FALSE)
cat("\nSaved: bayes_out/sim4_n1000_reml_vs_bayes.csv\n")
write.csv(cov_tab,   "bayes_out/sim4_n1000_bayes_coverage.csv", row.names = FALSE)
cat("Saved: bayes_out/sim4_n1000_bayes_coverage.csv\n")

cat("\nAll outputs in bayes_out/\n")
cat("  sim4_n1000_reml_marginals.pdf\n")
cat("  sim4_n1000_bayes_marginals.pdf\n")
cat("  sim4_n1000_bayes_trace.pdf\n")
cat("  sim4_n1000_reml_vs_bayes.csv\n")
cat("  sim4_n1000_bayes_coverage.csv\n")
