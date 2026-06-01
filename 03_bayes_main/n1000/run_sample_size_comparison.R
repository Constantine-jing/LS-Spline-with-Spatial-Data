# ============================================================
# run_sample_size_comparison.R
# Side-by-side n=400 vs n=1000 comparison
# on Sim3 and Sim4_plus
#
# Sim3:      spatial confounding (spatially correlated X, spatial b)
# Sim4_plus: same as Sim3 but with 6 extra garbage covariates (p=10)
#
# For each sim x sample size:
#   - Full Bayes (M=20, 10000 iter, tuned MH)
#   - Records: curve RMSE, 95% band coverage, sigma2/tau2/rho recovery, b cor
#
# Story: "Performance improves with sample size as expected"
#
# Run: Rscript run_sample_size_comparison.R
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
# Shared settings
# ============================================================
M_C    <- 20
n_iter <- 10000
n_burn <- 2500
n_thin <- 1
seed   <- 42
nu     <- 1.5

x_grid <- seq(0, 1, length.out = 101)

sigma2_true <- 0.8
rho_true    <- 0.2
tau2_true   <- 0.15
rho_X       <- 0.10
nu_X        <- 1.0

# Tuned MH proposals
mh_sd_sigma2 <- 0.5
mh_sd_tau2   <- 0.3
mh_sd_rho    <- 0.3

sample_sizes <- c(400, 1000)

dir.create("bayes_out", showWarnings = FALSE)
stopifnot(ls_tests())

# ============================================================
# Truth functions
# ============================================================
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

# p=4 truth (Sim3)
truth_f_list_4 <- default_truth_f_list()

# p=10 truth (Sim4_plus): X1-X4 real, X5-X10 zero
truth_f_list_10 <- c(
  default_truth_f_list(),
  replicate(6, function(x) rep(0, length(x)), simplify = FALSE)
)

# ============================================================
# Simulation functions
# ============================================================
simulate_sim3 <- function(n, mu = 0, rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords  <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D       <- pairdist(coords)
  R_X     <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X <- as.matrix(X); colnames(X) <- paste0("X", 1:4)
  eta <- eta_truth_additive(X, mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b   <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = eta + b + eps, eta = eta, b = b, eps = eps, X)
}

simulate_sim4_plus <- function(n, mu = 0, rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                               sigma2 = 0.8, rho = 0.2, nu = 1.5,
                               tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords  <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D       <- pairdist(coords)
  R_X     <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X5_10 <- matrix(runif(n * 6), nrow = n, ncol = 6)
  X <- cbind(X14, X5_10); colnames(X) <- paste0("X", 1:10)
  eta <- eta_truth_additive(X[, 1:4, drop = FALSE], mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b   <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = eta + b + eps, eta = eta, b = b, eps = eps, X)
}

# ============================================================
# Generic runner: one sim x one sample size
# ============================================================
run_one <- function(sim_name, dat, p, truth_f_list, var_names) {

  n      <- nrow(dat)
  coords <- as.matrix(dat[, c("x", "y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y
  b_true <- dat$b
  D      <- pairdist(coords)

  cat(sprintf("\n  --- %s  n=%d  p=%d ---\n", sim_name, n, p))

  des <- ls_additive_build(X_raw, M_vec = rep(M_C, p))
  H   <- cbind(1, des$W)

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

  t0 <- proc.time()
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
  elapsed <- (proc.time() - t0)[3]

  gsm         <- gibbs_summary(gs)
  rho_mean    <- mean(gs$rho_samples)
  tau2s_means <- colMeans(gs$tau2_s_samples)
  names(tau2s_means) <- var_names
  b_post <- colMeans(gs$b_samples)

  obj_tmp <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)
  cov_tab <- gibbs_curve_coverage(gs, obj_tmp, x_grid, truth_f_list,
                                   X_raw = X_raw, var_names = var_names)

  cat(sprintf("    sigma2=%.4f  tau2=%.4f  rho=%.4f  b_cor=%.4f  time=%.0fs\n",
              gsm$sigma2["mean"], gsm$tau2["mean"], rho_mean,
              cor(b_true, b_post), elapsed))
  cat(sprintf("    MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))
  cat("    Curve metrics:\n"); print(cov_tab, digits = 4)

  if (p > 4) {
    cat("    tau2_s (garbage = X5+):\n"); print(round(tau2s_means, 5))
  }

  # Marginal curves plot
  nrow_plot <- ceiling(p / 2)
  pdf(sprintf("bayes_out/ncomp_%s_n%d_marginals.pdf", sim_name, n),
      width = 12, height = nrow_plot * 3)
  par(mfrow = c(nrow_plot, 2))
  plot_gibbs_marginals(gs, obj_tmp, x_grid = x_grid,
                       truth_f_list = truth_f_list,
                       var_names = var_names, X_raw = X_raw)
  dev.off()

  list(gsm = gsm, rho_mean = rho_mean, tau2s_means = tau2s_means,
       cov_tab = cov_tab, b_cor = cor(b_true, b_post), elapsed = elapsed,
       accept = gs$accept_rate)
}

# ============================================================
# Run all 4 combinations
# ============================================================
all_results <- list()

for (n_val in sample_sizes) {
  cat(sprintf("\n\n####################################################\n"))
  cat(sprintf("# Sim3  n=%d\n", n_val))
  cat(sprintf("####################################################\n"))
  dat3 <- simulate_sim3(n = n_val, rho_X = rho_X, nu_X = nu_X,
                        sigma2 = sigma2_true, rho = rho_true, nu = nu,
                        tau2 = tau2_true, seed = seed)
  key <- sprintf("sim3_n%d", n_val)
  all_results[[key]] <- run_one(key, dat3, p = 4,
                                truth_f_list = truth_f_list_4,
                                var_names = paste0("X", 1:4))
}

for (n_val in sample_sizes) {
  cat(sprintf("\n\n####################################################\n"))
  cat(sprintf("# Sim4_plus  n=%d\n", n_val))
  cat(sprintf("####################################################\n"))
  dat4p <- simulate_sim4_plus(n = n_val, rho_X = rho_X, nu_X = nu_X,
                              sigma2 = sigma2_true, rho = rho_true, nu = nu,
                              tau2 = tau2_true, seed = seed)
  key <- sprintf("sim4plus_n%d", n_val)
  all_results[[key]] <- run_one(key, dat4p, p = 10,
                                truth_f_list = truth_f_list_10,
                                var_names = paste0("X", 1:10))
}

# ============================================================
# Grand comparison table
# ============================================================
cat("\n\n####################################################\n")
cat("# SAMPLE SIZE COMPARISON: n=400 vs n=1000\n")
cat("####################################################\n\n")

# --- Variance recovery ---
cat("--- Variance recovery ---\n")
cat(sprintf("%-20s  %8s  %8s  %8s  %8s\n", "Setting", "sigma2", "tau2", "rho", "b_cor"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (key in names(all_results)) {
  r <- all_results[[key]]
  cat(sprintf("%-20s  %8.4f  %8.4f  %8.4f  %8.4f\n",
              key, r$gsm$sigma2["mean"], r$gsm$tau2["mean"],
              r$rho_mean, r$b_cor))
}
cat(sprintf("%-20s  %8.4f  %8.4f  %8.4f\n", "TRUE",
            sigma2_true, tau2_true, rho_true))

# --- Curve RMSE for X1-X4 (comparable across sims) ---
cat("\n--- Curve RMSE: X1-X4 (real covariates) ---\n")
cat(sprintf("%-20s  %8s  %8s  %8s  %8s\n", "Setting", "X1", "X2", "X3", "X4"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (key in names(all_results)) {
  r <- all_results[[key]]
  rmse4 <- r$cov_tab$rmse_curve[1:4]
  cat(sprintf("%-20s  %8.5f  %8.5f  %8.5f  %8.5f\n",
              key, rmse4[1], rmse4[2], rmse4[3], rmse4[4]))
}

# --- Band coverage for X1-X4 ---
cat("\n--- 95% band coverage: X1-X4 ---\n")
cat(sprintf("%-20s  %8s  %8s  %8s  %8s\n", "Setting", "X1", "X2", "X3", "X4"))
cat(paste(rep("-", 60), collapse=""), "\n")
for (key in names(all_results)) {
  r <- all_results[[key]]
  cov4 <- r$cov_tab$band_coverage[1:4]
  cat(sprintf("%-20s  %8.3f  %8.3f  %8.3f  %8.3f\n",
              key, cov4[1], cov4[2], cov4[3], cov4[4]))
}

# --- tau2_s for Sim4_plus (garbage shrinkage) ---
cat("\n--- tau2_s shrinkage in Sim4_plus (X5-X10 = garbage) ---\n")
for (n_val in sample_sizes) {
  key <- sprintf("sim4plus_n%d", n_val)
  r   <- all_results[[key]]
  cat(sprintf("  n=%d  real(X1-X4): %s\n", n_val,
              paste(sprintf("%.4f", r$tau2s_means[1:4]), collapse=" ")))
  cat(sprintf("        garbage(X5-X10): %s\n",
              paste(sprintf("%.4f", r$tau2s_means[5:10]), collapse=" ")))
}

# --- Timing ---
cat("\n--- Timing ---\n")
for (key in names(all_results)) {
  cat(sprintf("  %-20s  %.0f seconds\n", key, all_results[[key]]$elapsed))
}

# ============================================================
# Save CSV
# ============================================================
rows <- lapply(names(all_results), function(key) {
  r   <- all_results[[key]]
  sim <- sub("_n[0-9]+", "", key)
  n_v <- as.integer(sub(".*_n", "", key))
  p_v <- nrow(r$cov_tab)
  rmse_all <- setNames(r$cov_tab$rmse_curve,  paste0("rmse_X",  1:p_v))
  cov_all  <- setNames(r$cov_tab$band_coverage, paste0("cov_X", 1:p_v))
  # Pad to 10 columns for consistent CSV
  rmse_pad <- c(rmse_all, rep(NA, 10 - p_v))
  cov_pad  <- c(cov_all,  rep(NA, 10 - p_v))
  names(rmse_pad) <- paste0("rmse_X", 1:10)
  names(cov_pad)  <- paste0("cov_X",  1:10)
  tau2s_pad <- c(r$tau2s_means, rep(NA, 10 - p_v))
  names(tau2s_pad) <- paste0("tau2s_X", 1:10)
  data.frame(setting  = key, sim = sim, n = n_v,
             sigma2   = r$gsm$sigma2["mean"],
             tau2     = r$gsm$tau2["mean"],
             rho      = r$rho_mean,
             b_cor    = r$b_cor,
             time_sec = r$elapsed,
             t(rmse_pad), t(cov_pad), t(tau2s_pad),
             row.names = NULL)
})
comp_df <- do.call(rbind, rows)
write.csv(comp_df, "bayes_out/sample_size_comparison.csv", row.names = FALSE)
cat("\nSaved: bayes_out/sample_size_comparison.csv\n")

cat("\nAll marginal plots in bayes_out/ncomp_*.pdf\n")
