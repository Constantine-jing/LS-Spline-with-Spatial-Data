# ============================================================
# run_gibbs_bayes.R
# Stage B: Test Gibbs sampler on Sim1 and Sim2
#
# Compares with Stage A (Baby Bayes) to show:
#   - credible bands widen (better coverage)
#   - variance components recovered
#   - spatial effect b recovered
#
# Run on hellbender:
#   Rscript run_gibbs_bayes.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
library(mvtnorm)

# ---- Settings ----
M      <- 6
seed   <- 42
nu     <- 1.5
x_grid <- seq(0, 1, length.out = 101)

truth_f_list <- default_truth_f_list()
var_names    <- paste0("X", 1:4)

dir.create("bayes_out", showWarnings = FALSE)

# ---- Simulation functions (same as before) ----
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2] - 0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim1 <- function(n = 400, p = 4, mu = 0, tau2 = 0.15, seed = 42) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  X <- matrix(runif(n*p, 0, 1), n, p); colnames(X) <- paste0("X", 1:p)
  eta <- eta_truth_additive(X, mu = mu)
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x = coords[,1], y_coord = coords[,2], Y = eta + eps,
             eta = eta, b = rep(0, n), eps = eps, X)
}

simulate_sim2 <- function(n = 400, p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  X <- matrix(runif(n*p, 0, 1), n, p); colnames(X) <- paste0("X", 1:p)
  eta <- eta_truth_additive(X, mu = mu)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R + diag(jitter, n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x = coords[,1], y_coord = coords[,2], Y = eta + b + eps,
             eta = eta, b = b, eps = eps, X)
}


# ============================================================
# run_gibbs_one()
# ============================================================
run_gibbs_one <- function(sim_type = "sim2", n = 400,
                           sigma2_true = 0.8, rho_true = 0.2,
                           tau2_true = 0.15,
                           n_iter = 5000, n_burn = 1000, n_thin = 1,
                           a_sigma = 1, b_sigma = 0.005,
                           a_tau = 1, b_tau = 0.005) {
  
  cat("\n========================================\n")
  cat(sprintf("Gibbs Stage B: %s, n=%d\n", sim_type, n))
  cat(sprintf("  n_iter=%d, n_burn=%d, n_thin=%d\n", n_iter, n_burn, n_thin))
  cat("========================================\n")
  
  # --- 1. Simulate ---
  if (sim_type == "sim1") {
    dat <- simulate_sim1(n = n, tau2 = tau2_true, seed = seed)
    sigma2_true_use <- 0
  } else {
    dat <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                         nu = nu, tau2 = tau2_true, seed = seed)
    sigma2_true_use <- sigma2_true
  }
  
  coords <- as.matrix(dat[, c("x", "y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  b_true <- dat$b
  
  # --- 2. REML fit (for initialization + comparison) ---
  rho_upper <- if (sim_type == "sim1") 1.0 * max(pairdist(coords)) else NULL
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho_true, lambda_init = tau2_true / max(sigma2_true, 0.01),
    rho_upper = rho_upper, verbose = TRUE
  )
  obj$X_raw_for_marginal <- X_raw
  
  # --- 3. Build R matrix (rho FIXED at true value) ---
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho_true, nu = nu)
  
  # --- 4. Initialize Gibbs at REML estimates ---
  init <- list(
    eta    = as.vector(obj$fit$beta),
    b      = obj$b_hat,
    sigma2 = max(obj$fit$sigma2, 0.01),
    tau2   = max(obj$fit$tau2, 0.01)
  )
  
  # --- 5. Run Gibbs ---
  cat("\nRunning Gibbs sampler...\n")
  t0 <- proc.time()
  
  gs <- gibbs_sampler(
    y = y, H = obj$X_fix, R = R,
    n_iter = n_iter, n_burn = n_burn, n_thin = n_thin,
    a_sigma = a_sigma, b_sigma = b_sigma,
    a_tau = a_tau, b_tau = b_tau,
    init = init, verbose = TRUE
  )
  
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("Gibbs done in %.1f seconds.\n", elapsed))
  
  # --- 6. Summaries ---
  gsm <- gibbs_summary(gs)
  
  cat("\n--- Variance component posteriors ---\n")
  cat(sprintf("  sigma2: mean=%.4f  sd=%.4f  95%% CI=[%.4f, %.4f]  (true=%.4f)\n",
              gsm$sigma2["mean"], gsm$sigma2["sd"],
              gsm$sigma2["2.5%"], gsm$sigma2["97.5%"], sigma2_true_use))
  cat(sprintf("  tau2:   mean=%.4f  sd=%.4f  95%% CI=[%.4f, %.4f]  (true=%.4f)\n",
              gsm$tau2["mean"], gsm$tau2["sd"],
              gsm$tau2["2.5%"], gsm$tau2["97.5%"], tau2_true))
  
  # --- 7. Curve coverage (Gibbs) ---
  cov_gibbs <- gibbs_curve_coverage(gs, obj, x_grid, truth_f_list,
                                     X_raw = X_raw, var_names = var_names)
  
  cat("\n--- Gibbs curve metrics ---\n")
  print(cov_gibbs, digits = 4)
  
  # --- 8. Baby Bayes comparison (Stage A with true params) ---
  if (sigma2_true_use > 0) {
    Sigma_true <- sigma2_true_use * R + tau2_true * diag(n)
  } else {
    Sigma_true <- tau2_true * diag(n)
  }
  bb <- baby_bayes_fit(y = y, H = obj$X_fix, Sigma = Sigma_true)
  cov_baby <- bayes_curve_coverage(bb, obj, x_grid, truth_f_list,
                                    X_raw = X_raw, var_names = var_names)
  
  cat("\n--- Baby Bayes (Stage A) curve metrics for comparison ---\n")
  print(cov_baby, digits = 4)
  
  # --- 9. Coverage improvement ---
  cat("\n--- Coverage improvement (Gibbs - Baby) ---\n")
  diff_cov <- cov_gibbs$band_coverage - cov_baby$band_coverage
  print(data.frame(var = var_names, baby = cov_baby$band_coverage,
                    gibbs = cov_gibbs$band_coverage, diff = diff_cov), digits = 4)
  
  # --- 10. Spatial effect recovery (Sim2 only) ---
  if (sigma2_true_use > 0) {
    b_post_mean <- colMeans(gs$b_samples)
    cat(sprintf("\n--- Spatial effect b recovery ---\n"))
    cat(sprintf("  cor(b_true, b_post_mean) = %.4f\n", cor(b_true, b_post_mean)))
    cat(sprintf("  RMSE(b_true, b_post_mean) = %.4f\n",
                sqrt(mean((b_true - b_post_mean)^2))))
    cat(sprintf("  cor(b_true, b_BLUP) = %.4f\n", cor(b_true, obj$b_hat)))
  }
  
  # --- 11. Save plots ---
  tag <- sprintf("%s_n%d", sim_type, n)
  
  # Marginal curves with Gibbs bands
  png(file.path("bayes_out", paste0("gibbs_marginals_", tag, ".png")),
      width = 1200, height = 900)
  par(mfrow = c(2, 2))
  plot_gibbs_marginals(gs, obj, x_grid = x_grid,
                        truth_f_list = truth_f_list,
                        var_names = var_names, X_raw = X_raw)
  dev.off()
  
  # Trace plots
  png(file.path("bayes_out", paste0("gibbs_trace_", tag, ".png")),
      width = 1000, height = 700)
  plot_gibbs_trace(gs, true_sigma2 = sigma2_true_use, true_tau2 = tau2_true)
  dev.off()
  
  # Save coverage table
  cov_compare <- data.frame(
    var   = var_names,
    rmse_baby  = cov_baby$rmse_curve,
    cov_baby   = cov_baby$band_coverage,
    rmse_gibbs = cov_gibbs$rmse_curve,
    cov_gibbs  = cov_gibbs$band_coverage
  )
  write.csv(cov_compare,
            file.path("bayes_out", paste0("gibbs_vs_baby_", tag, ".csv")),
            row.names = FALSE)
  
  cat("\nSaved to bayes_out/gibbs_*_", tag, ".*\n")
  
  list(gs = gs, gsm = gsm, obj = obj, bb = bb,
       cov_gibbs = cov_gibbs, cov_baby = cov_baby,
       elapsed = elapsed, tag = tag)
}


# ============================================================
# MAIN
# ============================================================

cat("\n####################################################\n")
cat("# Gibbs Stage B: Sampling variance components\n")
cat("####################################################\n")

stopifnot(ls_tests())

# ----------------------------------------------------------
# Test 1: Sim2 (spatial b), n=400
# Main test: should improve X1 coverage over Baby Bayes.
# ----------------------------------------------------------
res_sim2 <- run_gibbs_one(
  sim_type = "sim2", n = 400,
  sigma2_true = 0.8, rho_true = 0.2, tau2_true = 0.15,
  n_iter = 5000, n_burn = 1000, n_thin = 1
)

# ----------------------------------------------------------
# Test 2: Sim1 (no spatial b), n=400
# Stress test: sigma2 should concentrate near 0.
# ----------------------------------------------------------
res_sim1 <- run_gibbs_one(
  sim_type = "sim1", n = 400,
  sigma2_true = 0, rho_true = 0.2, tau2_true = 0.15,
  n_iter = 5000, n_burn = 1000, n_thin = 1
)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------
cat("\n\n####################################################\n")
cat("# Stage B Summary\n")
cat("####################################################\n\n")

cat("--- Sim2: Coverage comparison ---\n")
print(data.frame(
  var   = var_names,
  baby  = res_sim2$cov_baby$band_coverage,
  gibbs = res_sim2$cov_gibbs$band_coverage
), digits = 4)

cat("\n--- Sim2: Variance recovery ---\n")
cat(sprintf("  sigma2: posterior mean=%.4f (true=0.80)\n",
            res_sim2$gsm$sigma2["mean"]))
cat(sprintf("  tau2:   posterior mean=%.4f (true=0.15)\n",
            res_sim2$gsm$tau2["mean"]))

cat("\n--- Sim1: sigma2 near zero? ---\n")
cat(sprintf("  sigma2: posterior mean=%.4f  95%% upper=%.4f (true=0)\n",
            res_sim1$gsm$sigma2["mean"], res_sim1$gsm$sigma2["97.5%"]))

cat("\nDone. Check bayes_out/ for plots and tables.\n")
