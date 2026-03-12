# ============================================================
# run_baby_bayes.R
# Stage A: Test Baby Bayes on Sim1 (no spatial b) and Sim2 (spatial b)
#
# What this script does:
#   1. Simulate data (Sim1 or Sim2)
#   2. Fit REML model (your current pipeline)
#   3. Run Baby Bayes with fixed (rho, sigma2, tau2)
#   4. Sanity check: posterior mean ≈ GLS estimate
#   5. Compute credible bands on marginal curves
#   6. Compute coverage (truth inside band?)
#   7. Save plots and tables
#
# Run on hellbender:
#   Rscript run_baby_bayes.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
library(mvtnorm)

# ---- Settings ----
M      <- 6
seed   <- 42
nu     <- 1.5
x_grid <- seq(0, 1, length.out = 101)

truth_f_list <- default_truth_f_list()
var_names    <- paste0("X", 1:4)

dir.create("bayes_out", showWarnings = FALSE)

# ============================================================
# Simulation functions (same as your existing code)
# ============================================================

eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

simulate_sim1 <- function(n = 400, domain = c(0,1,0,1),
                          p = 4, mu = 0, tau2 = 0.15, seed = 42) {
  set.seed(seed)
  coords <- cbind(
    x = runif(n, domain[1], domain[2]),
    y = runif(n, domain[3], domain[4])
  )
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  eta <- eta_truth_additive(X, mu = mu)
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + eps
  data.frame(x = coords[,1], y_coord = coords[,2], Y = y, eta = eta,
             b = rep(0, n), eps = eps, X)
}

simulate_sim2 <- function(n = 400, domain = c(0,1,0,1),
                          p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
  set.seed(seed)
  coords <- cbind(
    x = runif(n, domain[1], domain[2]),
    y = runif(n, domain[3], domain[4])
  )
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  eta <- eta_truth_additive(X, mu = mu)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R + diag(jitter, n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps
  data.frame(x = coords[,1], y_coord = coords[,2], Y = y, eta = eta,
             b = b, eps = eps, X)
}


# ============================================================
# run_baby_bayes_one()
#
# For a single (sim_type, n, use_truth_params) combination:
#   fit REML → run Baby Bayes → compare → plot → save
# ============================================================

run_baby_bayes_one <- function(sim_type = "sim2", n = 400,
                                use_truth_params = TRUE,
                                sigma2_true = 0.8, rho_true = 0.2,
                                tau2_true = 0.15,
                                kappa2 = 1e6) {
  
  cat("\n========================================\n")
  cat(sprintf("Baby Bayes: %s, n=%d, truth_params=%s, kappa2=%.0e\n",
              sim_type, n, use_truth_params, kappa2))
  cat("========================================\n")
  
  # --- 1. Simulate ---
  if (sim_type == "sim1") {
    dat <- simulate_sim1(n = n, tau2 = tau2_true, seed = seed)
    sigma2_true_use <- 0      # Sim1 has no spatial effect
    rho_true_use    <- rho_true
  } else {
    dat <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                         nu = nu, tau2 = tau2_true, seed = seed)
    sigma2_true_use <- sigma2_true
    rho_true_use    <- rho_true
  }
  
  coords <- as.matrix(dat[, c("x", "y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  # --- 2. Fit REML (your current pipeline) ---
  rho_upper <- if (sim_type == "sim1") 1.0 * max(pairdist(coords)) else NULL
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho_true, lambda_init = tau2_true / max(sigma2_true, 0.01),
    rho_upper = rho_upper,
    verbose = TRUE
  )
  obj$X_raw_for_marginal <- X_raw
  
  cat("\nREML estimates:\n")
  cat("  rho    =", obj$fit$rho, "\n")
  cat("  sigma2 =", obj$fit$sigma2, "\n")
  cat("  tau2   =", obj$fit$tau2, "\n")
  
  # --- 3. Build Sigma ---
  D <- pairdist(coords)
  
  if (use_truth_params) {
    # Use true parameters (strongest sanity check)
    if (sigma2_true_use > 0) {
      R_use <- matern_cor(D, rho = rho_true_use, nu = nu)
      Sigma <- sigma2_true_use * R_use + tau2_true * diag(n)
    } else {
      # Sim1: no spatial effect, just nugget
      Sigma <- tau2_true * diag(n)
    }
    param_label <- "truth"
  } else {
    # Use REML estimates
    Sigma <- obj$fit$sigma2 * obj$fit$R + obj$fit$tau2 * diag(n)
    param_label <- "reml"
  }
  
  # --- 4. Run Baby Bayes ---
  bb <- baby_bayes_fit(y = y, H = obj$X_fix, Sigma = Sigma, kappa2 = kappa2)
  
  # --- 5. Sanity check: posterior mean vs GLS ---
  eta_gls <- obj$fit$beta
  max_diff <- max(abs(bb$eta_bar - eta_gls))
  
  cat("\n--- Sanity Check ---\n")
  cat("Max |eta_bar - eta_GLS| =", format(max_diff, digits = 6), "\n")
  
  if (use_truth_params) {
    # When using truth params, GLS uses different Sigma than REML,
    # so they won't match exactly. That's expected.
    cat("  (Using truth params → GLS and Bayes use different Sigma than REML;\n")
    cat("   difference is expected. Check that both are reasonable.)\n")
  } else {
    # When using REML params, should match closely
    if (max_diff > 1e-3) {
      cat("  WARNING: large difference! Check code.\n")
    } else {
      cat("  OK: posterior mean matches GLS.\n")
    }
  }
  
  # --- 6. Credible interval coverage for eta ---
  # We can't easily get true eta (it depends on the realization),
  # but we can check marginal curve coverage
  
  # --- 7. Curve coverage ---
  cov_tab <- bayes_curve_coverage(bb, obj, x_grid, truth_f_list,
                                   X_raw = X_raw, var_names = var_names)
  
  cat("\n--- Curve Metrics ---\n")
  print(cov_tab, digits = 4)
  
  # --- 8. Save plots ---
  tag <- sprintf("%s_n%d_%s_k%.0e", sim_type, n, param_label, kappa2)
  
  png(filename = file.path("bayes_out", paste0("baby_bayes_", tag, ".png")),
      width = 1200, height = 900)
  par(mfrow = c(2, 2))
  plot_bayes_marginals(bb, obj, x_grid = x_grid,
                        truth_f_list = truth_f_list,
                        var_names = var_names,
                        X_raw = X_raw)
  dev.off()
  
  # --- 9. Save tables ---
  write.csv(cov_tab,
            file = file.path("bayes_out", paste0("baby_bayes_coverage_", tag, ".csv")),
            row.names = FALSE)
  
  # --- 10. Return summary ---
  cat("\nSaved to bayes_out/baby_bayes_", tag, ".*\n")
  
  list(bb = bb, obj = obj, cov_tab = cov_tab,
       max_diff_gls = max_diff, tag = tag)
}


# ============================================================
# MAIN: Run the Baby Bayes tests
# ============================================================

cat("\n####################################################\n")
cat("# Baby Bayes Stage A: Sanity checks\n")
cat("####################################################\n")

# ----------------------------------------------------------
# Test 1: Sim2 (spatial b), n=400, TRUE params, kappa2=1e6
# This is the primary sanity check.
# ----------------------------------------------------------
res_sim2_truth <- run_baby_bayes_one(
  sim_type = "sim2", n = 400,
  use_truth_params = TRUE,
  sigma2_true = 0.8, rho_true = 0.2, tau2_true = 0.15,
  kappa2 = 1e6
)

# ----------------------------------------------------------
# Test 2: Sim2, n=400, REML params, kappa2=1e6
# Check that posterior mean matches GLS when same Sigma.
# ----------------------------------------------------------
res_sim2_reml <- run_baby_bayes_one(
  sim_type = "sim2", n = 400,
  use_truth_params = FALSE,
  sigma2_true = 0.8, rho_true = 0.2, tau2_true = 0.15,
  kappa2 = 1e6
)

# ----------------------------------------------------------
# Test 3: Sim1 (no spatial b), n=400, TRUE params
# Sigma = tau2 * I (no spatial correlation).
# ----------------------------------------------------------
res_sim1_truth <- run_baby_bayes_one(
  sim_type = "sim1", n = 400,
  use_truth_params = TRUE,
  sigma2_true = 0, rho_true = 0.2, tau2_true = 0.15,
  kappa2 = 1e6
)

# ----------------------------------------------------------
# Test 4: Prior sensitivity — Sim2, kappa2=1e3 vs 1e6
# Results should be nearly identical (data dominates).
# ----------------------------------------------------------
res_sim2_k3 <- run_baby_bayes_one(
  sim_type = "sim2", n = 400,
  use_truth_params = TRUE,
  sigma2_true = 0.8, rho_true = 0.2, tau2_true = 0.15,
  kappa2 = 1e3
)

# ----------------------------------------------------------
# Summary comparison
# ----------------------------------------------------------
cat("\n\n####################################################\n")
cat("# Summary of all Baby Bayes runs\n")
cat("####################################################\n\n")

cat("--- Test 1: Sim2, truth params, kappa2=1e6 ---\n")
print(res_sim2_truth$cov_tab, digits = 4)

cat("\n--- Test 2: Sim2, REML params, kappa2=1e6 ---\n")
cat("Max |eta_bar - eta_GLS| =", format(res_sim2_reml$max_diff_gls, digits = 6), "\n")
print(res_sim2_reml$cov_tab, digits = 4)

cat("\n--- Test 3: Sim1, truth params, kappa2=1e6 ---\n")
print(res_sim1_truth$cov_tab, digits = 4)

cat("\n--- Test 4: Prior sensitivity (kappa2=1e3 vs 1e6) ---\n")
cat("Max |eta_bar(k3) - eta_bar(k6)| =",
    format(max(abs(res_sim2_k3$bb$eta_bar - res_sim2_truth$bb$eta_bar)), digits = 6), "\n")
cat("kappa2=1e6 coverage:\n")
print(res_sim2_truth$cov_tab[, c("var", "band_coverage")], digits = 4)
cat("kappa2=1e3 coverage:\n")
print(res_sim2_k3$cov_tab[, c("var", "band_coverage")], digits = 4)

cat("\nDone. Check bayes_out/ for plots.\n")
