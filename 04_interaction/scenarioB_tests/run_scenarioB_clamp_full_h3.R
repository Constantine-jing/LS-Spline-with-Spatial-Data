# ============================================================
# run_scenarioB_clamp_full_h3.R
#
# PLAN C: clamp ALL variance components and smoothing variances
# -------------------------------------------------------------
# Refit Scenario B, seed 1, with:
#   - tau2_s_main[1,2,3] clamped at GLS-optimal (Plan A values)
#   - tau2_s_int[1_2, 1_3, 2_3] clamped at 0.001  (Plan B)
#   - sigma2, tau2 (nugget), rho clamped at REML estimates (NEW)
#
# Rationale:
#   Plan A and Plan B both gave RMSE_f1 ~ 0.283. Plan B falsified
#   per-iteration interaction coupling. The remaining suspect is
#   per-iteration Sigma fluctuation.
#
#   In the ortho=TRUE MCMC, sigma2 traced 0.5 to 1.4, tau2 traced
#   1.0 to 1.6, rho traced 0.05 to 0.19. The GLS test that gave
#   RMSE_f1 = 0.066 used FIXED REML Sigma. If we clamp sigma2,
#   tau2, rho at REML values too, the MCMC effectively has a fixed
#   Sigma. Only beta and the spatial RE b are stochastic.
#
# Decision rule:
#   - RMSE_f1 (full) drops to ~0.07 (matches GLS):
#       -> Sigma fluctuation CONFIRMED as source of bias gap.
#          Mechanism: MCMC averages over different Sigma_t,
#          producing a posterior mean beta_1 different from
#          GLS at REML Sigma.
#   - RMSE_f1 stays > 0.20:
#       -> Sigma fluctuation FALSIFIED. The bias has to come
#          from beta exploration (NOT Sigma exploration).
#          That would point to multimodality in the conditional
#          posterior of beta_1, which would be very odd given
#          it's Gaussian, but worth investigating.
#   - In between: partial.
#
# Settings:
#   - Same DGP (seed = 1), same MCMC budget
#   - orthogonalize = TRUE
#   - All clamps active
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")

out_dir <- file.path("comparison", "output", "scenarioB_h3_clamp_full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" PLAN C: clamp ALL variance components AND smoothing variances\n")
cat("================================================================\n\n")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n", sim$n, sim$p))

# Get REML Sigma values (same as the diagnostic GLS used)
cat("\nFitting REML to extract Sigma clamp values...\n")
loc <- as.matrix(sim$data[, c("lon", "lat")])
obj6 <- fit_ls_spatial(
  y     = sim$data$y,
  X_raw = cbind(sim$data$X1, sim$data$X2, sim$data$X3),
  coords = loc,
  M_vec = rep(6, 3), nu = 1.0,
  rho_init = 0.06, lambda_init = 1.0,
  verbose = FALSE
)
sigma2_clamp <- max(obj6$fit$sigma2, 0.01)
tau2_clamp   <- max(obj6$fit$tau2,   0.01)
rho_clamp    <- obj6$fit$rho

cat(sprintf("REML estimates (used as clamps):\n"))
cat(sprintf("  sigma2 = %.4f\n", sigma2_clamp))
cat(sprintf("  tau2   = %.4f  (nugget)\n", tau2_clamp))
cat(sprintf("  rho    = %.4f\n\n", rho_clamp))

# All clamp values
clamp_main <- c(0.00625, 0.00250, 0.00125)
clamp_int  <- c(0.001, 0.001, 0.001)

cat("Smoothing-variance clamps:\n")
for (j in 1:3) cat(sprintf("  tau2_s_main[%d] = %.5f\n", j, clamp_main[j]))
for (k in 1:3) cat(sprintf("  tau2_s_int[%d]  = %.5f\n", k, clamp_int[k]))
cat("\n")

cat("With sigma2/tau2/rho also clamped, MCMC has a FIXED Sigma.\n")
cat("Only beta and the spatial RE b are stochastic now.\n\n")

# Fit
t0 <- proc.time()
fit_C <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize     = TRUE,
    tau2_s_main_fixed = clamp_main,
    tau2_s_int_fixed  = clamp_int,
    sigma2_fixed      = sigma2_clamp,
    tau2_fixed        = tau2_clamp,
    rho_fixed         = rho_clamp,
    verbose           = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

out_rds <- file.path(out_dir, "fit_h3_clamp_full_seed1.rds")
saveRDS(fit_C, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Sanity: clamps held ----
cat("---- Sanity: posterior means of clamped values ----\n")
cat(sprintf("  sigma2  mean=%.5f  (clamp=%.5f) %s\n",
            mean(fit_C$var_comp$sigma2), sigma2_clamp,
            if (abs(mean(fit_C$var_comp$sigma2) - sigma2_clamp) < 1e-10) "OK" else "FAIL"))
cat(sprintf("  tau2    mean=%.5f  (clamp=%.5f) %s\n",
            mean(fit_C$var_comp$tau2_s),  tau2_clamp,
            if (abs(mean(fit_C$var_comp$tau2_s) - tau2_clamp) < 1e-10) "OK" else "FAIL"))
cat(sprintf("  rho     mean=%.5f  (clamp=%.5f) %s\n",
            mean(fit_C$var_comp$rho), rho_clamp,
            if (abs(mean(fit_C$var_comp$rho) - rho_clamp) < 1e-10) "OK" else "FAIL"))
for (j in 1:3) {
  cat(sprintf("  tau2_s_main[%d]  mean=%.5f  (clamp=%.5f) %s\n",
              j, mean(fit_C$var_comp$tau2[[j]]), clamp_main[j],
              if (abs(mean(fit_C$var_comp$tau2[[j]]) - clamp_main[j]) < 1e-10) "OK" else "FAIL"))
}
cat("\n")

# ---- Verdict ----
fit_y_path <- "fit_scenarioB_seed1_ortho.rds"
fit_a_path <- "comparison/output/scenarioB_h3_clamp_all/fit_h3_clamp_all_seed1.rds"
fit_b_path <- "comparison/output/scenarioB_h3_clamp_main_int/fit_h3_clamp_main_int_seed1.rds"

if (file.exists(fit_y_path) && file.exists(fit_a_path) && file.exists(fit_b_path)) {
  fit_y <- readRDS(fit_y_path)
  fit_A <- readRDS(fit_a_path)
  fit_B <- readRDS(fit_b_path)

  rmse_curve <- function(M, truth) {
    pm <- rowMeans(M); pm <- pm - mean(pm)
    tt <- truth - mean(truth); sqrt(mean((pm - tt)^2))
  }
  cov_curve <- function(M, truth, level = 0.95) {
    a <- (1 - level) / 2
    Mc <- sweep(M, 2, colMeans(M), FUN = "-")
    lo <- apply(Mc, 1, quantile, probs = a)
    hi <- apply(Mc, 1, quantile, probs = 1 - a)
    tt <- truth - mean(truth); mean(tt >= lo & tt <= hi)
  }

  cat("---- Plan C verdict (vs ortho, A, B) ----\n")
  cat(sprintf("  %-22s  %8s  %8s  %8s  %8s\n",
              "metric", "ortho", "A", "B", "C"))
  cat(paste0(rep("-", 64), collapse = ""), "\n")
  for (j in 1:3) {
    vy <- rmse_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vA <- rmse_curve(fit_A$f_main[[j]], sim$truth_f_grid[[j]])
    vB <- rmse_curve(fit_B$f_main[[j]], sim$truth_f_grid[[j]])
    vC <- rmse_curve(fit_C$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d (full)         %8.4f  %8.4f  %8.4f  %8.4f\n",
                j, vy, vA, vB, vC))
  }
  for (j in 1:3) {
    vy <- cov_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vA <- cov_curve(fit_A$f_main[[j]], sim$truth_f_grid[[j]])
    vB <- cov_curve(fit_B$f_main[[j]], sim$truth_f_grid[[j]])
    vC <- cov_curve(fit_C$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  cov_f%d                 %8.4f  %8.4f  %8.4f  %8.4f\n",
                j, vy, vA, vB, vC))
  }
  cat("\n")

  rmse_f1_C <- rmse_curve(fit_C$f_main[[1]], sim$truth_f_grid[[1]])
  cat(sprintf("Headline: RMSE_f1 (full) = %.4f\n", rmse_f1_C))
  cat(sprintf("  GLS prediction (same Sigma): 0.066\n"))
  cat(sprintf("  Plan A actual:               0.283\n"))
  cat(sprintf("  Plan B actual:               0.283\n\n"))

  cat("Decision rule:\n")
  cat("  - RMSE_f1 ~ 0.07-0.10 -> Sigma fluctuation CONFIRMED.\n")
  cat("                            MCMC averaging over Sigma_t was the bias source.\n")
  cat("  - RMSE_f1 ~ 0.15-0.25 -> Sigma fluctuation explains some, not all.\n")
  cat("  - RMSE_f1 >= 0.25     -> Sigma fluctuation FALSIFIED.\n")
  cat("                            Look at the spatial RE b's posterior.\n\n")
}

cat("Done.\n")
