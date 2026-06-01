# ============================================================
# run_scenarioB_clamp_all_h3.R
#
# PLAN A: Confirm/falsify smoothing-strength diagnosis
# -----------------------------------------------------
# Refit Scenario B, seed 1, with ALL THREE main-effect smoothing
# variances tau2_s_main[1], [2], [3] clamped at values that the
# noise-free GLS sweep predicted optimal:
#
#   tau2_s_main[1] = 0.00625  (sp ~ 200, GLS RMSE_f1 ~ 0.068)
#   tau2_s_main[2] = 0.00250  (sp ~ 500, GLS RMSE_f2 ~ 0.097)
#   tau2_s_main[3] = 0.00125  (sp ~ 1000, GLS RMSE_f3 ~ 0.091)
#
# These are derived as tau2_s = sigma2 / sp, with sigma2 = 1.25
# (current ortho=TRUE posterior mean of sigma2). Note sigma2 is
# still resampled in this run; only tau2_s_main is fixed.
#
# Logic:
#   * If smoothing-strength is the dominant cause of poor f1 fit,
#     RMSE_f1 in this run should drop to roughly 0.10-0.15
#     (matching the GLS prediction up to MCMC variability).
#   * If RMSE_f1 stays near 0.30, the smoothing-strength
#     diagnosis has a hole and we need to find what else is
#     driving MCMC to underperform GLS at the same sp.
#
# Output:
#   comparison/output/scenarioB_h3_clamp_all/fit_h3_clamp_all_seed1.rds
#
# Usage:
#   Rscript run_scenarioB_clamp_all_h3.R
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_h3_clamp_all")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" PLAN A: clamp ALL THREE main-effect smoothing variances\n")
cat(" at GLS-optimal values\n")
cat("================================================================\n\n")

# DGP: same seed-1 Scenario B
sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n", sim$n, sim$p))

# Clamp values from GLS sweep (see prior diagnostics)
clamp_values <- c(
  tau2_s_main_1 = 0.00625,   # f1: sp ~ 200
  tau2_s_main_2 = 0.00250,   # f2: sp ~ 500
  tau2_s_main_3 = 0.00125    # f3: sp ~ 1000
)
cat("Clamping main-effect smoothing variances at:\n")
for (j in seq_along(clamp_values)) {
  cat(sprintf("  tau2_s_main[%d] = %.5f  (target sp = %.0f)\n",
              j, clamp_values[j], 1.25 / clamp_values[j]))
}
cat("Interaction smoothing variances tau2_s_int[*] are NOT clamped.\n")
cat("(They still get conjugate IG draws as in production.)\n\n")

# Fit
t0 <- proc.time()
fit_clamp <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize     = TRUE,
    tau2_s_main_fixed = as.numeric(clamp_values),
    verbose           = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

# Save
out_rds <- file.path(out_dir, "fit_h3_clamp_all_seed1.rds")
saveRDS(fit_clamp, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Sanity check: clamp held ----
cat("---- Sanity: tau2_s_main posterior means in this run ----\n")
for (j in 1:3) {
  pm  <- mean(fit_clamp$var_comp$tau2[[j]])
  psd <- sd(fit_clamp$var_comp$tau2[[j]])
  cat(sprintf("  tau2_s_main[%d]  mean=%.5f  sd=%.5f  (clamp=%.5f)\n",
              j, pm, psd, clamp_values[j]))
}
cat("\n")

# ---- Quick verdict against ortho=TRUE baseline ----
fit_y_path <- "fit_scenarioB_seed1_ortho.rds"
if (file.exists(fit_y_path)) {
  fit_y <- readRDS(fit_y_path)

  rmse_curve <- function(M, truth) {
    pm <- rowMeans(M); pm <- pm - mean(pm)
    tt <- truth - mean(truth)
    sqrt(mean((pm - tt)^2))
  }
  cov_curve <- function(M, truth, level = 0.95) {
    a <- (1 - level) / 2
    Mc <- sweep(M, 2, colMeans(M), FUN = "-")
    lo <- apply(Mc, 1, quantile, probs = a)
    hi <- apply(Mc, 1, quantile, probs = 1 - a)
    tt <- truth - mean(truth)
    mean(tt >= lo & tt <= hi)
  }

  xg <- sim$x_grid_1d
  interior <- which(xg >= 0.10 & xg <= 0.90)
  rmse_curve_int <- function(M, truth) {
    pm <- rowMeans(M); pm <- pm - mean(pm)
    tt <- truth - mean(truth)
    sqrt(mean((pm[interior] - tt[interior])^2))
  }

  cat("---- Plan A verdict ----\n")
  cat(sprintf("  %-22s  %12s  %12s\n",
              "metric", "ortho=TRUE", "clamp_all"))
  cat(paste0(rep("-", 52), collapse = ""), "\n")
  for (j in 1:3) {
    vy <- rmse_curve    (fit_y    $f_main[[j]], sim$truth_f_grid[[j]])
    vc <- rmse_curve    (fit_clamp$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d (full)         %12.4f  %12.4f\n", j, vy, vc))
  }
  for (j in 1:3) {
    vy <- rmse_curve_int(fit_y    $f_main[[j]], sim$truth_f_grid[[j]])
    vc <- rmse_curve_int(fit_clamp$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d ([0.1,0.9])    %12.4f  %12.4f\n", j, vy, vc))
  }
  for (j in 1:3) {
    vy <- cov_curve(fit_y    $f_main[[j]], sim$truth_f_grid[[j]])
    vc <- cov_curve(fit_clamp$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  cov_f%d (full)          %12.4f  %12.4f\n", j, vy, vc))
  }

  cat("\nGLS predictions for comparison:\n")
  cat("  RMSE_f1 (full): GLS predicts ~0.068 at sp=200\n")
  cat("  RMSE_f2 (full): GLS predicts ~0.097 at sp=500\n")
  cat("  RMSE_f3 (full): GLS predicts ~0.091 at sp=1000\n")
  cat("  mgcv reference:  f1=0.057, f2~0.04 (edf=1.00), f3~0.05 (edf=1.89)\n\n")

  cat("Decision rule:\n")
  cat("  - If RMSE_f1 (full) drops to 0.07-0.15:\n")
  cat("       -> smoothing-strength CONFIRMED. Move to half-Cauchy MCMC.\n")
  cat("  - If RMSE_f1 (full) stays >0.25:\n")
  cat("       -> smoothing-strength insufficient. MCMC underperforms\n")
  cat("          GLS even at clamped sp. Need to investigate the gap.\n")
  cat("  - In between (0.15-0.25):\n")
  cat("       -> partial confirmation. Smoothing helps but other\n")
  cat("          factors also matter.\n\n")
} else {
  cat(sprintf("(skipping verdict: '%s' not found)\n", fit_y_path))
}

cat("Done.\n")
