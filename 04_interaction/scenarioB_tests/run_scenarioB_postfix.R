# ============================================================
# run_scenarioB_postfix.R
#
# Rerun Scenario B seed 1 with the patched wrapper (off-by-one
# bug fixed) and the production setting:
#   - orthogonalize = TRUE
#   - all priors normal (no clamps, no fixes)
#   - full MCMC budget
#
# Compares to:
#   - mgcv seed 1 (RMSE_f1 = 0.057)
#   - old buggy ortho=TRUE result (RMSE_f1 = 0.330)
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_postfix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" POST-FIX RUN: Scenario B seed 1, ortho=TRUE, normal MCMC\n")
cat(" (off-by-one bug in col_map_main fixed)\n")
cat("================================================================\n\n")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n\n", sim$n, sim$p))

t0 <- proc.time()
fit_postfix <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize = TRUE,
    verbose       = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

out_rds <- file.path(out_dir, "fit_postfix_seed1.rds")
saveRDS(fit_postfix, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Compare to mgcv and old buggy ortho fit ----
cat("---- Headline comparison (Scenario B seed 1) ----\n")
xg <- sim$x_grid_1d

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

# Load old buggy fit
buggy_path <- "fit_scenarioB_seed1_ortho.rds"
have_buggy <- file.exists(buggy_path)
if (have_buggy) fit_buggy <- readRDS(buggy_path)

# Load mgcv
mgcv_path <- "fit_scenarioB_seed1_mgcv_quick.rds"
have_mgcv <- file.exists(mgcv_path)

cat(sprintf("  %-22s  %12s  %12s  %12s\n",
            "metric", "buggy ortho", "POST-FIX", "mgcv"))
cat(paste0(rep("-", 64), collapse = ""), "\n")

# RMSE for f1, f2, f3
for (j in 1:3) {
  vbug <- if (have_buggy) rmse_curve(fit_buggy$f_main[[j]], sim$truth_f_grid[[j]]) else NA
  vfix <- rmse_curve(fit_postfix$f_main[[j]], sim$truth_f_grid[[j]])
  vmg <- if (have_mgcv && j == 1) {
    fm <- readRDS(mgcv_path)
    sqrt(mean((fm$f1_mgcv_on_grid - (fm$truth1))^2))
  } else NA
  cat(sprintf("  RMSE_f%d (full)         %12.4f  %12.4f  %12s\n",
              j, vbug, vfix,
              if (is.na(vmg)) "NA" else sprintf("%.4f", vmg)))
}

# Coverage
for (j in 1:3) {
  vbug <- if (have_buggy) cov_curve(fit_buggy$f_main[[j]], sim$truth_f_grid[[j]]) else NA
  vfix <- cov_curve(fit_postfix$f_main[[j]], sim$truth_f_grid[[j]])
  cat(sprintf("  cov_f%d (full)          %12.4f  %12.4f  %12s\n",
              j, vbug, vfix, "NA"))
}

# Interaction surfaces
int_keys <- names(fit_postfix$f_int)
for (k in seq_along(int_keys)) {
  key <- int_keys[k]
  vbug <- if (have_buggy && key %in% names(fit_buggy$f_int)) {
    rmse_curve(fit_buggy$f_int[[key]], sim$truth_f_int[[key]])
  } else NA
  vfix <- rmse_curve(fit_postfix$f_int[[key]], sim$truth_f_int[[key]])
  cat(sprintf("  RMSE_f%-3s (full)       %12.4f  %12.4f  %12s\n",
              key, vbug, vfix, "NA"))
}

# Posterior mean of variance components (using CORRECT labels from sampler)
cat(sprintf("\nVariance components (post-fix):\n"))
cat(sprintf("  E[sigma2 | y] (spatial):  %.4f  (truth = 1.0)\n",
            mean(fit_postfix$var_comp$tau2_s)))   # NOTE: labels still swapped in wrapper
cat(sprintf("  E[tau2 | y]   (nugget) :  %.4f  (truth = 1.0)\n",
            mean(fit_postfix$var_comp$sigma2)))
cat(sprintf("  E[rho | y]:                %.4f  (truth = 0.06)\n",
            mean(fit_postfix$var_comp$rho)))
cat(sprintf("  E[tau2_s_main_1]:         %.5f\n", mean(fit_postfix$var_comp$tau2[[1]])))
cat(sprintf("  E[tau2_s_main_2]:         %.5f\n", mean(fit_postfix$var_comp$tau2[[2]])))
cat(sprintf("  E[tau2_s_main_3]:         %.5f\n", mean(fit_postfix$var_comp$tau2[[3]])))

cat("\n")
rmse_f1_postfix <- rmse_curve(fit_postfix$f_main[[1]], sim$truth_f_grid[[1]])
cat(sprintf("Headline: RMSE_f1 (post-fix) = %.4f\n", rmse_f1_postfix))
cat(sprintf("  Was (buggy):  0.330\n"))
cat(sprintf("  mgcv:         0.057\n"))
if (rmse_f1_postfix < 0.10) {
  cat("  -> FULLY SOLVED.\n")
} else if (rmse_f1_postfix < 0.20) {
  cat("  -> Major improvement; small residual gap (likely prior strength).\n")
} else {
  cat("  -> Surprising result; investigate further.\n")
}

cat("\nDone.\n")
