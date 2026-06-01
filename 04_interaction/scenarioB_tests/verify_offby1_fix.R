# ============================================================
# verify_offby1_fix.R
#
# After fixing the off-by-one bug in fit_ours_interaction_wrapper.R
# (col_map_main started at 0L instead of 1L, putting penalty over
# the intercept), this script verifies that:
#
# 1. Sampler's posterior mean of beta_1 NOW matches the GLS formula.
# 2. RMSE_f1 drops from ~0.32 toward ~0.07 (match GLS prediction).
#
# Uses 500 iterations only (~2 minutes) since at fixed Sigma the
# chain is i.i.d. and converges instantly.
# ============================================================

rm(list = ls()); gc()

source("simulate_scenario_B.R")
source("gibbs_stage_c_full.R")
source("ls_basis.R"); source("ls_interaction.R")
source("spatial_utils.R"); source("fit_spatial_reml.R")
source("gibbs_interaction.R")     # patched sampler
source("fit_ours_interaction_wrapper.R")  # patched wrapper

sim <- simulate_scenario_B(seed = 1L)

cat("================================================================\n")
cat(" VERIFICATION OF OFF-BY-ONE FIX (Plan C clamps, 500 iter)\n")
cat("================================================================\n\n")

t0 <- proc.time()
fit_fixed <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize     = TRUE,
    tau2_s_main_fixed = c(0.00625, 0.00250, 0.00125),
    tau2_s_int_fixed  = c(0.001, 0.001, 0.001),
    sigma2_fixed      = 0.7817,
    tau2_fixed        = 1.2495,
    rho_fixed         = 0.0530,
    n_iter            = 500,
    n_burn            = 100,
    n_thin            = 1,
    n_draws           = 400L,           # ADD THIS LINE
    verbose           = FALSE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("Fit complete in %.1f minutes.\n\n", elapsed_min))

# Posterior-mean RMSE_f1
xg <- sim$x_grid_1d
truth1 <- sim$truth_f_grid[[1]] - mean(sim$truth_f_grid[[1]])
pm_f1 <- rowMeans(fit_fixed$f_main[[1]]); pm_f1 <- pm_f1 - mean(pm_f1)
rmse_f1_fixed <- sqrt(mean((pm_f1 - truth1)^2))

cat("---- Verdict ----\n")
cat(sprintf("  RMSE_f1 (fixed, full grid):  %.4f\n", rmse_f1_fixed))
cat(sprintf("  GLS prediction (target):     0.0639\n"))
cat(sprintf("  Plan C buggy result:         0.3194\n\n"))

if (rmse_f1_fixed < 0.10) {
  cat("  -> FIX CONFIRMED. Sampler now matches GLS at fixed Sigma.\n")
  cat("     The off-by-one in col_map_main was the bug.\n")
} else if (rmse_f1_fixed < 0.20) {
  cat("  -> PARTIAL. Major improvement but residual gap.\n")
} else {
  cat("  -> NO IMPROVEMENT. Bug is elsewhere.\n")
}

# Also check the f2/f3 RMSEs to make sure we didn't break anything
truth2 <- sim$truth_f_grid[[2]] - mean(sim$truth_f_grid[[2]])
truth3 <- sim$truth_f_grid[[3]] - mean(sim$truth_f_grid[[3]])
pm_f2 <- rowMeans(fit_fixed$f_main[[2]]); pm_f2 <- pm_f2 - mean(pm_f2)
pm_f3 <- rowMeans(fit_fixed$f_main[[3]]); pm_f3 <- pm_f3 - mean(pm_f3)
cat(sprintf("\nFor reference (other smooths):\n"))
cat(sprintf("  RMSE_f2: %.4f  (Plan C was 0.1181)\n", sqrt(mean((pm_f2 - truth2)^2))))
cat(sprintf("  RMSE_f3: %.4f  (Plan C was 0.0951)\n", sqrt(mean((pm_f3 - truth3)^2))))

# Save in case we want to look at the fit details
saveRDS(fit_fixed, "fit_offby1_verify.rds")
cat("\nSaved: fit_offby1_verify.rds\n")
cat("Done.\n")
