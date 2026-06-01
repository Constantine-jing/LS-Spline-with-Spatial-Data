# ============================================================
# run_scenarioB_postfix_orthoFALSE.R
#
# Post-fix Scenario B, seed 1, with orthogonalize = FALSE.
#
# Purpose: the off-by-one bug is now fixed. Every prior comparison
# of orthogonalize TRUE vs FALSE was on the buggy sampler and is
# not informative. This run produces the bug-free ortho=FALSE
# baseline so it can be compared against the existing post-fix
# ortho=TRUE fit (fit_postfix_seed1.rds, RMSE_f1 = 0.073).
#
# This tells us how to frame sec:int_design_anova:
#   - if ortho=FALSE ~ ortho=TRUE   -> orthogonalisation adopted
#                                       for identifiability cleanliness,
#                                       negligible effect on accuracy
#   - if ortho=FALSE notably worse  -> orthogonalisation materially
#                                       helps; can say so
#
# Settings are identical to run_scenarioB_postfix.R EXCEPT
# orthogonalize = FALSE.
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_postfix_orthoFALSE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" POST-FIX RUN: Scenario B seed 1, orthogonalize = FALSE\n")
cat(" (patched wrapper; for comparison vs the ortho=TRUE post-fix fit)\n")
cat("================================================================\n\n")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n\n", sim$n, sim$p))

t0 <- proc.time()
fit_orthoF <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize = FALSE,
    verbose       = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

out_rds <- file.path(out_dir, "fit_postfix_orthoFALSE_seed1.rds")
saveRDS(fit_orthoF, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Compare to the post-fix ortho=TRUE fit ----
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

orthoT_path <- "comparison/output/scenarioB_postfix/fit_postfix_seed1.rds"
have_T <- file.exists(orthoT_path)
if (have_T) fit_orthoT <- readRDS(orthoT_path)

cat("---- ortho=FALSE vs ortho=TRUE (both post-fix, seed 1) ----\n")
cat(sprintf("  %-22s  %14s  %14s\n", "metric", "ortho=FALSE", "ortho=TRUE"))
cat(paste0(rep("-", 56), collapse = ""), "\n")

for (j in 1:3) {
  vF <- rmse_curve(fit_orthoF$f_main[[j]], sim$truth_f_grid[[j]])
  vT <- if (have_T) rmse_curve(fit_orthoT$f_main[[j]], sim$truth_f_grid[[j]]) else NA
  cat(sprintf("  RMSE_f%d (full)         %14.4f  %14.4f\n", j, vF, vT))
}
for (j in 1:3) {
  vF <- cov_curve(fit_orthoF$f_main[[j]], sim$truth_f_grid[[j]])
  vT <- if (have_T) cov_curve(fit_orthoT$f_main[[j]], sim$truth_f_grid[[j]]) else NA
  cat(sprintf("  cov_f%d                 %14.4f  %14.4f\n", j, vF, vT))
}
for (key in names(fit_orthoF$f_int)) {
  vF <- rmse_curve(fit_orthoF$f_int[[key]], sim$truth_f_int[[key]])
  vT <- if (have_T && key %in% names(fit_orthoT$f_int))
          rmse_curve(fit_orthoT$f_int[[key]], sim$truth_f_int[[key]]) else NA
  cat(sprintf("  RMSE_f%-3s (full)       %14.4f  %14.4f\n", key, vF, vT))
}

# variance components (wrapper labels: var_comp$tau2_s is spatial sigma2,
# var_comp$sigma2 is nugget tau2 -- labels swapped in wrapper, values correct)
cat(sprintf("\nVariance components (ortho=FALSE):\n"))
cat(sprintf("  sigma2_spatial = %.4f  (truth 1.0)\n", mean(fit_orthoF$var_comp$tau2_s)))
cat(sprintf("  tau2_nugget    = %.4f  (truth 1.0)\n", mean(fit_orthoF$var_comp$sigma2)))
cat(sprintf("  rho            = %.4f  (truth 0.06)\n", mean(fit_orthoF$var_comp$rho)))

cat("\n---- Interpretation ----\n")
if (have_T) {
  rF1 <- rmse_curve(fit_orthoF$f_main[[1]], sim$truth_f_grid[[1]])
  rT1 <- rmse_curve(fit_orthoT$f_main[[1]], sim$truth_f_grid[[1]])
  delta <- rF1 - rT1
  cat(sprintf("  RMSE_f1: ortho=FALSE %.4f vs ortho=TRUE %.4f  (diff %+.4f)\n",
              rF1, rT1, delta))
  if (abs(delta) < 0.02) {
    cat("  -> Negligible difference. Frame sec:int_design_anova as a\n")
    cat("     principled identifiability property adopted for cleanliness,\n")
    cat("     not as an accuracy fix.\n")
  } else if (delta > 0.02) {
    cat("  -> ortho=TRUE is meaningfully better. Orthogonalisation\n")
    cat("     materially helps; sec:int_design_anova can say so.\n")
  } else {
    cat("  -> ortho=FALSE is better (unexpected). Investigate before\n")
    cat("     finalising sec:int_design_anova.\n")
  }
} else {
  cat("  ortho=TRUE post-fix fit not found; cannot compare.\n")
}

cat("\nDone.\n")
