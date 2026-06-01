# ============================================================
# run_scenarioB_confirm_h3.R
#
# HYPOTHESIS 3 CONFIRMATION RUN
# -----------------------------
# Refit Scenario B, seed 1, with tau2_s_main[1] CLAMPED at 0.005
# (matching the posterior mean of tau2_s_main[2] from the existing
# fit_scenarioB_seed1_ortho.rds run, where f2 is well-recovered).
#
# Logic of the test:
#   * If undersmoothing on f1 is the cause of poor recovery, then
#     forcing tau2_s_1 to a "small enough" value (0.005) should
#     - drop RMSE_f1 from ~0.337 toward ~0.08 (mgcv level),
#     - eliminate the wiggly posterior mean,
#     - bring cov_f1 closer to 0.95.
#   * If clamping tau2_s_1 does NOT clean up f1, the diagnosis is
#     wrong and we look for a different mechanism (sampler mixing,
#     residual basis correlation against b, etc.) BEFORE tuning the
#     prior.
#
# Settings:
#   - same DGP (seed = 1)
#   - same MCMC budget (n_iter = 10500, n_burn = 2500, n_thin = 4)
#   - orthogonalize = TRUE  (we are testing on top of the H1 fix)
#   - tau2_s_2 and tau2_s_3 still drawn from IG (only f1 clamped)
#
# Output:
#   comparison/output/scenarioB_h3_confirm/fit_h3_confirm_seed1.rds
#
# Usage:
#   Rscript run_scenarioB_confirm_h3.R
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

# Output directory (kept separate from the existing seed-1 outputs
# so the original ortho=TRUE/FALSE artifacts are untouched).
out_dir <- file.path("comparison", "output", "scenarioB_h3_confirm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" HYPOTHESIS 3 CONFIRMATION: clamp tau2_s_main[1] = 0.005\n")
cat("================================================================\n\n")

# ---- DGP: Scenario B seed 1, same as the orthogonalized run ----
sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n", sim$n, sim$p))

# ---- Fit with f1 smoothing clamped, f2/f3 unconstrained ----
clamp_value <- 0.005
cat(sprintf("Clamping tau2_s_main[1] = %.4f (f2 and f3 left to IG sampler)\n\n",
            clamp_value))

t0 <- proc.time()
fit_h3 <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize     = TRUE,
    tau2_s_main_fixed = c(clamp_value, NA_real_, NA_real_),
    verbose           = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

# ---- Save full fit object ----
out_rds <- file.path(out_dir, "fit_h3_confirm_seed1.rds")
saveRDS(fit_h3, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Quick on-the-fly comparison vs the existing ortho=TRUE fit ----
# (Detailed side-by-side reporting is in compare_h3_metrics.R.)
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

  cat("---- Quick verdict (full grid) ----\n")
  cat(sprintf("  %-18s  %12s  %12s  %12s\n",
              "metric", "ortho=TRUE", "H3 clamp", "ratio"))
  cat(paste0(rep("-", 60), collapse = ""), "\n")
  for (j in 1:3) {
    vy <- rmse_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vh <- rmse_curve(fit_h3$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d            %12.4f  %12.4f  %12.3f\n",
                j, vy, vh, vh / vy))
  }
  for (j in 1:3) {
    vy <- cov_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vh <- cov_curve(fit_h3$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  cov_f%d             %12.4f  %12.4f  %12.3f\n",
                j, vy, vh, vh / vy))
  }
  cat("\nReference: mgcv RMSE_f1 = 0.076 (seed 1).\n")
  cat("Verdict rule: H3 confirmed if RMSE_f1 drops by >=50% AND cov_f1 rises.\n\n")
} else {
  cat(sprintf("(skipping quick verdict: '%s' not found)\n",
              fit_y_path))
}

cat("Done. Run compare_h3_metrics.R for the full side-by-side table.\n")
