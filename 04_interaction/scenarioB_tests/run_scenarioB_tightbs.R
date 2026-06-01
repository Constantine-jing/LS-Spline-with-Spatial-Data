# ============================================================
# run_scenarioB_tightbs.R
#
# HYPOTHESIS 3 PRIOR-FIX RUN (candidate 1: tighten b_s)
# ------------------------------------------------------
# Refit Scenario B, seed 1, with the IG prior on each tau2_s_j
# tightened from b_smooth = 0.005 (default) to b_smooth = 0.0005.
#
# Mechanism:
#   The IG full conditional on tau2_s_j has scale parameter
#       b_smooth + (1/2) * beta_j' K_j beta_j.
#   The data term scales with f_j's amplitude and dominates the
#   prior anchor b_smooth = 0.005 for large-amplitude smooths
#   like f_1. Cutting b_smooth by 10x moves the prior anchor
#   small enough to NOT dominate, but lets the prior contribute
#   meaningfully when the data term is small (which it is in
#   stationary distribution if tau2_s_1 has been pulled down).
#   This is the cheapest test of Hypothesis 3 that does not
#   require refactoring the sampler.
#
# Caveats / what to watch:
#   - Risk: over-shrinks f_2 and f_3 (their amplitudes are smaller,
#     so even the cut prior could dominate). Watch RMSE_f2, RMSE_f3
#     for a regression vs ortho=TRUE baseline.
#   - The prior-sensitivity study (run_prior_sensitivity_sim4.R)
#     already showed the additive sampler is fairly robust to
#     b_smooth in [0.001, 0.01]; we are extrapolating slightly
#     below that range.
#
# Settings:
#   - same DGP (seed = 1)
#   - same MCMC budget
#   - orthogonalize = TRUE
#   - b_smooth = 0.0005 (everything else identical to ortho=TRUE)
#
# Output:
#   comparison/output/scenarioB_h3_tightbs/fit_h3_tightbs_seed1.rds
#
# Usage:
#   Rscript run_scenarioB_tightbs.R
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_h3_tightbs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" HYPOTHESIS 3 PRIOR FIX: tighten b_smooth from 0.005 to 0.0005\n")
cat("================================================================\n\n")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n", sim$n, sim$p))

b_s_new <- 0.0005
cat(sprintf("b_smooth: 0.005 (default) -> %.4f (this run)\n", b_s_new))
cat(sprintf("a_smooth: 1.0 (unchanged); applied to both main and interaction blocks\n\n"))

t0 <- proc.time()
fit_tb <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize = TRUE,
    b_smooth      = b_s_new,
    verbose       = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

out_rds <- file.path(out_dir, "fit_h3_tightbs_seed1.rds")
saveRDS(fit_tb, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Quick verdict ----
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
              "metric", "ortho=TRUE", "tight b_s", "ratio"))
  cat(paste0(rep("-", 60), collapse = ""), "\n")
  for (j in 1:3) {
    vy <- rmse_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vt <- rmse_curve(fit_tb$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d            %12.4f  %12.4f  %12.3f\n",
                j, vy, vt, vt / vy))
  }
  for (j in 1:3) {
    vy <- cov_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vt <- cov_curve(fit_tb$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  cov_f%d             %12.4f  %12.4f  %12.3f\n",
                j, vy, vt, vt / vy))
  }
  cat("\nReference: mgcv RMSE_f1 = 0.076 (seed 1).\n")
  cat("Decision rule:\n")
  cat("  - SUCCESS  : RMSE_f1 drops materially AND RMSE_f2, RMSE_f3\n")
  cat("               do not regress (>10% worse) vs ortho=TRUE.\n")
  cat("  - PARTIAL  : RMSE_f1 drops but f2/f3 over-shrunk -> need\n")
  cat("               half-Cauchy / PC prior on sigma_s instead.\n")
  cat("  - FAILURE  : RMSE_f1 unchanged -> H3 wrong, look elsewhere.\n\n")
} else {
  cat(sprintf("(skipping quick verdict: '%s' not found)\n",
              fit_y_path))
}

cat("Done. Run compare_h3_metrics.R for the full side-by-side table.\n")
