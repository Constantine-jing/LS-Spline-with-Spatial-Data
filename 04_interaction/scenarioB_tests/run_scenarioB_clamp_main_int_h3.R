# ============================================================
# run_scenarioB_clamp_main_int_h3.R
#
# PLAN B: clamp interaction smoothing variances on top of Plan A
# --------------------------------------------------------------
# Refit Scenario B, seed 1, with:
#   - tau2_s_main[1,2,3] clamped at GLS-optimal (Plan A values)
#   - tau2_s_int[1_2, 1_3, 2_3] also clamped at small values to
#     pin down beta_int and remove its per-iteration sampling
#     variance as a coupling channel for beta_1.
#
# Rationale:
#   In Plan A, GLS at clamp values gave RMSE_f1 = 0.066.
#   MCMC at the same clamp values gave RMSE_f1 = 0.283.
#   Every MCMC draw of f1 had the same downward bump near
#   x ~ 0.05; not a sampling-variability artifact.
#
#   The hypothesis: the conditional posterior of beta_1 in MCMC
#   is biased by per-iteration coupling with beta_int. Each
#   iteration, beta_int gets a fresh draw conditioned on the
#   current Sigma; the H1+H2 orthogonalization was set up at
#   one fixed (REML) Sigma and doesn't hold under the time-
#   varying Sigma sampled by MCMC.
#
#   If we tightly clamp tau2_s_int as well, beta_int's draws are
#   tightly concentrated around their conditional mean. Under
#   strong shrinkage on beta_int, even with imperfect
#   orthogonality at any single iteration, the leakage into
#   beta_1 should be substantially reduced.
#
# Decision rule:
#   - RMSE_f1 (full) drops to ~0.10 or below
#       -> per-iteration interaction coupling CONFIRMED.
#          Fix is iterative (Lang-Brezger style) re-orthogonalization
#          of W_uv against W_main using the current MCMC Sigma at
#          each sweep.
#   - RMSE_f1 (full) stays > 0.25
#       -> interaction-coupling falsified. Look at sigma2 or rho
#          fluctuation as the coupling channel.
#   - In between
#       -> partial mechanism; need to test more.
#
# Settings:
#   - Same DGP (seed = 1)
#   - Same MCMC budget
#   - orthogonalize = TRUE
#   - tau2_s_main_fixed = c(0.00625, 0.00250, 0.00125)  (Plan A)
#   - tau2_s_int_fixed  = c(0.001, 0.001, 0.001)         (NEW)
#
# Note: clamping tau2_s_int *will* over-smooth f_1_2 (true surface
# has nontrivial structure). RMSE_f_1_2 will likely degrade in
# this run. That's a known trade-off — the test is whether
# pinning beta_int down rescues f_1.
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_h3_clamp_main_int")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat(" PLAN B: clamp main AND interaction smoothing variances\n")
cat("================================================================\n\n")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("Simulated Scenario B at seed=1: n=%d, p=%d\n", sim$n, sim$p))

clamp_main <- c(0.00625, 0.00250, 0.00125)   # same as Plan A
clamp_int  <- c(0.001, 0.001, 0.001)         # tight on all 3 interactions

cat("Clamp values:\n")
cat("  Main effects (same as Plan A):\n")
for (j in 1:3) {
  cat(sprintf("    tau2_s_main[%d] = %.5f  (target sp = %.0f)\n",
              j, clamp_main[j], 1.0 / clamp_main[j]))
}
cat("  Interaction blocks (NEW):\n")
for (k in 1:3) {
  cat(sprintf("    tau2_s_int[%d]  = %.5f  (tight)\n", k, clamp_int[k]))
}
cat("\nNote: tight interaction clamp WILL over-smooth f_1_2.\n")
cat("This is intentional. We are testing whether pinning beta_int\n")
cat("rescues beta_1, accepting f_1_2 degradation as the cost.\n\n")

t0 <- proc.time()
fit_B <- fit_ours_interaction(
  sim,
  settings = list(
    orthogonalize     = TRUE,
    tau2_s_main_fixed = clamp_main,
    tau2_s_int_fixed  = clamp_int,
    verbose           = TRUE
  )
)
elapsed_min <- (proc.time() - t0)[3] / 60
cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

out_rds <- file.path(out_dir, "fit_h3_clamp_main_int_seed1.rds")
saveRDS(fit_B, out_rds)
cat(sprintf("Saved: %s\n\n", out_rds))

# ---- Sanity check: clamps held ----
cat("---- Sanity: posterior means of clamped variances ----\n")
for (j in 1:3) {
  pm  <- mean(fit_B$var_comp$tau2[[j]])
  cat(sprintf("  tau2_s_main[%d]  mean=%.5f  (clamp=%.5f) %s\n",
              j, pm, clamp_main[j],
              if (abs(pm - clamp_main[j]) < 1e-10) "OK" else "FAIL"))
}
int_keys <- names(fit_B$var_comp$tau2_int)
for (k in seq_along(int_keys)) {
  pm <- mean(fit_B$var_comp$tau2_int[[k]])
  cat(sprintf("  tau2_s_int[%s]  mean=%.5f  (clamp=%.5f) %s\n",
              int_keys[k], pm, clamp_int[k],
              if (abs(pm - clamp_int[k]) < 1e-10) "OK" else "FAIL"))
}
cat("\n")

# ---- Verdict ----
fit_y_path <- "fit_scenarioB_seed1_ortho.rds"
fit_a_path <- "comparison/output/scenarioB_h3_clamp_all/fit_h3_clamp_all_seed1.rds"
if (file.exists(fit_y_path) && file.exists(fit_a_path)) {
  fit_y <- readRDS(fit_y_path)
  fit_A <- readRDS(fit_a_path)

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
  xg <- sim$x_grid_1d
  interior <- which(xg >= 0.10 & xg <= 0.90)
  rmse_int <- function(M, truth) {
    pm <- rowMeans(M); pm <- pm - mean(pm)
    tt <- truth - mean(truth); sqrt(mean((pm[interior] - tt[interior])^2))
  }

  cat("---- Plan B verdict (vs ortho=TRUE and Plan A) ----\n")
  cat(sprintf("  %-22s  %10s  %10s  %10s\n",
              "metric", "ortho", "Plan A", "Plan B"))
  cat(paste0(rep("-", 58), collapse = ""), "\n")
  for (j in 1:3) {
    vy <- rmse_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vA <- rmse_curve(fit_A$f_main[[j]], sim$truth_f_grid[[j]])
    vB <- rmse_curve(fit_B$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d (full)         %10.4f  %10.4f  %10.4f\n",
                j, vy, vA, vB))
  }
  for (j in 1:3) {
    vy <- rmse_int(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vA <- rmse_int(fit_A$f_main[[j]], sim$truth_f_grid[[j]])
    vB <- rmse_int(fit_B$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  RMSE_f%d ([0.1,0.9])    %10.4f  %10.4f  %10.4f\n",
                j, vy, vA, vB))
  }
  for (j in 1:3) {
    vy <- cov_curve(fit_y$f_main[[j]], sim$truth_f_grid[[j]])
    vA <- cov_curve(fit_A$f_main[[j]], sim$truth_f_grid[[j]])
    vB <- cov_curve(fit_B$f_main[[j]], sim$truth_f_grid[[j]])
    cat(sprintf("  cov_f%d  (full)         %10.4f  %10.4f  %10.4f\n",
                j, vy, vA, vB))
  }

  # interaction RMSEs (these will degrade in Plan B by design)
  for (k in seq_along(int_keys)) {
    key <- int_keys[k]
    vy <- rmse_curve(fit_y$f_int[[key]], sim$truth_f_int[[key]])
    vA <- rmse_curve(fit_A$f_int[[key]], sim$truth_f_int[[key]])
    vB <- rmse_curve(fit_B$f_int[[key]], sim$truth_f_int[[key]])
    cat(sprintf("  RMSE_f%-3s  (full)      %10.4f  %10.4f  %10.4f\n",
                key, vy, vA, vB))
  }
  cat("\n")

  rmse_f1_B <- rmse_curve(fit_B$f_main[[1]], sim$truth_f_grid[[1]])
  cat(sprintf("Headline: RMSE_f1 (full) = %.4f\n", rmse_f1_B))
  cat(sprintf("  GLS prediction (clamp_main only): 0.066\n"))
  cat(sprintf("  Plan A MCMC actual:               0.283\n\n"))
  cat("Decision rule:\n")
  cat("  - RMSE_f1 < 0.15  -> interaction-coupling CONFIRMED.\n")
  cat("                       Next: implement iterative\n")
  cat("                       re-orthogonalization in MCMC.\n")
  cat("  - 0.15-0.25       -> partial; coupling is a factor but\n")
  cat("                       not the only one.\n")
  cat("  - >= 0.25         -> interaction-coupling FALSIFIED.\n")
  cat("                       Look at sigma2 or rho time-variation.\n\n")
} else {
  cat("(skipping verdict: missing baseline files)\n")
}

cat("Done.\n")
