# ============================================================
# test_interaction_on_additive_data.R
# 
# Diagnostic: does fitting our interaction model on data WITHOUT
# a real interaction also degrade f1 recovery?
# 
# Logic:
#   - "absorption from a real f12" hypothesis: f1 is bad only
#     when there's a real f12 to compete with. Fitting interaction
#     model on additive data should give f1 RMSE ~0.10.
#   - "interaction columns hurt regardless" hypothesis: just
#     having interaction columns degrades f1. f1 RMSE will
#     stay around 0.30 even with no real f12 in truth.
# 
# Reference points:
#   - ours additive on Scenario A:                       f1 RMSE = 0.09
#   - ours interaction on Scenario B (real f12):         f1 RMSE = 0.337
#   - ours interaction on Scenario A (this test):        f1 RMSE = ?
# 
# Decision rule:
#   - if RMSE drops to ~0.10:  absorption confirmed.
#                              ANOVA constraints would help by
#                              orthogonalizing f12 vs f1.
#   - if RMSE stays ~0.25-0.35: not absorption per se. Issue is
#                              something else (conditioning,
#                              prior interaction, identifiability
#                              with spatial RE, etc.). ANOVA
#                              constraints might not be the fix.
# 
# Expected runtime: ~2 hours (same Gibbs as Scenario B fits).
# ============================================================

source("simulate_scenario_A.R")   # additive truth (no interaction in DGP)
source("simulate_scenario_B.R")   # for the truth-grid format
source("fit_ours_interaction_wrapper.R")

# ---- 1. simulate Scenario A data ----
sim_A <- simulate_scenario_A(seed = 1L)
cat(sprintf("[sim]   scenario A, seed=%d, n=%d, p=%d (NO interaction in truth)\n",
            sim_A$seed, sim_A$n, sim_A$p))

# The interaction wrapper expects sim$int_keys, sim$truth_f_int, sim$flat_grid_2d.
# Scenario A has truth_f_int = list() (empty) and no int_keys.
# We need to "upgrade" sim_A to look like a Scenario B object so the
# wrapper can fit interaction blocks. The data are unchanged; only the
# bookkeeping is added (with truth_f_int = all zero).
cat("\n[note]  upgrading sim_A to interaction-schema layout (truth f_int = 0)\n")
sim_A$int_keys <- c("1_2", "1_3", "2_3")
sim_A$truth_f_int <- list(
  "1_2" = numeric(nrow(sim_A$flat_grid_2d)),  # all zeros
  "1_3" = numeric(nrow(sim_A$flat_grid_2d)),
  "2_3" = numeric(nrow(sim_A$flat_grid_2d))
)

# ---- 2. fit ours WITH interactions on this no-interaction data ----
cat("\n[fit]   running interaction Gibbs on Scenario A data...\n")
fit_intA <- fit_ours_interaction(sim_A, settings = list(verbose = FALSE))

out_path <- file.path("..", "output", "fits",
                      "fit_ours_interaction_on_scenarioA_seed1.rds")
saveRDS(fit_intA, out_path)

cat(sprintf("\n[time]   total = %.1fs (%.1f min)\n",
            fit_intA$timing$total_sec, fit_intA$timing$total_sec / 60))

# ---- 3. variance components ----
cat("\n--- variance components ---\n")
report_vc <- function(label, x, true_val = NA) {
  ci <- quantile(x, c(0.025, 0.975))
  if (is.na(true_val)) {
    cat(sprintf("  %-16s: %.4f  [%.4f, %.4f]\n",
                label, mean(x), ci[1], ci[2]))
  } else {
    inside <- (true_val >= ci[1]) && (true_val <= ci[2])
    cat(sprintf("  %-16s: %.4f  [%.4f, %.4f]   (true=%.3f, %s)\n",
                label, mean(x), ci[1], ci[2], true_val,
                if (inside) "in CI" else "OUT"))
  }
}
report_vc("sigma2",  fit_intA$var_comp$sigma2,  sim_A$truth_params$sigma2)
report_vc("tau2_s",  fit_intA$var_comp$tau2_s,  sim_A$truth_params$tau2_s)
report_vc("rho",     fit_intA$var_comp$rho,     sim_A$truth_params$rho)
for (j in seq_len(sim_A$p)) {
  report_vc(sprintf("tau2_%d", j), fit_intA$var_comp$tau2[[j]])
}
cat("\n  Interaction smoothing variances (truth: ALL NULL)\n")
for (k in names(fit_intA$var_comp$tau2_int)) {
  report_vc(sprintf("tau2_int_%s", k), fit_intA$var_comp$tau2_int[[k]])
}

# ---- 4. main-effect recovery ----
cat("\n--- main-effect recovery ---\n")
for (j in seq_len(sim_A$p)) {
  truth <- sim_A$truth_f_grid[[j]]
  draws <- fit_intA$f_main[[j]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  cat(sprintf("  f%d:  RMSE = %.3f   cov95 = %.2f\n", j, rmse, cov95))
}

# ---- 5. interaction "recovery" (truth is zero, so RMSE measures how much
#       spurious mass the model put in the interactions) ----
cat("\n--- interaction recovery (truth is ZERO for all) ---\n")
for (k in sim_A$int_keys) {
  truth <- sim_A$truth_f_int[[k]]
  draws <- fit_intA$f_int[[k]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  width <- mean(qq[2, ] - qq[1, ])
  cat(sprintf("  f_%-3s [NULL]:  RMSE = %.3f   cov95 = %.2f   width = %.3f\n",
              k, rmse, cov95, width))
}

# ---- 6. side-by-side diagnostic ----
cat("\n\n--- DIAGNOSTIC: f1 RMSE across configurations ---\n")
cat("  ours additive on Scenario A (no f12 in truth, no f12 in model): 0.091\n")
cat("  ours interaction on Scenario B (real f12 in truth, f12 in model): 0.337\n")
cat(sprintf("  ours interaction on Scenario A (no f12 in truth, f12 in model): %.3f\n",
            sqrt(mean((rowMeans(fit_intA$f_main[[1]]) - sim_A$truth_f_grid[[1]])^2))))

cat("\nReading the result:\n")
cat("  if ours-int-on-A is ~0.10  -> ABSORPTION confirmed.\n")
cat("                                f1 only suffers when real f12 is present.\n")
cat("                                ANOVA constraints could help.\n")
cat("  if ours-int-on-A is ~0.30  -> NOT ABSORPTION (or not just absorption).\n")
cat("                                Having f12 columns at all degrades f1.\n")
cat("                                ANOVA constraints might not fix it.\n")
cat("  if intermediate (~0.20)    -> mixed. Some absorption, some basis pollution.\n")

cat("\nOK.\n")
