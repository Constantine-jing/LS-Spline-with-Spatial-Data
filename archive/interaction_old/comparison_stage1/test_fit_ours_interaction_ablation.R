# ============================================================
# test_fit_ours_interaction_ablation.R
# 
# Diagnose the Scenario B problem (f1 RMSE = 0.337, coverage =
# 0.56 at the default n=500, M=20 setup) by running three small
# ablation experiments at seed=1 only:
# 
#   T1: n=500,  M=8,  MCMC=10500/2500/4   (was M=20)
#   T2: n=1000, M=8,  MCMC=10500/2500/4   (adds +n)
#   T3: n=1000, M=8,  MCMC=3000/1000/1    (matches validated work)
# 
# Comparison anchor: validated n=1000, M=8 result reports
#   f1 RMSE ~ 0.17, f12 RMSE ~ 0.058, f12 coverage ~ 0.94.
# 
# Decision rule for next steps:
#   - If T1 already gives f1 RMSE ~ 0.17 and f12 OK, the issue
#     was M=20. Stay at n=500.
#   - If T1 doesn't help but T2 does, n=1000 is required.
#   - If T2 doesn't help but T3 does (unlikely), MCMC matters.
# 
# Total wall time estimate: ~2 hours.
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

# helper to quickly summarize a fit on the truth
summarize_fit <- function(fit, sim, label) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  total_sec = %.1f  (%.1f min)\n",
              fit$timing$total_sec, fit$timing$total_sec / 60))
  
  # main-effect RMSE
  for (j in seq_len(sim$p)) {
    truth <- sim$truth_f_grid[[j]]
    draws <- fit$f_main[[j]]
    pmean <- rowMeans(draws)
    qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
    rmse  <- sqrt(mean((pmean - truth)^2))
    cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
    cat(sprintf("  f%d:    RMSE = %.3f   cov95 = %.2f\n",
                j, rmse, cov95))
  }
  # interaction RMSE
  for (k in sim$int_keys) {
    truth <- sim$truth_f_int[[k]]
    draws <- fit$f_int[[k]]
    pmean <- rowMeans(draws)
    qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
    rmse  <- sqrt(mean((pmean - truth)^2))
    cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
    width <- mean(qq[2, ] - qq[1, ])
    is_null <- isTRUE(all(truth == 0))
    tag <- if (is_null) "[NULL]" else "[TRUE]"
    cat(sprintf("  f_%-3s %s:  RMSE = %.3f   cov95 = %.2f   width = %.3f\n",
                k, tag, rmse, cov95, width))
  }
  cat(sprintf("  ESS tau2_int_1_2 = %d\n",
              round(fit$convergence$ess["tau2_int_1_2"])))
}

# -----------------------------------------------------------------
# Reference summary from existing seed=1 fit (n=500, M=20)
# -----------------------------------------------------------------
ref_path <- file.path("..", "output", "fits",
                      "fit_ours_scenarioB_seed1.rds")
if (file.exists(ref_path)) {
  cat("\n--- existing reference fit (n=500, M=20, MCMC=10500/2500/4) ---\n")
  ref_fit <- readRDS(ref_path)
  ref_sim <- simulate_scenario_B(seed = 1L, n = 500L)
  summarize_fit(ref_fit, ref_sim, "REF: n=500 M=20 MCMC=10500/2500/4")
} else {
  cat("\n[note] no existing reference fit found; running ablations only\n")
}

# -----------------------------------------------------------------
# T1: n=500, M=8
# -----------------------------------------------------------------
cat("\n\n##### Test 1: n=500, M=8, MCMC=10500/2500/4 #####\n")
sim_T1 <- simulate_scenario_B(seed = 1L, n = 500L)
fit_T1 <- fit_ours_interaction(
  sim_T1,
  settings = list(M = 8L, verbose = FALSE)
)
saveRDS(fit_T1, file.path("..", "output", "fits",
                          "fit_ours_scenarioB_seed1_T1_n500_M8.rds"))
summarize_fit(fit_T1, sim_T1, "T1: n=500 M=8")

# -----------------------------------------------------------------
# T2: n=1000, M=8
# -----------------------------------------------------------------
cat("\n\n##### Test 2: n=1000, M=8, MCMC=10500/2500/4 #####\n")
sim_T2 <- simulate_scenario_B(seed = 1L, n = 1000L)
fit_T2 <- fit_ours_interaction(
  sim_T2,
  settings = list(M = 8L, verbose = FALSE)
)
saveRDS(fit_T2, file.path("..", "output", "fits",
                          "fit_ours_scenarioB_seed1_T2_n1000_M8.rds"))
summarize_fit(fit_T2, sim_T2, "T2: n=1000 M=8")

# -----------------------------------------------------------------
# T3: n=1000, M=8, MCMC=3000/1000/1 (your validated config)
# -----------------------------------------------------------------
cat("\n\n##### Test 3: n=1000, M=8, MCMC=3000/1000/1 (validated cfg) #####\n")
sim_T3 <- simulate_scenario_B(seed = 1L, n = 1000L)
fit_T3 <- fit_ours_interaction(
  sim_T3,
  settings = list(
    M       = 8L,
    n_iter  = 3000L,
    n_burn  = 1000L,
    n_thin  = 1L,
    n_draws = 2000L,
    verbose = FALSE
  )
)
saveRDS(fit_T3, file.path("..", "output", "fits",
                          "fit_ours_scenarioB_seed1_T3_n1000_M8_validated.rds"))
summarize_fit(fit_T3, sim_T3, "T3: n=1000 M=8 validated MCMC")

cat("\n\nALL ABLATIONS DONE.\n")
cat("Compare RMSE(f1) and coverage to the reference (0.337, 0.56)\n")
cat("and the validated anchor (~0.17, 0.94).\n")
