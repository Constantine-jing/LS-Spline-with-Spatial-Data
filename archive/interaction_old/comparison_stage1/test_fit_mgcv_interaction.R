# ============================================================
# test_fit_mgcv_interaction.R
# 
# Fit mgcv on Scenario B seed=1 and compare to our existing
# REF fit (fit_ours_scenarioB_seed1.rds). Goal: figure out
# whether ours is bad in absolute terms or just bad relative to
# the comparator.
# 
# Expected runtime: a few seconds to maybe a minute (mgcv with
# k_main=20 and three ti(20,20) terms is heavier than Scenario A).
# ============================================================

source("simulate_scenario_B.R")
source("fit_mgcv_interaction_wrapper.R")

sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("[sim]   scenario B, seed=%d, n=%d, p=%d, c_int=%g\n",
            sim$seed, sim$n, sim$p, sim$truth_params$c_int))

cat("\n[fit]   running mgcv with ti() interactions...\n")
fit <- fit_mgcv_interaction(sim)

cat(sprintf("\n[time] total = %.1fs   fit = %.1fs   post = %.2fs\n",
            fit$timing$total_sec, fit$timing$fit_sec, fit$timing$post_sec))
cat(sprintf("       converged = %s\n", fit$convergence$converged))

# --------- recovery summary, mgcv ---------
cat("\n--- mgcv main-effect recovery ---\n")
for (j in seq_len(sim$p)) {
  truth <- sim$truth_f_grid[[j]]
  draws <- fit$f_main[[j]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov   <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  cat(sprintf("  f%d:  RMSE = %.3f   cov95 = %.2f\n", j, rmse, cov))
}

cat("\n--- mgcv interaction recovery ---\n")
for (k in sim$int_keys) {
  truth <- sim$truth_f_int[[k]]
  draws <- fit$f_int[[k]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov   <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  width <- mean(qq[2, ] - qq[1, ])
  is_null <- isTRUE(all(truth == 0))
  tag <- if (is_null) "[NULL]" else "[TRUE]"
  cat(sprintf("  f_%-3s %s:  RMSE = %.3f   cov95 = %.2f   width = %.3f\n",
              k, tag, rmse, cov, width))
}

cat("\n--- mgcv spatial RE recovery ---\n")
spm    <- rowMeans(fit$s_obs)
s_rmse <- sqrt(mean((spm - sim$truth_s_obs)^2))
s_cor  <- cor(spm, sim$truth_s_obs)
cat(sprintf("  s(loc):  RMSE = %.3f   cor = %.3f\n", s_rmse, s_cor))

# ---- save and side-by-side comparison ----
out_path <- file.path("..", "output", "fits", "fit_mgcv_scenarioB_seed1.rds")
saveRDS(fit, out_path)
cat(sprintf("\n[save]  %s\n", out_path))

ours_path <- file.path("..", "output", "fits", "fit_ours_scenarioB_seed1.rds")
if (file.exists(ours_path)) {
  cat("\n--- side-by-side: scenario B, seed=1 ---\n")
  fit_ours <- readRDS(ours_path)
  
  cat("Method  total_s   f1 RMSE  f1 cov  f2 RMSE  f3 RMSE  f12 RMSE  f12 cov  f12 wid\n")
  for (lab in c("ours", "mgcv")) {
    f <- if (lab == "ours") fit_ours else fit
    
    rmses <- numeric(sim$p); covs <- numeric(sim$p)
    for (j in seq_len(sim$p)) {
      pm <- rowMeans(f$f_main[[j]])
      qq <- apply(f$f_main[[j]], 1, quantile, probs = c(0.025, 0.975))
      rmses[j] <- sqrt(mean((pm - sim$truth_f_grid[[j]])^2))
      covs[j]  <- mean(sim$truth_f_grid[[j]] >= qq[1, ] &
                       sim$truth_f_grid[[j]] <= qq[2, ])
    }
    pm12 <- rowMeans(f$f_int[["1_2"]])
    qq12 <- apply(f$f_int[["1_2"]], 1, quantile, probs = c(0.025, 0.975))
    rmse12 <- sqrt(mean((pm12 - sim$truth_f_int[["1_2"]])^2))
    cov12  <- mean(sim$truth_f_int[["1_2"]] >= qq12[1, ] &
                   sim$truth_f_int[["1_2"]] <= qq12[2, ])
    wid12  <- mean(qq12[2, ] - qq12[1, ])
    
    cat(sprintf("%-7s %7.1f   %.3f   %.2f    %.3f    %.3f    %.3f      %.2f     %.3f\n",
                lab, f$timing$total_sec, rmses[1], covs[1],
                rmses[2], rmses[3],
                rmse12, cov12, wid12))
  }
  
  # null interactions side-by-side
  cat("\n--- null interactions (smaller is better) ---\n")
  cat("Method   f_1_3 RMSE  f_1_3 wid    f_2_3 RMSE  f_2_3 wid\n")
  for (lab in c("ours", "mgcv")) {
    f <- if (lab == "ours") fit_ours else fit
    summary_pair <- function(k) {
      pm <- rowMeans(f$f_int[[k]])
      qq <- apply(f$f_int[[k]], 1, quantile, probs = c(0.025, 0.975))
      r  <- sqrt(mean((pm - sim$truth_f_int[[k]])^2))
      w  <- mean(qq[2, ] - qq[1, ])
      sprintf("%.3f      %.3f", r, w)
    }
    cat(sprintf("%-7s  %s    %s\n", lab,
                summary_pair("1_3"), summary_pair("2_3")))
  }
}

cat("\nOK.\n")
