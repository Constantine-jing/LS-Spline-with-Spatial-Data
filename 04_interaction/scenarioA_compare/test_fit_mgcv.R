# ============================================================
# test_fit_mgcv.R
# 
# End-to-end check: simulate Scenario A, fit mgcv,
# print quick diagnostics, save the canonical fit_result.
# 
# Expected runtime: ~5-30 seconds.
# 
# Run this AFTER test_fit_ours.R has succeeded; the simulator
# is identical, but it's good to confirm both wrappers see the
# same data when seed=1.
# ============================================================

source("simulate_scenario_A.R")
source("fit_mgcv_wrapper.R")

# ---- 1. simulate ----
sim <- simulate_scenario_A(seed = 1L)
cat(sprintf("[sim]   scenario A, seed=%d, n=%d, p=%d\n",
            sim$seed, sim$n, sim$p))

# ---- 2. fit mgcv ----
cat("\n[fit]   running mgcv...\n")
fit <- fit_mgcv(sim)

# ---- 3. timing ----
cat("\n--- timing ---\n")
cat(sprintf("  total_sec = %.2f\n", fit$timing$total_sec))
cat(sprintf("  fit_sec   = %.2f\n", fit$timing$fit_sec))
cat(sprintf("  post_sec  = %.2f\n", fit$timing$post_sec))

# ---- 4. variance components ----
cat("\n--- variance components ---\n")
cat(sprintf("  sigma2 (REML pt est)    : %.3f   (true=%.3f)\n",
            fit$var_comp$sigma2[1], sim$truth_params$sigma2))
cat(sprintf("  tau2_s (sig2/lambda_sp) : %.3f   (true=%.3f)\n",
            fit$var_comp$tau2_s[1], sim$truth_params$tau2_s))
for (j in seq_len(sim$p)) {
  cat(sprintf("  tau2_%d                  : %.4f\n",
              j, fit$var_comp$tau2[[j]][1]))
}
cat(sprintf("  rho_hat (gp smooth)     : %s\n",
            if (is.na(fit$settings$rho_hat)) "n/a (mgcv version-dependent)"
            else sprintf("%.4f   (true=%.4f)",
                         fit$settings$rho_hat, sim$truth_params$rho)))
cat(sprintf("  REML converged          : %s\n",
            fit$convergence$converged))

# ---- 5. recovery sanity ----
cat("\n--- main-effect recovery on x_grid_1d ---\n")
for (j in seq_len(sim$p)) {
  truth <- sim$truth_f_grid[[j]]
  draws <- fit$f_main[[j]]
  pmean <- rowMeans(draws)
  qlo   <- apply(draws, 1, quantile, probs = 0.025)
  qhi   <- apply(draws, 1, quantile, probs = 0.975)
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qlo & truth <= qhi)
  cat(sprintf("  f%d:  RMSE = %.3f   95%%coverage = %.2f\n",
              j, rmse, cov95))
}

cat("\n--- spatial RE recovery at observed locs ---\n")
s_pmean <- rowMeans(fit$s_obs)
s_rmse  <- sqrt(mean((s_pmean - sim$truth_s_obs)^2))
cat(sprintf("  s(loc):  RMSE = %.3f   cor(pmean, truth) = %.3f\n",
            s_rmse, cor(s_pmean, sim$truth_s_obs)))

# ---- 6. save ----
out_path <- file.path("..", "output", "fits",
                      sprintf("fit_mgcv_scenarioA_seed%d.rds", sim$seed))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, out_path)
cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))

# ---- 7. side-by-side comparison with ours, if available ----
ours_path <- file.path("..", "output", "fits",
                       sprintf("fit_ours_scenarioA_seed%d.rds", sim$seed))
if (file.exists(ours_path)) {
  cat("\n--- comparison with ours (seed=1) ---\n")
  fit_ours_loaded <- readRDS(ours_path)
  
  cat(sprintf("  Method  total_sec   RMSE(f1)  RMSE(f2)  RMSE(f3)  cov(f1) cov(f2) cov(f3)  RMSE(s)\n"))
  for (lab in c("ours", "mgcv")) {
    f <- if (lab == "ours") fit_ours_loaded else fit
    rmses <- numeric(sim$p); covs <- numeric(sim$p)
    for (j in seq_len(sim$p)) {
      pm <- rowMeans(f$f_main[[j]])
      qq <- apply(f$f_main[[j]], 1, quantile, probs = c(0.025, 0.975))
      rmses[j] <- sqrt(mean((pm - sim$truth_f_grid[[j]])^2))
      covs[j]  <- mean(sim$truth_f_grid[[j]] >= qq[1, ] &
                       sim$truth_f_grid[[j]] <= qq[2, ])
    }
    spm <- rowMeans(f$s_obs)
    s_rmse <- sqrt(mean((spm - sim$truth_s_obs)^2))
    cat(sprintf("  %-7s %8.1f   %.3f     %.3f     %.3f     %.2f    %.2f    %.2f     %.3f\n",
                lab, f$timing$total_sec, rmses[1], rmses[2], rmses[3],
                covs[1], covs[2], covs[3], s_rmse))
  }
}

# ---- 8. optional plot ----
if (interactive()) {
  par(mfrow = c(1, 3))
  for (j in seq_len(sim$p)) {
    truth <- sim$truth_f_grid[[j]]
    draws <- fit$f_main[[j]]
    pmean <- rowMeans(draws)
    qlo   <- apply(draws, 1, quantile, probs = 0.025)
    qhi   <- apply(draws, 1, quantile, probs = 0.975)
    yr <- range(c(truth, pmean, qlo, qhi))
    plot(sim$x_grid_1d, pmean, type = "l", lwd = 2, col = "darkgreen",
         ylim = yr, main = sprintf("f%d  (mgcv)", j),
         xlab = "x", ylab = "")
    polygon(c(sim$x_grid_1d, rev(sim$x_grid_1d)),
            c(qlo, rev(qhi)),
            col = adjustcolor("darkgreen", alpha.f = 0.2), border = NA)
    lines(sim$x_grid_1d, truth, lwd = 2, col = "red", lty = 2)
    abline(h = 0, lty = 3)
  }
  par(mfrow = c(1, 1))
}

cat("\nOK.\n")
