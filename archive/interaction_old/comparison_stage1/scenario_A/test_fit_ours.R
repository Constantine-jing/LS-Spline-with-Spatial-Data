# ============================================================
# test_fit_ours.R
# 
# End-to-end check: simulate Scenario A, fit our method,
# print quick diagnostics, save the canonical fit_result.
# 
# Run from the comparison/src/ folder (with all dissertation
# .R files copied alongside).
# 
# Time budget: with n=500, M=20, n_iter=10500, this should
# take ~30-90s on a modern laptop. If it's much slower, the
# bottleneck is the per-iteration Cholesky in the marginal
# likelihood -- expected, no action needed.
# ============================================================

source("simulate_scenario_A.R")
source("fit_ours_wrapper.R")

# ---- 1. simulate ----
sim <- simulate_scenario_A(seed = 1L)
cat(sprintf("[sim]   scenario A, seed=%d, n=%d, p=%d\n",
            sim$seed, sim$n, sim$p))

# ---- 2. fit ours ----
cat("\n[fit]   running gibbs (this may take a minute)...\n")
fit <- fit_ours(sim)

# ---- 3. quick diagnostics ----
cat("\n--- timing ---\n")
cat(sprintf("  total_sec = %.1f\n", fit$timing$total_sec))
cat(sprintf("  fit_sec   = %.1f\n", fit$timing$fit_sec))
cat(sprintf("  post_sec  = %.2f\n", fit$timing$post_sec))

cat("\n--- variance components (posterior mean, 95%% CI) ---\n")
report_vc <- function(label, x, true_val = NA) {
  ci <- quantile(x, c(0.025, 0.975))
  if (is.na(true_val)) {
    cat(sprintf("  %-10s: %.3f  [%.3f, %.3f]\n",
                label, mean(x), ci[1], ci[2]))
  } else {
    inside <- (true_val >= ci[1]) && (true_val <= ci[2])
    cat(sprintf("  %-10s: %.3f  [%.3f, %.3f]   (true=%.3f, %s)\n",
                label, mean(x), ci[1], ci[2], true_val,
                if (inside) "in CI" else "OUT OF CI"))
  }
}
report_vc("sigma2",  fit$var_comp$sigma2,  sim$truth_params$sigma2)
report_vc("tau2_s",  fit$var_comp$tau2_s,  sim$truth_params$tau2_s)
report_vc("rho",     fit$var_comp$rho,     sim$truth_params$rho)
for (j in seq_len(sim$p)) {
  report_vc(sprintf("tau2_%d", j), fit$var_comp$tau2[[j]])
}

cat("\n--- MH acceptance rates (target ~0.20-0.50) ---\n")
print(round(fit$convergence$mh_accept, 3))

cat("\n--- ESS (target > 400 on 2000 draws) ---\n")
print(round(fit$convergence$ess, 0))

# ---- 4. recovery sanity check (rough RMSE on the eval grid) ----
cat("\n--- main-effect recovery on x_grid_1d ---\n")
for (j in seq_len(sim$p)) {
  truth <- sim$truth_f_grid[[j]]
  draws <- fit$f_main[[j]]               # (n_grid x n_draws)
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

# ---- 5. save fit ----
out_path <- file.path("..", "output", "fits",
                      sprintf("fit_ours_scenarioA_seed%d.rds", sim$seed))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, out_path)
cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))

# ---- 6. optional plots ----
if (interactive()) {
  par(mfrow = c(1, 3))
  for (j in seq_len(sim$p)) {
    truth <- sim$truth_f_grid[[j]]
    draws <- fit$f_main[[j]]
    pmean <- rowMeans(draws)
    qlo   <- apply(draws, 1, quantile, probs = 0.025)
    qhi   <- apply(draws, 1, quantile, probs = 0.975)
    yr <- range(c(truth, pmean, qlo, qhi))
    plot(sim$x_grid_1d, pmean, type = "l", lwd = 2, col = "navy",
         ylim = yr, main = sprintf("f%d  (ours)", j),
         xlab = "x", ylab = "")
    polygon(c(sim$x_grid_1d, rev(sim$x_grid_1d)),
            c(qlo, rev(qhi)),
            col = adjustcolor("navy", alpha.f = 0.2), border = NA)
    lines(sim$x_grid_1d, truth, lwd = 2, col = "red", lty = 2)
    abline(h = 0, lty = 3)
    legend("topright", c("post. mean", "truth"),
           col = c("navy", "red"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  par(mfrow = c(1, 1))
}

cat("\nOK.\n")
