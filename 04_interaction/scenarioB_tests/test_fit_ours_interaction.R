# ============================================================
# test_fit_ours_interaction.R
# 
# End-to-end check: simulate Scenario B, fit ours-with-interactions,
# print quick diagnostics, save the canonical fit_result.
# 
# Expected runtime: ~40-90 minutes. Same per-iter cost as the
# Scenario A fit, just with more parameters (3 interaction blocks
# of 19x19 = 361 coefficients each).
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

# ---- 1. simulate ----
sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("[sim]   scenario B, seed=%d, n=%d, p=%d, c_int=%g\n",
            sim$seed, sim$n, sim$p, sim$truth_params$c_int))

# ---- 2. fit ----
cat("\n[fit]   running interaction Gibbs (this will take a while)...\n")
fit <- fit_ours_interaction(sim)

# ---- 3. timing ----
cat("\n--- timing ---\n")
cat(sprintf("  total_sec = %.1f  (%.1f min)\n",
            fit$timing$total_sec, fit$timing$total_sec / 60))
cat(sprintf("  fit_sec   = %.1f\n", fit$timing$fit_sec))
cat(sprintf("  post_sec  = %.2f\n", fit$timing$post_sec))

# ---- 4. variance components ----
cat("\n--- variance components (posterior mean, 95%% CI) ---\n")
report_vc <- function(label, x, true_val = NA) {
  ci <- quantile(x, c(0.025, 0.975))
  if (is.na(true_val)) {
    cat(sprintf("  %-16s: %.4f  [%.4f, %.4f]\n",
                label, mean(x), ci[1], ci[2]))
  } else {
    inside <- (true_val >= ci[1]) && (true_val <= ci[2])
    cat(sprintf("  %-16s: %.4f  [%.4f, %.4f]   (true=%.3f, %s)\n",
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
cat("\n  Interaction smoothing variances (small => model thinks pair is null)\n")
for (k in names(fit$var_comp$tau2_int)) {
  report_vc(sprintf("tau2_int_%s", k), fit$var_comp$tau2_int[[k]])
}

# ---- 5. acceptance + ESS ----
cat("\n--- MH acceptance rates ---\n")
print(round(fit$convergence$mh_accept, 3))
cat("\n--- ESS (target > 400) ---\n")
print(round(fit$convergence$ess, 0))

# ---- 6. main-effect recovery (same as scenario A) ----
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

# ---- 7. interaction recovery ----
cat("\n--- interaction recovery on flat_grid_2d ---\n")
for (k in sim$int_keys) {
  truth <- sim$truth_f_int[[k]]
  draws <- fit$f_int[[k]]
  pmean <- rowMeans(draws)
  qlo   <- apply(draws, 1, quantile, probs = 0.025)
  qhi   <- apply(draws, 1, quantile, probs = 0.975)
  
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qlo & truth <= qhi)
  width <- mean(qhi - qlo)
  
  is_null <- isTRUE(all(truth == 0))
  tag <- if (is_null) "[NULL]" else "[TRUE]"
  cat(sprintf("  f_%-3s %s:  RMSE = %.3f   cov95 = %.2f   mean width = %.3f\n",
              k, tag, rmse, cov95, width))
}

cat("\n--- spatial RE recovery at observed locs ---\n")
s_pmean <- rowMeans(fit$s_obs)
s_rmse  <- sqrt(mean((s_pmean - sim$truth_s_obs)^2))
cat(sprintf("  s(loc):  RMSE = %.3f   cor(pmean, truth) = %.3f\n",
            s_rmse, cor(s_pmean, sim$truth_s_obs)))

# ---- 8. save ----
out_path <- file.path("..", "output", "fits",
                      sprintf("fit_ours_scenarioB_seed%d.rds", sim$seed))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, out_path)
cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))

# ---- 9. optional plots ----
if (interactive()) {
  par(mfrow = c(2, 2))
  
  # the true f12 surface
  ng <- length(sim$x_grid_2d$u)
  Z_truth <- matrix(sim$truth_f_int[["1_2"]], ng, ng)
  image(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_truth,
        col = hcl.colors(64, "Blue-Red 3"),
        zlim = range(c(Z_truth)),
        xlab = "X1", ylab = "X2", main = "f_12 truth")
  contour(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_truth, add = TRUE)
  
  # the posterior-mean f12
  Z_hat <- matrix(rowMeans(fit$f_int[["1_2"]]), ng, ng)
  image(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_hat,
        col = hcl.colors(64, "Blue-Red 3"),
        zlim = range(c(Z_truth)),
        xlab = "X1", ylab = "X2", main = "f_12 posterior mean")
  contour(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_hat, add = TRUE)
  
  # null surface 1_3 (should be ~0)
  Z_null13 <- matrix(rowMeans(fit$f_int[["1_3"]]), ng, ng)
  image(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_null13,
        col = hcl.colors(64, "Blue-Red 3"),
        zlim = range(c(Z_truth)),
        xlab = "X1", ylab = "X3", main = "f_13 posterior mean (null)")
  
  # null surface 2_3 (should be ~0)
  Z_null23 <- matrix(rowMeans(fit$f_int[["2_3"]]), ng, ng)
  image(sim$x_grid_2d$u, sim$x_grid_2d$v, Z_null23,
        col = hcl.colors(64, "Blue-Red 3"),
        zlim = range(c(Z_truth)),
        xlab = "X2", ylab = "X3", main = "f_23 posterior mean (null)")
  par(mfrow = c(1, 1))
}

cat("\nOK.\n")
