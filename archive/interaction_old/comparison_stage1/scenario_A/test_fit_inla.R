# ============================================================
# test_fit_inla.R
# 
# End-to-end check: simulate Scenario A, fit INLA,
# print quick diagnostics, save the canonical fit_result.
# 
# Expected runtime: 30-180 seconds depending on mesh size and
# whether INLA's caches are warm.
# ============================================================

source("simulate_scenario_A.R")
source("fit_inla_wrapper.R")

# ---- 1. simulate ----
sim <- simulate_scenario_A(seed = 1L)
cat(sprintf("[sim]   scenario A, seed=%d, n=%d, p=%d\n",
            sim$seed, sim$n, sim$p))

# ---- 2. fit INLA ----
cat("\n[fit]   running INLA (mesh build + SPDE fit + 2000 draws)...\n")
fit <- fit_inla(sim)

# ---- 3. timing ----
cat("\n--- timing ---\n")
cat(sprintf("  total_sec  = %.1f\n", fit$timing$total_sec))
cat(sprintf("  fit_sec    = %.1f\n", fit$timing$fit_sec))
cat(sprintf("  post_sec   = %.1f\n", fit$timing$post_sec))
cat(sprintf("  mesh nodes = %d\n",   fit$convergence$n_mesh))
cat(sprintf("  mlik       = %.2f\n", fit$convergence$mlik))

# ---- 4. variance components ----
cat("\n--- variance components (posterior mean, 95%% CI) ---\n")
report_vc <- function(label, x, true_val = NA, note = "") {
  ci <- quantile(x, c(0.025, 0.975))
  if (is.na(true_val)) {
    cat(sprintf("  %-10s: %.3f  [%.3f, %.3f]   %s\n",
                label, mean(x), ci[1], ci[2], note))
  } else {
    inside <- (true_val >= ci[1]) && (true_val <= ci[2])
    cat(sprintf("  %-10s: %.3f  [%.3f, %.3f]   (true=%.3f, %s) %s\n",
                label, mean(x), ci[1], ci[2], true_val,
                if (inside) "in CI" else "OUT OF CI", note))
  }
}
report_vc("sigma2", fit$var_comp$sigma2, sim$truth_params$sigma2)
report_vc("tau2_s", fit$var_comp$tau2_s, sim$truth_params$tau2_s,
          "[INLA spatial sd^2; not directly comparable to ours -- different param.]")
report_vc("rho",    fit$var_comp$rho,    NA,
          "[INLA practical range, NOT our matern_cor rho]")
for (j in seq_len(sim$p)) {
  report_vc(sprintf("tau2_%d", j), fit$var_comp$tau2[[j]])
}

# ---- 5. recovery ----
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
                      sprintf("fit_inla_scenarioA_seed%d.rds", sim$seed))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, out_path)
cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))

# ---- 7. three-way comparison ----
ours_path <- file.path("..", "output", "fits",
                       sprintf("fit_ours_scenarioA_seed%d.rds", sim$seed))
mgcv_path <- file.path("..", "output", "fits",
                       sprintf("fit_mgcv_scenarioA_seed%d.rds", sim$seed))

avail <- list()
if (file.exists(ours_path)) avail[["ours"]] <- readRDS(ours_path)
if (file.exists(mgcv_path)) avail[["mgcv"]] <- readRDS(mgcv_path)
avail[["inla"]] <- fit

if (length(avail) >= 2) {
  cat("\n--- comparison (seed=1) ---\n")
  cat(sprintf("  Method  total_sec   RMSE(f1)  RMSE(f2)  RMSE(f3)  cov(f1) cov(f2) cov(f3)  RMSE(s)\n"))
  for (lab in names(avail)) {
    f <- avail[[lab]]
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

cat("\nOK.\n")
