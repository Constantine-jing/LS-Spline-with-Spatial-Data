# ============================================================
# test_b_smooth_loose.R
# 
# Test the looser-smoothing-prior hypothesis on Scenario B.
# 
# Hypothesis: the f1 RMSE problem (ours 0.34 vs mgcv 0.08) is
# at least partly due to b_smooth=0.005 being too aggressive,
# shrinking f1 toward zero. The diagnostic showed ours' joint
# (f1+f12) RMSE was 0.40 vs mgcv 0.26 -- ours is *losing*
# signal in the X1 direction, not just shifting it.
# 
# Test: rerun with b_smooth = 0.05, same seed, same n=500, M=20.
# 
# Decision rule:
#   - If joint RMSE drops to ~0.26 or below: prior was dominant.
#     Use b_smooth=0.05 going forward, multi-seed Scenario B.
#   - If joint stays around 0.35+: basis issue too. Need ANOVA
#     constraints (Path B).
# 
# Expected runtime: ~2 hours.
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

# ---- 1. simulate ----
sim <- simulate_scenario_B(seed = 1L)
cat(sprintf("[sim]   scenario B, seed=%d, n=%d, p=%d, c_int=%g\n",
            sim$seed, sim$n, sim$p, sim$truth_params$c_int))

# ---- 2. fit ours with looser smoothing prior ----
cat("\n[fit]   running interaction Gibbs with b_smooth=0.05...\n")
cat("        (was 0.005; prior median tau2 ~0.07 vs ~0.007)\n")

fit_loose <- fit_ours_interaction(
  sim,
  settings = list(
    b_smooth = 0.05,    # was 0.005
    verbose  = FALSE
  )
)

# ---- 3. save and report ----
out_path <- file.path("..", "output", "fits",
                      "fit_ours_scenarioB_seed1_b005.rds")
saveRDS(fit_loose, out_path)
cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))

cat("\n--- timing ---\n")
cat(sprintf("  total = %.1fs (%.1f min)\n",
            fit_loose$timing$total_sec, fit_loose$timing$total_sec / 60))

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
report_vc("sigma2",  fit_loose$var_comp$sigma2,  sim$truth_params$sigma2)
report_vc("tau2_s",  fit_loose$var_comp$tau2_s,  sim$truth_params$tau2_s)
report_vc("rho",     fit_loose$var_comp$rho,     sim$truth_params$rho)
for (j in seq_len(sim$p)) {
  report_vc(sprintf("tau2_%d", j), fit_loose$var_comp$tau2[[j]])
}
for (k in names(fit_loose$var_comp$tau2_int)) {
  report_vc(sprintf("tau2_int_%s", k), fit_loose$var_comp$tau2_int[[k]])
}

cat("\n--- main-effect recovery ---\n")
for (j in seq_len(sim$p)) {
  truth <- sim$truth_f_grid[[j]]
  draws <- fit_loose$f_main[[j]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  cat(sprintf("  f%d:  RMSE = %.3f   cov95 = %.2f\n", j, rmse, cov95))
}

cat("\n--- interaction recovery ---\n")
for (k in sim$int_keys) {
  truth <- sim$truth_f_int[[k]]
  draws <- fit_loose$f_int[[k]]
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

# ---- 4. absorption diagnostic ----
cat("\n\n--- ABSORPTION DIAGNOSTIC: ours_loose vs mgcv (vs ours_tight ref) ---\n")

ours_tight_path <- file.path("..", "output", "fits",
                             "fit_ours_scenarioB_seed1.rds")
mgcv_path       <- file.path("..", "output", "fits",
                             "fit_mgcv_scenarioB_seed1.rds")
fit_tight <- readRDS(ours_tight_path)
fit_mgcv  <- readRDS(mgcv_path)

# build truth on the 2D grid
flat2d   <- sim$flat_grid_2d
u_vec    <- flat2d[, "u"]
v_vec    <- flat2d[, "v"]
.f1   <- function(x) 2.0 * (sin(pi * x) - 2 / pi)
.f2   <- function(x) 1.5 * (exp(x - 0.5) - (exp(0.5) - exp(-0.5)))
.f12  <- function(x1, x2) 1.5 * (sin(pi * x1) - 2 / pi) *
                                (exp(x2 - 0.5) - (exp(0.5) - exp(-0.5)))
truth_f1_2d   <- .f1(u_vec)
truth_f2_2d   <- .f2(v_vec)
truth_f12_2d  <- .f12(u_vec, v_vec)
truth_joint   <- truth_f1_2d + truth_f12_2d
truth_joint_v <- truth_f2_2d + truth_f12_2d

idx_u <- sapply(u_vec, function(u) which.min(abs(sim$x_grid_1d - u)))
idx_v <- sapply(v_vec, function(v) which.min(abs(sim$x_grid_1d - v)))

joint_summary <- function(fit, label) {
  f1_at_u <- fit$f_main[[1]][idx_u, , drop = FALSE]
  f2_at_v <- fit$f_main[[2]][idx_v, , drop = FALSE]
  f12     <- fit$f_int[["1_2"]]
  joint_u <- f1_at_u + f12
  joint_v <- f2_at_v + f12
  
  pm_f1   <- rowMeans(f1_at_u)
  pm_f12  <- rowMeans(f12)
  pm_ju   <- rowMeans(joint_u)
  pm_jv   <- rowMeans(joint_v)
  
  cat(sprintf(
    "  %-13s  f1(2D)=%.3f   f12=%.3f   joint(f1+f12)=%.3f   joint(f2+f12)=%.3f\n",
    label,
    sqrt(mean((pm_f1  - truth_f1_2d )^2)),
    sqrt(mean((pm_f12 - truth_f12_2d)^2)),
    sqrt(mean((pm_ju  - truth_joint )^2)),
    sqrt(mean((pm_jv  - truth_joint_v)^2))))
}

joint_summary(fit_tight, "ours (b=0.005)")
joint_summary(fit_loose, "ours (b=0.05) ")
joint_summary(fit_mgcv,  "mgcv         ")

cat("\nReading the joint(f1+f12):\n")
cat("  - ours_tight = 0.398 (was the problem)\n")
cat("  - mgcv       = 0.260 (the target)\n")
cat("  - ours_loose = (just printed above)\n")
cat("  - if ours_loose drops to ~0.30 or below: prior was the issue\n")
cat("  - if it stays around 0.35-0.40:           basis issue too\n")
