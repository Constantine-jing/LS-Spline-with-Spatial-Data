# ============================================================
# diagnose_f12_absorption.R
# 
# Test the hypothesis: ours' apparent f12 advantage (RMSE 0.105
# vs mgcv 0.248) is actually mass-shifting between f1 and f12
# rather than better recovery of the true interaction.
# 
# Logic:
#   - The data only sees the sum. If ours is "swapping" f1 mass
#     into f12, the JOINT  surface  f1(u) + f12(u,v)  should
#     still be recovered well, while individual f1 and f12 are
#     wrong in compensating directions.
#   - If both ours and mgcv have similar JOINT RMSE, the
#     decompositions differ but the total signal is matched
#     equally well -- absorption hypothesis confirmed.
#   - If ours has clearly LOWER joint RMSE than mgcv, then ours
#     genuinely captures more signal; the f12 advantage is real.
# 
# Reads existing seed=1 fits, no new MCMC.
# ============================================================

source("simulate_scenario_B.R")

sim <- simulate_scenario_B(seed = 1L)

ours_path <- file.path("..", "output", "fits", "fit_ours_scenarioB_seed1.rds")
mgcv_path <- file.path("..", "output", "fits", "fit_mgcv_scenarioB_seed1.rds")

stopifnot(file.exists(ours_path), file.exists(mgcv_path))
fit_ours <- readRDS(ours_path)
fit_mgcv <- readRDS(mgcv_path)

# ---- Build truth_joint on the flattened 2D grid ----
# Each row of flat_grid_2d is (u, v). truth f1 depends only on u,
# truth f12 depends on (u, v). We need f1(u) + f12(u, v)
# evaluated at every (u, v) row of the grid.
flat2d <- sim$flat_grid_2d
u_vec  <- flat2d[, "u"]
v_vec  <- flat2d[, "v"]

# truth f1 on every grid row (using u coordinate only)
.scenarioB_f1   <- function(x) 2.0 * (sin(pi * x) - 2 / pi)
.scenarioB_f12  <- function(x1, x2) 1.5 * (sin(pi * x1) - 2 / pi) *
                                          (exp(x2 - 0.5) - (exp(0.5) - exp(-0.5)))

truth_f1_2d   <- .scenarioB_f1(u_vec)               # depends only on u
truth_f12_2d  <- .scenarioB_f12(u_vec, v_vec)
truth_joint   <- truth_f1_2d + truth_f12_2d         # length nrow(flat2d)

# Also test f2 + f12 (symmetric test, weaker because f2 amplitude is smaller)
.scenarioB_f2 <- function(x) 1.5 * (exp(x - 0.5) - (exp(0.5) - exp(-0.5)))
truth_f2_2d   <- .scenarioB_f2(v_vec)               # f2 depends on v
truth_joint_v <- truth_f2_2d + truth_f12_2d

# ---- For each method, build per-draw joint = f1_at_u + f12 ----
# fit$f_main[[1]] is (n_grid_1d x n_draws) on x_grid_1d.
# We need its values at u_vec. We do this by interpolation:
# x_grid_1d is uniform 0..1 with 101 points, u_vec is 30 unique
# values from the 30x30 grid. Find the closest x_grid_1d index
# for each u value, take that row.
match_u_to_grid <- function(u_vec, x_grid_1d) {
  # return integer index into x_grid_1d
  sapply(u_vec, function(u) which.min(abs(x_grid_1d - u)))
}

idx_u <- match_u_to_grid(u_vec, sim$x_grid_1d)
idx_v <- match_u_to_grid(v_vec, sim$x_grid_1d)

joint_summary <- function(fit, sim, label) {
  # f1 draws at u-coordinates of the 2D grid
  f1_at_u <- fit$f_main[[1]][idx_u, , drop = FALSE]   # (n_grid_2d x n_draws)
  f12     <- fit$f_int[["1_2"]]                       # (n_grid_2d x n_draws)
  joint   <- f1_at_u + f12
  
  pmean   <- rowMeans(joint)
  rmse_joint <- sqrt(mean((pmean - truth_joint)^2))
  
  # Also pull individual RMSEs against truth, on the SAME 2D grid,
  # for clean side-by-side reporting.
  pm_f1 <- rowMeans(f1_at_u)
  pm_f12 <- rowMeans(f12)
  rmse_f1_2d  <- sqrt(mean((pm_f1  - truth_f1_2d )^2))
  rmse_f12_2d <- sqrt(mean((pm_f12 - truth_f12_2d)^2))
  
  # f2 + f12 symmetric test
  f2_at_v <- fit$f_main[[2]][idx_v, , drop = FALSE]
  joint_v <- f2_at_v + f12
  pm_jv   <- rowMeans(joint_v)
  rmse_joint_v <- sqrt(mean((pm_jv - truth_joint_v)^2))
  
  cat(sprintf(
    "  %-7s   f1 RMSE = %.3f   f12 RMSE = %.3f   |   f1+f12 joint = %.3f   f2+f12 joint = %.3f\n",
    label, rmse_f1_2d, rmse_f12_2d, rmse_joint, rmse_joint_v))
  
  list(rmse_f1_2d = rmse_f1_2d,
       rmse_f12_2d = rmse_f12_2d,
       rmse_joint = rmse_joint,
       rmse_joint_v = rmse_joint_v)
}

cat("\n--- absorption diagnostic on 30x30 2D grid ---\n")
cat("All RMSEs computed on the SAME 900-point grid for fair comparison.\n\n")

cat(sprintf("  %-7s   %-21s   %-21s   |   joint metrics\n", "method",
            "f1 (2D)", "f12 (true)"))
cat(strrep("-", 92), "\n", sep = "")

ours_summ <- joint_summary(fit_ours, sim, "ours")
mgcv_summ <- joint_summary(fit_mgcv, sim, "mgcv")

cat("\n--- interpretation ---\n")
delta_f1    <- ours_summ$rmse_f1_2d  - mgcv_summ$rmse_f1_2d
delta_f12   <- ours_summ$rmse_f12_2d - mgcv_summ$rmse_f12_2d
delta_joint <- ours_summ$rmse_joint  - mgcv_summ$rmse_joint

cat(sprintf("  ours - mgcv:  f1 = %+.3f   f12 = %+.3f   joint = %+.3f\n",
            delta_f1, delta_f12, delta_joint))

cat("\nReading the joint RMSE:\n")
cat("  - if joint is similar between methods (|delta| small),\n")
cat("    decomposition differs but total signal recovery is the same.\n")
cat("    => ours is shifting mass between f1 and f12, not\n")
cat("       genuinely capturing more signal (ABSORPTION).\n")
cat("  - if ours has clearly lower joint RMSE,\n")
cat("    ours genuinely captures more signal (REAL advantage).\n")
cat("  - if mgcv has clearly lower joint RMSE,\n")
cat("    ours is just worse overall.\n")
