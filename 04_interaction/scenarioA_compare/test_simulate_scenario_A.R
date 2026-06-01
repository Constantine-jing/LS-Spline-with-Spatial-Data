# ============================================================
# test_simulate_scenario_A.R
# 
# Quick sanity check before moving on to method wrappers.
# Run from the comparison/src/ directory after copying
# ls_basis.R, spatial_utils.R, etc. alongside (or by editing
# the source paths below to point at the dissertation
# project root).
# 
# What this script verifies:
#   1. simulator runs without error on seed=1
#   2. truth functions on the eval grid are centered to <1e-3
#   3. sample variances of components are roughly as
#      calibrated (sigma2~1, var(b)~1)
#   4. quick visual: plot truth f_j on the eval grid
# ============================================================

source("simulate_scenario_A.R")

sim <- simulate_scenario_A(seed = 1L)

cat("\n--- structure ---\n")
cat("scenario       :", sim$scenario, "\n")
cat("seed           :", sim$seed, "\n")
cat("n              :", sim$n, "\n")
cat("p              :", sim$p, "\n")
cat("data dim       :", paste(dim(sim$data), collapse = " x "), "\n")
cat("data colnames  :", paste(colnames(sim$data), collapse = ", "), "\n")
cat("x_grid_1d len  :", length(sim$x_grid_1d), "\n")
cat("flat_grid_2d   :", paste(dim(sim$flat_grid_2d), collapse = " x "), "\n")

cat("\n--- truth function centering on x_grid_1d ---\n")
for (j in seq_along(sim$truth_f_grid)) {
  m <- mean(sim$truth_f_grid[[j]])
  s <- sd(sim$truth_f_grid[[j]])
  cat(sprintf("  f%d: mean = %+.4f, sd = %.4f\n", j, m, s))
}

cat("\n--- sampled variance components ---\n")
cat(sprintf("  var(b)         = %.3f   (target tau2_s = %.2f)\n",
            var(sim$truth_s_obs), sim$truth_params$tau2_s))
# residual variance check: y - eta_main - b should have variance ~ sigma2
eta_main <- numeric(sim$n)
for (j in seq_len(sim$p)) {
  fj <- scenarioA_truth_f_list()[[j]]
  eta_main <- eta_main + fj(sim$data[[paste0("X", j)]])
}
resid <- sim$data$y - eta_main - sim$truth_s_obs
cat(sprintf("  var(residual)  = %.3f   (target sigma2 = %.2f)\n",
            var(resid), sim$truth_params$sigma2))

cat("\n--- truth params ---\n")
print(sim$truth_params)

# optional plot
if (interactive()) {
  par(mfrow = c(1, 3))
  for (j in seq_along(sim$truth_f_grid)) {
    plot(sim$x_grid_1d, sim$truth_f_grid[[j]], type = "l",
         lwd = 2, col = "navy",
         main = sprintf("f%d (centered)", j),
         xlab = "x", ylab = sprintf("f%d(x)", j))
    abline(h = 0, lty = 3)
  }
  par(mfrow = c(1, 1))
}

cat("\nOK.\n")
