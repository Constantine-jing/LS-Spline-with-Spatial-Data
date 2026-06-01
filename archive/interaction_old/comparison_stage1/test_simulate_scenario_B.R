# ============================================================
# test_simulate_scenario_B.R
# 
# Sanity check on the Scenario B DGP before fitting anything.
# ============================================================

source("simulate_scenario_B.R")

sim <- simulate_scenario_B(seed = 1L)

cat("\n--- structure ---\n")
cat("scenario       :", sim$scenario, "\n")
cat("seed           :", sim$seed, "\n")
cat("n              :", sim$n, "\n")
cat("p              :", sim$p, "\n")
cat("c_int          :", sim$truth_params$c_int, "\n")
cat("data dim       :", paste(dim(sim$data), collapse = " x "), "\n")
cat("x_grid_1d len  :", length(sim$x_grid_1d), "\n")
cat("flat_grid_2d   :", paste(dim(sim$flat_grid_2d), collapse = " x "), "\n")
cat("int_keys       :", paste(sim$int_keys, collapse = ", "), "\n")

cat("\n--- truth function centering on x_grid_1d ---\n")
for (j in seq_along(sim$truth_f_grid)) {
  m <- mean(sim$truth_f_grid[[j]])
  s <- sd(sim$truth_f_grid[[j]])
  cat(sprintf("  f%d: mean = %+.4f, sd = %.4f\n", j, m, s))
}

cat("\n--- truth interaction surfaces on the 2D grid ---\n")
for (k in sim$int_keys) {
  fk <- sim$truth_f_int[[k]]
  cat(sprintf("  f_%s: mean = %+.4f, sd = %.4f, range = [%.3f, %.3f]\n",
              k, mean(fk), sd(fk), min(fk), max(fk)))
}

cat("\n--- sampled variance components ---\n")
cat(sprintf("  var(b)         = %.3f   (target tau2_s = %.2f)\n",
            var(sim$truth_s_obs), sim$truth_params$tau2_s))
# residual: y - eta_main - f12 - b
eta_main <- numeric(sim$n)
for (j in seq_len(sim$p)) {
  fj <- scenarioB_truth_f_list()[[j]]
  eta_main <- eta_main + fj(sim$data[[paste0("X", j)]])
}
resid <- sim$data$y - eta_main - sim$truth_f12_obs - sim$truth_s_obs
cat(sprintf("  var(residual)  = %.3f   (target sigma2 = %.2f)\n",
            var(resid), sim$truth_params$sigma2))
cat(sprintf("  var(f12 obs)   = %.3f   (signal magnitude of interaction)\n",
            var(sim$truth_f12_obs)))

cat("\n--- truth params ---\n")
print(sim$truth_params)

if (interactive()) {
  par(mfrow = c(2, 2))
  for (j in seq_along(sim$truth_f_grid)) {
    plot(sim$x_grid_1d, sim$truth_f_grid[[j]], type = "l",
         lwd = 2, col = "navy",
         main = sprintf("f%d (centered)", j),
         xlab = "x", ylab = sprintf("f%d(x)", j))
    abline(h = 0, lty = 3)
  }
  # 2D plot of f12 truth
  ng <- length(sim$x_grid_2d$u)
  Z <- matrix(sim$truth_f_int[["1_2"]], ng, ng)
  image(sim$x_grid_2d$u, sim$x_grid_2d$v, Z,
        col = hcl.colors(64, "Blue-Red 3"),
        xlab = "X1", ylab = "X2",
        main = "f_{12}(X1, X2)  truth")
  contour(sim$x_grid_2d$u, sim$x_grid_2d$v, Z, add = TRUE)
  par(mfrow = c(1, 1))
}

cat("\nOK.\n")
