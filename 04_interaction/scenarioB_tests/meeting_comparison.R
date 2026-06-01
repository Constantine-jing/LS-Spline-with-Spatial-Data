# ============================================================
# meeting_comparison.R
#
# Side-by-side comparison of the two seed-1 fits:
#   ortho=FALSE  (original behavior, basis-correlation present)
#   ortho=TRUE   (closed-form Lang-Brezger fix)
#
# Produces:
#   1. A consolidated metrics table (RMSE + coverage + variance components)
#   2. Per-smooth tau2_s posterior means (the smoothing-variance evidence)
#   3. Side-by-side plots of f1, f2, f3 truth vs both fits
#   4. Side-by-side plots of f12 truth vs both fits
# ============================================================

source("simulate_scenario_B.R")

sim <- simulate_scenario_B(seed = 1L)
fy  <- readRDS("fit_scenarioB_seed1_ortho.rds")
fn  <- readRDS("fit_scenarioB_seed1_noortho.rds")

# ============================================================
# (A) Consolidated metrics table
# ============================================================
rmse_curve <- function(draws_matrix, truth) {
  pm <- rowMeans(draws_matrix)
  pm <- pm - mean(pm)
  tt <- truth - mean(truth)
  sqrt(mean((pm - tt)^2))
}
cov_curve <- function(draws_matrix, truth, level = 0.95) {
  alpha <- (1 - level) / 2
  pm_c <- sweep(draws_matrix, 2, colMeans(draws_matrix), FUN = "-")
  lo <- apply(pm_c, 1, quantile, probs = alpha)
  hi <- apply(pm_c, 1, quantile, probs = 1 - alpha)
  tt <- truth - mean(truth)
  mean(tt >= lo & tt <= hi)
}

cat("\n================================================================\n")
cat(" SCENARIO B, seed=1: paired comparison\n")
cat(" ortho=FALSE (orig)  vs  ortho=TRUE (Lang-Brezger fix)\n")
cat("================================================================\n\n")

cat(sprintf("%-26s  %12s  %12s  %12s\n",
            "metric", "ortho=FALSE", "ortho=TRUE", "ratio (Y/N)"))
cat(paste0(rep("-", 70), collapse = ""), "\n")

print_row <- function(label, vn, vy) {
  ratio <- if (abs(vn) > 1e-12) vy / vn else NA_real_
  cat(sprintf("  %-24s  %12.4f  %12.4f  %12.3f\n", label, vn, vy, ratio))
}

# --- main effects ---
cat("\n  -- Main effect RMSE --\n")
for (j in 1:3) {
  vn <- rmse_curve(fn$f_main[[j]], sim$truth_f_grid[[j]])
  vy <- rmse_curve(fy$f_main[[j]], sim$truth_f_grid[[j]])
  print_row(sprintf("RMSE_f%d (full grid)", j), vn, vy)
}

# Interior-only (this is the fairer comparison)
interior_idx <- which(sim$x_grid_1d >= 0.10 & sim$x_grid_1d <= 0.90)
cat("\n  -- Main effect RMSE (interior x in [0.1, 0.9]) --\n")
for (j in 1:3) {
  pm_n <- rowMeans(fn$f_main[[j]]); pm_n <- pm_n - mean(pm_n)
  pm_y <- rowMeans(fy$f_main[[j]]); pm_y <- pm_y - mean(pm_y)
  tt   <- sim$truth_f_grid[[j]] - mean(sim$truth_f_grid[[j]])
  vn <- sqrt(mean((pm_n[interior_idx] - tt[interior_idx])^2))
  vy <- sqrt(mean((pm_y[interior_idx] - tt[interior_idx])^2))
  print_row(sprintf("RMSE_f%d (interior)", j), vn, vy)
}

# Coverage
cat("\n  -- Main effect 95% pointwise coverage --\n")
for (j in 1:3) {
  vn <- cov_curve(fn$f_main[[j]], sim$truth_f_grid[[j]])
  vy <- cov_curve(fy$f_main[[j]], sim$truth_f_grid[[j]])
  print_row(sprintf("cov_f%d", j), vn, vy)
}

# --- interaction surfaces ---
cat("\n  -- Interaction surface RMSE --\n")
for (k in names(fn$f_int)) {
  vn <- rmse_curve(fn$f_int[[k]], sim$truth_f_int[[k]])
  vy <- rmse_curve(fy$f_int[[k]], sim$truth_f_int[[k]])
  print_row(sprintf("RMSE_f%s", k), vn, vy)
}

cat("\n  -- Interaction surface coverage --\n")
for (k in names(fn$f_int)) {
  vn <- cov_curve(fn$f_int[[k]], sim$truth_f_int[[k]])
  vy <- cov_curve(fy$f_int[[k]], sim$truth_f_int[[k]])
  print_row(sprintf("cov_f%s", k), vn, vy)
}

# --- spatial RE ---
cat("\n  -- Spatial random effect --\n")
vn_rs <- sqrt(mean((rowMeans(fn$s_obs) - sim$truth_s_obs)^2))
vy_rs <- sqrt(mean((rowMeans(fy$s_obs) - sim$truth_s_obs)^2))
print_row("RMSE_s",  vn_rs, vy_rs)
vn_cs <- cor(rowMeans(fn$s_obs), sim$truth_s_obs)
vy_cs <- cor(rowMeans(fy$s_obs), sim$truth_s_obs)
print_row("cor_s",   vn_cs, vy_cs)

# --- variance components ---
cat("\n  -- Variance components (truth: sigma2=1.0, tau2_s=1.0, rho=0.06) --\n")
print_row("sigma2 (residual)",   mean(fn$var_comp$sigma2), mean(fy$var_comp$sigma2))
print_row("tau2_s (spatial)",    mean(fn$var_comp$tau2_s), mean(fy$var_comp$tau2_s))
print_row("rho",                 mean(fn$var_comp$rho),    mean(fy$var_comp$rho))

# --- THE KEY DIAGNOSTIC: per-smooth smoothing variances ---
cat("\n  -- Per-smooth tau2_s (smaller = more shrinkage = smoother fit) --\n")
for (j in 1:3) {
  print_row(sprintf("tau2_s_%d (main)", j),
            mean(fn$var_comp$tau2[[j]]),
            mean(fy$var_comp$tau2[[j]]))
}
for (k in names(fn$var_comp$tau2_int)) {
  print_row(sprintf("tau2_s_int_%s", k),
            mean(fn$var_comp$tau2_int[[k]]),
            mean(fy$var_comp$tau2_int[[k]]))
}

cat("\n  -- Timing --\n")
print_row("total_min",
          fn$timing$total_sec / 60,
          fy$timing$total_sec / 60)


# ============================================================
# (B) Plot: side-by-side main effects
# ============================================================
cat("\n\nPlotting main effects...\n")

png("compare_main_effects.png", width = 1400, height = 900, res = 110)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

xg <- sim$x_grid_1d

for (label_fit in c("noortho", "ortho")) {
  fit <- if (label_fit == "noortho") fn else fy
  title_prefix <- if (label_fit == "noortho") "ortho=FALSE" else "ortho=TRUE"
  for (j in 1:3) {
    pm <- rowMeans(fit$f_main[[j]])
    pm <- pm - mean(pm)
    qq <- apply(sweep(fit$f_main[[j]], 2, colMeans(fit$f_main[[j]]), FUN = "-"),
                1, quantile, probs = c(0.025, 0.975))
    truth_c <- sim$truth_f_grid[[j]] - mean(sim$truth_f_grid[[j]])
    rmse <- sqrt(mean((pm - truth_c)^2))

    plot(xg, truth_c, type = "l", lwd = 2, col = "black",
         ylim = range(c(truth_c, qq)),
         main = sprintf("%s | f%d (RMSE=%.3f)", title_prefix, j, rmse),
         xlab = sprintf("X%d", j), ylab = "f")
    polygon(c(xg, rev(xg)), c(qq[1,], rev(qq[2,])),
            col = adjustcolor("red", 0.2), border = NA)
    lines(xg, pm, col = "red", lwd = 2)
    lines(xg, truth_c, col = "black", lwd = 2)
    legend("bottomleft", c("truth", "post.mean", "95% CI"),
           col = c("black", "red", adjustcolor("red", 0.4)),
           lwd = c(2, 2, 6), cex = 0.7, bty = "n")
  }
}
dev.off()
cat("Saved: compare_main_effects.png\n")


# ============================================================
# (C) Plot: side-by-side interaction surfaces (just the nonzero pair 1_2)
# ============================================================
cat("\nPlotting interaction surfaces...\n")

png("compare_interaction_f12.png", width = 1400, height = 500, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 5))

u_grid <- sim$x_grid_2d$u
v_grid <- sim$x_grid_2d$v

# Truth
truth_mat <- matrix(sim$truth_f_int[["1_2"]], length(u_grid), length(v_grid))
zlim <- range(c(truth_mat,
                rowMeans(fn$f_int[["1_2"]]),
                rowMeans(fy$f_int[["1_2"]])))

image(u_grid, v_grid, truth_mat, zlim = zlim,
      col = hcl.colors(50, "Blue-Red"),
      main = "Truth f12", xlab = "X1", ylab = "X2")

pmean_n_mat <- matrix(rowMeans(fn$f_int[["1_2"]]), length(u_grid), length(v_grid))
rmse_n <- sqrt(mean((pmean_n_mat - truth_mat)^2))
image(u_grid, v_grid, pmean_n_mat, zlim = zlim,
      col = hcl.colors(50, "Blue-Red"),
      main = sprintf("ortho=FALSE (RMSE=%.3f)", rmse_n),
      xlab = "X1", ylab = "X2")

pmean_y_mat <- matrix(rowMeans(fy$f_int[["1_2"]]), length(u_grid), length(v_grid))
rmse_y <- sqrt(mean((pmean_y_mat - truth_mat)^2))
image(u_grid, v_grid, pmean_y_mat, zlim = zlim,
      col = hcl.colors(50, "Blue-Red"),
      main = sprintf("ortho=TRUE (RMSE=%.3f)", rmse_y),
      xlab = "X1", ylab = "X2")

dev.off()
cat("Saved: compare_interaction_f12.png\n")


# ============================================================
# (D) Plot: tau2_s_j posterior densities (the undersmoothing evidence)
# ============================================================
cat("\nPlotting per-smooth tau2_s densities...\n")

png("compare_tau2_main.png", width = 1200, height = 400, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (j in 1:3) {
  dn <- density(fn$var_comp$tau2[[j]])
  dy <- density(fy$var_comp$tau2[[j]])
  xlim <- range(c(dn$x, dy$x))
  ylim <- c(0, max(c(dn$y, dy$y)))
  plot(dn, xlim = xlim, ylim = ylim,
       main = sprintf("tau2_s_%d posterior (smaller = smoother)", j),
       xlab = sprintf("tau2_s_%d", j), col = "blue", lwd = 2)
  lines(dy, col = "red", lwd = 2)
  legend("topright", c("ortho=FALSE", "ortho=TRUE"),
         col = c("blue", "red"), lwd = 2, cex = 0.8, bty = "n")
  abline(v = mean(fn$var_comp$tau2[[j]]), col = "blue", lty = 2)
  abline(v = mean(fy$var_comp$tau2[[j]]), col = "red",  lty = 2)
}
dev.off()
cat("Saved: compare_tau2_main.png\n")

cat("\n=== Done. Three PNG files saved in working directory. ===\n")
