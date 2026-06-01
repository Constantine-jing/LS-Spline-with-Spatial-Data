# ============================================================
# plot_seed1_comparison.R
#
# Three-panel comparison of f1 (Scenario B, seed 1):
#   Panel 1: Buggy ortho=TRUE  (RMSE_f1 = 0.330)
#   Panel 2: Post-fix          (RMSE_f1 = 0.073)
#   Panel 3: mgcv              (RMSE_f1 = 0.057)
#
# Each panel: truth (black), posterior mean / point estimate (red),
# 95% pointwise band (red shading).
# ============================================================

source("simulate_scenario_B.R")
sim <- simulate_scenario_B(seed = 1L)
xg <- sim$x_grid_1d
truth1 <- sim$truth_f_grid[[1]] - mean(sim$truth_f_grid[[1]])

# Load all three fits
fit_buggy <- readRDS("fit_scenarioB_seed1_ortho.rds")
fit_fixed <- readRDS("comparison/output/scenarioB_postfix/fit_postfix_seed1.rds")
fit_mgcv  <- readRDS("fit_scenarioB_seed1_mgcv_quick.rds")

# Helper: center curve, compute RMSE, get 95% band
center <- function(f) f - mean(f)
rmse <- function(pm, truth) sqrt(mean((pm - truth)^2))

# Buggy fit summary
pm_buggy <- center(rowMeans(fit_buggy$f_main[[1]]))
qq_buggy <- apply(sweep(fit_buggy$f_main[[1]], 2, colMeans(fit_buggy$f_main[[1]]), FUN = "-"),
                  1, quantile, probs = c(0.025, 0.975))

# Fixed fit summary
pm_fixed <- center(rowMeans(fit_fixed$f_main[[1]]))
qq_fixed <- apply(sweep(fit_fixed$f_main[[1]], 2, colMeans(fit_fixed$f_main[[1]]), FUN = "-"),
                  1, quantile, probs = c(0.025, 0.975))

# Refit mgcv fresh instead of using the saved object
library(mgcv)
df <- sim$data
fit_m_new <- gam(
  y ~ s(X1, bs = "ps", k = 20) +
    s(X2, bs = "ps", k = 20) +
    s(X3, bs = "ps", k = 20) +
    ti(X1, X2, bs = "ps", k = c(10, 10)) +
    ti(X1, X3, bs = "ps", k = c(10, 10)) +
    ti(X2, X3, bs = "ps", k = c(10, 10)) +
    s(lon, lat, bs = "gp"),
  data = df, method = "REML"
)

df_pred <- data.frame(X1 = xg,
                      X2 = mean(sim$data$X2),
                      X3 = mean(sim$data$X3),
                      lon = mean(sim$data$lon),
                      lat = mean(sim$data$lat))
pred_m <- predict(fit_m_new, newdata = df_pred, type = "terms", terms = "s(X1)")
f1_mgcv_pred <- as.numeric(pred_m) - mean(as.numeric(pred_m))

rmse_buggy <- rmse(pm_buggy, truth1)
rmse_fixed <- rmse(pm_fixed, truth1)
rmse_mgcv  <- rmse(f1_mgcv_pred, truth1)

# Range for consistent y-axis across panels
yrange <- range(c(truth1, qq_buggy, qq_fixed, f1_mgcv_pred))

# ---- Plot ----
png("seed1_comparison.png", width = 1500, height = 500, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# Panel 1: buggy
plot(xg, truth1, type = "l", lwd = 2, col = "black",
     ylim = yrange,
     main = sprintf("BUGGY ortho=TRUE\nRMSE_f1 = %.3f, cov = %.2f",
                    rmse_buggy,
                    mean(truth1 >= qq_buggy[1,] & truth1 <= qq_buggy[2,])),
     xlab = "X1", ylab = "f1(X1)")
polygon(c(xg, rev(xg)), c(qq_buggy[1,], rev(qq_buggy[2,])),
        col = adjustcolor("red", 0.20), border = NA)
lines(xg, pm_buggy, col = "red", lwd = 2.5)
lines(xg, truth1, col = "black", lwd = 2)
abline(v = c(0.10, 0.90), lty = 3, col = "gray60")
legend("bottomleft",
       c("truth", "posterior mean", "95% band"),
       col = c("black", "red", adjustcolor("red", 0.4)),
       lwd = c(2, 2.5, 6), lty = 1, cex = 0.85, bty = "n")

# Panel 2: fixed
plot(xg, truth1, type = "l", lwd = 2, col = "black",
     ylim = yrange,
     main = sprintf("POST-FIX (ortho=TRUE)\nRMSE_f1 = %.3f, cov = %.2f",
                    rmse_fixed,
                    mean(truth1 >= qq_fixed[1,] & truth1 <= qq_fixed[2,])),
     xlab = "X1", ylab = "f1(X1)")
polygon(c(xg, rev(xg)), c(qq_fixed[1,], rev(qq_fixed[2,])),
        col = adjustcolor("forestgreen", 0.20), border = NA)
lines(xg, pm_fixed, col = "forestgreen", lwd = 2.5)
lines(xg, truth1, col = "black", lwd = 2)
abline(v = c(0.10, 0.90), lty = 3, col = "gray60")

# Panel 3: mgcv
plot(xg, truth1, type = "l", lwd = 2, col = "black",
     ylim = yrange,
     main = sprintf("mgcv (REML)\nRMSE_f1 = %.3f",
                    rmse_mgcv),
     xlab = "X1", ylab = "f1(X1)")
polygon(c(xg, rev(xg)), c(mgcv_lo, rev(mgcv_hi)),
        col = adjustcolor("blue", 0.20), border = NA)
lines(xg, f1_mgcv_pred, col = "blue", lwd = 2.5)
lines(xg, truth1, col = "black", lwd = 2)
abline(v = c(0.10, 0.90), lty = 3, col = "gray60")

mtext("Scenario B, seed = 1: f1 recovery before/after off-by-one fix vs mgcv",
      outer = TRUE, cex = 1.1, font = 2)

dev.off()
cat("Saved: seed1_comparison.png\n")

# Also save a wider zoom-in version showing left boundary 0 to 0.2
png("seed1_comparison_zoom.png", width = 1500, height = 500, res = 120)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
zoom_idx <- which(xg <= 0.2)

for (lab in c("buggy", "fixed", "mgcv")) {
  if (lab == "buggy") {
    pm <- pm_buggy; lo <- qq_buggy[1,]; hi <- qq_buggy[2,]
    col_main <- "red"; title <- sprintf("BUGGY  RMSE=%.3f", rmse_buggy)
  } else if (lab == "fixed") {
    pm <- pm_fixed; lo <- qq_fixed[1,]; hi <- qq_fixed[2,]
    col_main <- "forestgreen"; title <- sprintf("POST-FIX  RMSE=%.3f", rmse_fixed)
  } else {
    pm <- f1_mgcv_pred; lo <- NULL; hi <- NULL
    col_main <- "blue"; title <- sprintf("mgcv  RMSE=%.3f", rmse_mgcv)
  }
  yr_zoom <- if (is.null(lo)) range(c(truth1[zoom_idx], pm[zoom_idx]))
  else range(c(truth1[zoom_idx], lo[zoom_idx], hi[zoom_idx], pm[zoom_idx]))
  plot(xg[zoom_idx], truth1[zoom_idx], type = "l", lwd = 2, col = "black",
       ylim = yr_zoom, main = title, xlab = "X1", ylab = "f1")
  if (!is.null(lo)) {
    polygon(c(xg[zoom_idx], rev(xg[zoom_idx])),
            c(lo[zoom_idx], rev(hi[zoom_idx])),
            col = adjustcolor(col_main, 0.20), border = NA)
  }
  lines(xg[zoom_idx], pm[zoom_idx], col = col_main, lwd = 2.5)
  lines(xg[zoom_idx], truth1[zoom_idx], col = "black", lwd = 2)
}
mtext("Left-boundary zoom: X1 in [0, 0.2]  (boundary 'spike' check)",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()
cat("Saved: seed1_comparison_zoom.png\n")

cat("\n--- Numerical summary ---\n")
cat(sprintf("BUGGY ortho:  RMSE_f1 = %.4f  cov = %.3f\n",
            rmse_buggy, mean(truth1 >= qq_buggy[1,] & truth1 <= qq_buggy[2,])))
cat(sprintf("POST-FIX:     RMSE_f1 = %.4f  cov = %.3f\n",
            rmse_fixed, mean(truth1 >= qq_fixed[1,] & truth1 <= qq_fixed[2,])))
cat(sprintf("mgcv:         RMSE_f1 = %.4f\n", rmse_mgcv))

cat(sprintf("\nFirst 5 grid points (showing the boundary 'spike' problem):\n"))
cat(sprintf("  %-6s %8s %8s %8s %8s\n", "x", "truth", "buggy", "fixed", "mgcv"))
for (i in 1:5) {
  cat(sprintf("  %.3f  %+8.3f %+8.3f %+8.3f %+8.3f\n",
              xg[i], truth1[i], pm_buggy[i], pm_fixed[i], f1_mgcv_pred[i]))
}
