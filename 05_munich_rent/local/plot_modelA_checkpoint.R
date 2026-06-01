# ============================================================
# plot_modelA_checkpoint.R
# Run in a SEPARATE R session while Model B is still running.
# Loads Model A checkpoint and generates ALL 12 plots.
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("ls_interaction.R")

load("rent99.rda")
load("rent99_polys.rda")
load("munich_rent_full/Model_A_experts_excluded.rda")

outdir <- "munich_rent_plots"
dir.create(outdir, showWarnings = FALSE)
M_C <- 20

# ---- Rebuild data and basis ----
area_raw <- rent99$area; yearc_raw <- rent99$yearc
area_min <- min(area_raw); area_max <- max(area_raw)
yearc_min <- min(yearc_raw); yearc_max <- max(yearc_raw)
area_01 <- (area_raw - area_min)/(area_max - area_min)
yearc_01 <- (yearc_raw - yearc_min)/(yearc_max - yearc_min)
smooth_names <- c("area", "yearc")
y <- rent99$rentsqm

centroid_df <- data.frame(
  district = as.integer(names(rent99.polys)),
  cx = sapply(rent99.polys, function(p) mean(p[,1])),
  cy = sapply(rent99.polys, function(p) mean(p[,2])))
idx <- match(rent99$district, centroid_df$district)
if (any(is.na(idx))) {
  keep <- !is.na(idx); rent99 <- rent99[keep,]; y <- y[keep]
  area_01 <- area_01[keep]; yearc_01 <- yearc_01[keep]
  area_raw <- area_raw[keep]; yearc_raw <- yearc_raw[keep]
  idx <- idx[keep]
}
coords <- as.matrix(centroid_df[idx, c("cx","cy")]); n <- nrow(coords)
coord_min <- apply(coords,2,min); coord_max <- apply(coords,2,max)
coord_range <- coord_max - coord_min
coords_01 <- sweep(coords,2,coord_min)/matrix(coord_range,n,2,byrow=TRUE)
set.seed(42)
coords_01 <- coords_01 + matrix(rnorm(n*2,sd=0.001),n,2)
coords_01 <- pmin(pmax(coords_01,0),1)

obj_area  <- ls_build_one_full(area_01, M = M_C)
obj_yearc <- ls_build_one_full(yearc_01, M = M_C)
int_12 <- ls_build_interaction(obj_area, obj_yearc, orthogonalize = TRUE)

x_grid <- seq(0, 1, length.out = 101)
eta_samp <- gs$eta_samples
n_keep <- nrow(eta_samp)

bath_bin <- as.integer(rent99$bath=="premium")
kitchen_bin <- as.integer(rent99$kitchen=="premium")
cheat_bin <- as.integer(rent99$cheating=="yes")
V_lin <- cbind(bath_bin, kitchen_bin, cheat_bin)
W_area <- obj_area$W; W_yearc <- obj_yearc$W; W_int <- int_12$W_uv
H <- cbind(1, V_lin, W_area, W_yearc, W_int)


# ============================================================
# 1. Main effect curves with 95% credible bands
# ============================================================
pdf(file.path(outdir, "marginals_modelA.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
for (j in 1:2) {
  idx_j <- 1 + cm_main[[j]]
  if (j == 1) {
    W_g <- obj_area$design_new(x_grid, "W", TRUE)
    x_orig <- x_grid * (area_max - area_min) + area_min; xlab <- "Floor space (m^2)"
  } else {
    W_g <- obj_yearc$design_new(x_grid, "W", TRUE)
    x_orig <- x_grid * (yearc_max - yearc_min) + yearc_min; xlab <- "Year of construction"
  }
  curve_mat <- t(apply(eta_samp[, idx_j, drop = FALSE], 1,
                        function(b) as.vector(W_g %*% b)))
  f_hat <- colMeans(curve_mat)
  f_lo  <- apply(curve_mat, 2, quantile, 0.025)
  f_hi  <- apply(curve_mat, 2, quantile, 0.975)
  plot(x_orig, f_hat, type = "n", ylim = range(c(f_lo, f_hi)),
       main = paste0("f(", smooth_names[j], ")"), xlab = xlab, ylab = paste0("f_", j, "(x)"))
  polygon(c(x_orig, rev(x_orig)), c(f_lo, rev(f_hi)),
          col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
  lines(x_orig, f_hat, lwd = 2, col = "blue")
  abline(h = 0, lty = 3, col = "gray50")
  legend("topleft", legend = c("posterior mean", "95% band"),
         col = c("blue", rgb(0.2, 0.4, 0.8, 0.4)),
         lty = c(1, NA), lwd = c(2, 10), bty = "n", cex = 0.85)
}
dev.off()
cat("1. marginals_modelA.pdf\n")


# ============================================================
# 2. Main effect curves overlaid with raw data
# ============================================================
pdf(file.path(outdir, "marginals_with_data_modelA.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
y_centered <- y - mean(y)
for (j in 1:2) {
  idx_j <- 1 + cm_main[[j]]
  if (j == 1) {
    W_g <- obj_area$design_new(x_grid, "W", TRUE)
    x_orig <- x_grid * (area_max - area_min) + area_min
    x_data <- area_raw; xlab <- "Floor space (m^2)"
  } else {
    W_g <- obj_yearc$design_new(x_grid, "W", TRUE)
    x_orig <- x_grid * (yearc_max - yearc_min) + yearc_min
    x_data <- yearc_raw; xlab <- "Year of construction"
  }
  curve_mat <- t(apply(eta_samp[, idx_j, drop = FALSE], 1,
                        function(b) as.vector(W_g %*% b)))
  f_hat <- colMeans(curve_mat)
  f_lo  <- apply(curve_mat, 2, quantile, 0.025)
  f_hi  <- apply(curve_mat, 2, quantile, 0.975)
  plot(x_data, y_centered, pch = ".", col = rgb(0, 0, 0, 0.15),
       main = paste0("f(", smooth_names[j], ") with data"),
       xlab = xlab, ylab = "centered rentsqm",
       ylim = range(c(y_centered, f_lo, f_hi)))
  polygon(c(x_orig, rev(x_orig)), c(f_lo, rev(f_hi)),
          col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
  lines(x_orig, f_hat, lwd = 2, col = "blue")
  abline(h = 0, lty = 3, col = "gray50")
}
dev.off()
cat("2. marginals_with_data_modelA.pdf\n")


# ============================================================
# 3. Interaction surface — 3D persp
# ============================================================
grid2d <- expand.grid(area = x_grid, yearc = x_grid)
Wgu <- obj_area$design_new(grid2d$area, "W", TRUE)
Wgv <- obj_yearc$design_new(grid2d$yearc, "W", TRUE)
Wgi <- khatri_rao_rowwise_R(Wgu, Wgv)
if (int_12$recipe$orthogonalize && !is.null(int_12$recipe$A_uv))
  Wgi <- Wgi - cbind(1, Wgu, Wgv) %*% int_12$recipe$A_uv
idx_int <- 1 + cm_int[[1]]
beta_int_mean <- eta_post[idx_int]
f12_grid <- as.vector(Wgi %*% beta_int_mean)
f12_mat <- matrix(f12_grid, length(x_grid), length(x_grid))
area_orig  <- x_grid * (area_max - area_min) + area_min
yearc_orig <- x_grid * (yearc_max - yearc_min) + yearc_min

pdf(file.path(outdir, "interaction_persp_modelA.pdf"), width = 8, height = 7)
persp(area_orig, yearc_orig, f12_mat,
      theta = -40, phi = 25, expand = 0.6,
      xlab = "Floor space (m^2)", ylab = "Year of construction",
      zlab = "f_12(area, yearc)",
      main = "Interaction: floor space x year of construction",
      col = "lightblue", shade = 0.3, border = NA, ticktype = "detailed")
dev.off()
cat("3. interaction_persp_modelA.pdf\n")


# ============================================================
# 4. Interaction surface — multiple angles
# ============================================================
pdf(file.path(outdir, "interaction_angles_modelA.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))
angles <- list(c(-40, 25), c(-140, 25), c(-40, 45), c(30, 15))
titles <- c("View 1 (default)", "View 2 (rotated)", "View 3 (top-down)", "View 4 (side)")
for (i in 1:4) {
  persp(area_orig, yearc_orig, f12_mat,
        theta = angles[[i]][1], phi = angles[[i]][2], expand = 0.6,
        xlab = "Floor space", ylab = "Year", zlab = "f_12", main = titles[i],
        col = "lightblue", shade = 0.3, border = NA, ticktype = "detailed")
}
dev.off()
cat("4. interaction_angles_modelA.pdf\n")


# ============================================================
# 5. Interaction contour/heatmap
# ============================================================
pdf(file.path(outdir, "interaction_contour_modelA.pdf"), width = 8, height = 7)
filled.contour(area_orig, yearc_orig, f12_mat,
               xlab = "Floor space (m^2)", ylab = "Year of construction",
               main = "Interaction f_12(area, yearc)",
               color.palette = colorRampPalette(c("blue", "white", "red")),
               nlevels = 20)
dev.off()
cat("5. interaction_contour_modelA.pdf\n")


# ============================================================
# 6. Spatial effect map (posterior mean)
# ============================================================
b_by_district <- tapply(b_post, rent99$district, mean)

pdf(file.path(outdir, "spatial_mean_modelA.pdf"), width = 8, height = 8)
b_range <- range(b_by_district, na.rm = TRUE)
nc <- 100; pal <- colorRampPalette(c("blue","white","red"))(nc)
bc <- function(b) { i <- round((b-b_range[1])/diff(b_range)*(nc-1))+1; pal[pmin(pmax(i,1),nc)] }
plot(NULL, xlim = range(coords_01[,1]), ylim = range(coords_01[,2]),
     asp = 1, main = "Spatial effect b(s) — posterior mean", xlab = "x", ylab = "y")
for (d in names(rent99.polys)) {
  p <- rent99.polys[[d]]
  p01 <- sweep(p,2,coord_min)/matrix(coord_range,nrow(p),2,byrow=TRUE)
  v <- b_by_district[d]; if (!is.na(v)) polygon(p01[,1],p01[,2],col=bc(v),border="gray30",lwd=0.3)
}
lv <- seq(b_range[1], b_range[2], length.out = 5)
legend("bottomleft", legend = sprintf("%.2f", lv), fill = bc(lv), title = "b(s)", bty = "n", cex = 0.8)
dev.off()
cat("6. spatial_mean_modelA.pdf\n")


# ============================================================
# 7. Spatial uncertainty map (posterior sd)
# ============================================================
b_sd_obs <- apply(gs$b_samples, 2, sd)
b_sd_by_district <- tapply(b_sd_obs, rent99$district, mean)

pdf(file.path(outdir, "spatial_uncertainty_modelA.pdf"), width = 8, height = 8)
sd_range <- range(b_sd_by_district, na.rm = TRUE)
nc <- 100; pal_sd <- colorRampPalette(c("white","orange","darkred"))(nc)
bc_sd <- function(v) { i <- round((v-sd_range[1])/diff(sd_range)*(nc-1))+1; pal_sd[pmin(pmax(i,1),nc)] }
plot(NULL, xlim = range(coords_01[,1]), ylim = range(coords_01[,2]),
     asp = 1, main = "Spatial effect — posterior SD (uncertainty)", xlab = "x", ylab = "y")
for (d in names(rent99.polys)) {
  p <- rent99.polys[[d]]
  p01 <- sweep(p,2,coord_min)/matrix(coord_range,nrow(p),2,byrow=TRUE)
  v <- b_sd_by_district[d]; if (!is.na(v)) polygon(p01[,1],p01[,2],col=bc_sd(v),border="gray30",lwd=0.3)
}
lv <- seq(sd_range[1], sd_range[2], length.out = 5)
legend("bottomleft", legend = sprintf("%.3f", lv), fill = bc_sd(lv), title = "sd(b)", bty = "n", cex = 0.8)
dev.off()
cat("7. spatial_uncertainty_modelA.pdf\n")


# ============================================================
# 8. Trace plots
# ============================================================
pdf(file.path(outdir, "trace_modelA.pdf"), width = 10, height = 10)
par(mfrow = c(3, 2))
plot(gs$sigma2_samples, type = "l", main = expression(sigma^2), ylab = "", xlab = "iteration", col = "steelblue")
plot(density(gs$sigma2_samples), main = expression(paste("Density: ", sigma^2)), col = "steelblue", lwd = 2)
plot(gs$tau2_samples, type = "l", main = expression(tau^2), ylab = "", xlab = "iteration", col = "darkgreen")
plot(density(gs$tau2_samples), main = expression(paste("Density: ", tau^2)), col = "darkgreen", lwd = 2)
plot(gs$rho_samples, type = "l", main = expression(rho), ylab = "", xlab = "iteration", col = "darkorange")
plot(density(gs$rho_samples), main = expression(paste("Density: ", rho)), col = "darkorange", lwd = 2)
dev.off()
cat("8. trace_modelA.pdf\n")


# ============================================================
# 9. Smoothing variance traces
# ============================================================
pdf(file.path(outdir, "tau2s_trace_modelA.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))
plot(gs$tau2_s_main_samp[,1], type = "l", main = expression(paste(tau[s]^2, " (area)")), ylab = "", col = "purple")
plot(gs$tau2_s_main_samp[,2], type = "l", main = expression(paste(tau[s]^2, " (yearc)")), ylab = "", col = "purple")
plot(gs$tau2_s_int_samp[,1], type = "l", main = expression(paste(tau[s]^2, " (interaction)")), ylab = "", col = "darkred")
dev.off()
cat("9. tau2s_trace_modelA.pdf\n")


# ============================================================
# 10. Linear coefficient posterior densities
# ============================================================
pdf(file.path(outdir, "linear_coefs_modelA.pdf"), width = 10, height = 4)
par(mfrow = c(1, p_linear))
for (k in 1:p_linear) {
  samp_k <- gs$eta_samples[, 1 + k]
  ci <- quantile(samp_k, c(0.025, 0.975))
  plot(density(samp_k), main = linear_names[k], xlab = "Euro/m^2", col = "steelblue", lwd = 2)
  abline(v = 0, lty = 2, col = "red")
  abline(v = ci, lty = 3, col = "gray40")
  abline(v = mean(samp_k), lty = 1, col = "blue")
  legend("topright", legend = c(sprintf("mean=%.2f", mean(samp_k)),
         sprintf("95%%CI=[%.2f,%.2f]", ci[1], ci[2])), bty = "n", cex = 0.8)
}
dev.off()
cat("10. linear_coefs_modelA.pdf\n")


# ============================================================
# 11. Residual diagnostics
# ============================================================
y_hat <- as.vector(H %*% eta_post) + b_post
resid <- y - y_hat

pdf(file.path(outdir, "residuals_modelA.pdf"), width = 12, height = 8)
par(mfrow = c(2, 2))
plot(y_hat, y, pch = ".", col = rgb(0,0,0,0.2),
     main = "Fitted vs Observed", xlab = "Fitted", ylab = "Observed rentsqm")
abline(0, 1, col = "red", lwd = 2)
hist(resid, breaks = 60, main = "Residual histogram", xlab = "Residual",
     col = "lightblue", border = "white", freq = FALSE)
curve(dnorm(x, mean(resid), sd(resid)), add = TRUE, col = "red", lwd = 2)
plot(y_hat, resid, pch = ".", col = rgb(0,0,0,0.2),
     main = "Residuals vs Fitted", xlab = "Fitted", ylab = "Residual")
abline(h = 0, col = "red", lwd = 2)
lo <- loess(resid ~ y_hat, span = 0.3)
ord <- order(y_hat); lines(y_hat[ord], predict(lo)[ord], col = "blue", lwd = 2)
qqnorm(resid, pch = ".", col = rgb(0,0,0,0.3), main = "Normal Q-Q")
qqline(resid, col = "red", lwd = 2)
dev.off()
cat("11. residuals_modelA.pdf\n")


# ============================================================
# 12. Residual spatial map
# ============================================================
resid_by_district <- tapply(resid, rent99$district, mean)

pdf(file.path(outdir, "residual_map_modelA.pdf"), width = 8, height = 8)
r_abs <- max(abs(resid_by_district), na.rm = TRUE)
sym_range <- c(-r_abs, r_abs)
nc <- 100; pal_r <- colorRampPalette(c("blue","white","red"))(nc)
bc_r <- function(v) { i <- round((v-sym_range[1])/diff(sym_range)*(nc-1))+1; pal_r[pmin(pmax(i,1),nc)] }
plot(NULL, xlim = range(coords_01[,1]), ylim = range(coords_01[,2]),
     asp = 1, main = "Mean residuals by district", xlab = "x", ylab = "y")
for (d in names(rent99.polys)) {
  p <- rent99.polys[[d]]
  p01 <- sweep(p,2,coord_min)/matrix(coord_range,nrow(p),2,byrow=TRUE)
  v <- resid_by_district[d]; if (!is.na(v)) polygon(p01[,1],p01[,2],col=bc_r(v),border="gray30",lwd=0.3)
}
lv <- seq(sym_range[1], sym_range[2], length.out = 5)
legend("bottomleft", legend = sprintf("%.2f", lv), fill = bc_r(lv), title = "resid", bty = "n", cex = 0.8)
dev.off()
cat("12. residual_map_modelA.pdf\n")


# ============================================================
# Summary
# ============================================================
cat("\n=== Model A Summary (n=3082) ===\n")
cat(sprintf("  Runtime: %.1f hours\n", t_el/3600))
cat(sprintf("  sigma2: %.4f  tau2: %.4f  rho: %.4f\n",
            mean(gs$sigma2_samples), mean(gs$tau2_samples), mean(gs$rho_samples)))
cat(sprintf("  MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))
cat(sprintf("  tau2_s: area=%.5f  yearc=%.5f  int=%.5f\n",
            mean(gs$tau2_s_main_samp[,1]), mean(gs$tau2_s_main_samp[,2]),
            mean(gs$tau2_s_int_samp[,1])))
cat(sprintf("  Spatial b: [%.3f, %.3f]  sd=%.3f\n", min(b_post), max(b_post), sd(b_post)))
cat(sprintf("  Residual: mean=%.4f  sd=%.4f\n", mean(resid), sd(resid)))
cat("\n  Linear coefficients:\n")
for (k in 1:p_linear) {
  pm <- mean(gs$eta_samples[,1+k]); psd <- sd(gs$eta_samples[,1+k])
  ci <- quantile(gs$eta_samples[,1+k], c(0.025,0.975))
  cat(sprintf("    %s: %.3f (%.3f)  95%%CI=[%.3f, %.3f] %s\n",
              linear_names[k], pm, psd, ci[1], ci[2],
              ifelse(ci[1]>0 | ci[2]<0, "***", "")))
}
cat("\nAll 12 plots saved in munich_rent_full/\n")
