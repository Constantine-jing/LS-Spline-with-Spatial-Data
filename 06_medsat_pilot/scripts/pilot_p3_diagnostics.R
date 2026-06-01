## ============================================================
## scripts/pilot_p3_diagnostics.R
##
## Full diagnostics for the n=1500 Hellbender all-pairs fit
## (output/fit_medsat_n1500_seed1_allpairs_hb.rds).
##
## Inputs:
##   - output/sim_medsat_london_asthma_n1500_seed1_v2.rds
##   - output/fit_medsat_n1500_seed1_allpairs_hb.rds
##   - output/sim_medsat_london_asthma_n500_seed1_v2.rds          (for side-by-side)
##   - output/fit_medsat_p1plus_allpairs_n500_seed1_local.rds     (for side-by-side)
##
## Outputs (output/plots/p3/):
##   - p3_smooths.pdf              : f_j(X_j) on raw NO2 / NDVI / IMD axes,
##                                   re-centered, 95% credible bands
##   - p3_spatial_re.pdf           : posterior mean spatial RE map at raw OSGB coords
##   - p3_interaction_surfaces.pdf : three fitted f_{u,v} interaction surfaces on
##                                   raw axes, with a hatch where the 95% CI
##                                   crosses zero
##   - p3_traces.pdf               : trace + running-mean panels for sigma2,
##                                   tau2_s (spatial), rho, tau2_s_main[1..3],
##                                   tau2_int[1_2|1_3|2_3] with ESS annotated
##   - p3_vs_p1plus_smooths.pdf    : side-by-side credible bands for n=500
##                                   p1plus vs n=1500 p3
##   - p3_vs_p1plus_varcomp.txt    : variance-component comparison table
##
## Reads only RDS files; does NOT refit anything.
## ============================================================

suppressPackageStartupMessages({
  library(coda)
})

OUT_DIR <- "output/plots/p3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FIT_P3_RDS    <- "output/fit_medsat_n1500_seed1_allpairs_hb.rds"
SIM_P3_RDS    <- "output/sim_medsat_london_asthma_n1500_seed1_v2.rds"
FIT_P1P_RDS   <- "output/fit_medsat_p1plus_allpairs_n500_seed1_local.rds"
SIM_P1P_RDS   <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"

stopifnot(file.exists(FIT_P3_RDS),  file.exists(SIM_P3_RDS),
          file.exists(FIT_P1P_RDS), file.exists(SIM_P1P_RDS))

cat("[load] reading RDS files...\n")
fit3 <- readRDS(FIT_P3_RDS)
sim3 <- readRDS(SIM_P3_RDS)
fit5 <- readRDS(FIT_P1P_RDS)
sim5 <- readRDS(SIM_P1P_RDS)

stopifnot(isTRUE(fit3$settings$fit_interactions))
# p1plus predates the `fit_interactions` setting being recorded; assert via
# f_int instead (must have all 3 pairwise interaction surfaces).
stopifnot(length(fit5$f_int) == 3L,
          setequal(names(fit5$f_int), c("1_2","1_3","2_3")))
stopifnot(sim3$n == 1500L, sim5$n == 500L)
stopifnot(setequal(sim3$int_keys, c("1_2","1_3","2_3")),
          setequal(sim5$int_keys, c("1_2","1_3","2_3")))

covnames <- c("NO2", "NDVI", "IMD")
int_pair_idx <- list("1_2" = c(1L, 2L),
                     "1_3" = c(1L, 3L),
                     "2_3" = c(2L, 3L))

# ---------- small helpers ----------
ess_one <- function(x) as.numeric(coda::effectiveSize(coda::as.mcmc(x)))

post_summary <- function(x, name) {
  q <- quantile(x, c(0.025, 0.5, 0.975))
  data.frame(name = name,
             mean = mean(x), sd = sd(x),
             q025 = q[1], q500 = q[2], q975 = q[3],
             ess  = ess_one(x),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Compute pointwise summaries of an (n_grid x n_draws) matrix.
post_pointwise <- function(M, probs = c(0.025, 0.975)) {
  list(mean = rowMeans(M),
       lo   = apply(M, 1, quantile, probs[1]),
       hi   = apply(M, 1, quantile, probs[2]))
}

# Re-centering anchor: posterior-mean curve evaluated at observed X_j (matches
# the metric in pilot_p1_diagnostics.R; centering-invariant for shape comparison
# and consistent with the way mgcv reports iterms).
recenter_smooth <- function(curve_mean, curve_lo, curve_hi, xs, obs_X) {
  anchor <- mean(stats::approx(x = xs, y = curve_mean, xout = obs_X, rule = 2)$y)
  list(mean = curve_mean - anchor,
       lo   = curve_lo   - anchor,
       hi   = curve_hi   - anchor,
       anchor = anchor)
}

# =============================================================
# 1. p3_smooths.pdf  -- f_j(X_j) on raw axes, re-centered, 95% CI
# =============================================================
cat("\n[plot 1/5] p3_smooths.pdf\n")
xs <- sim3$x_grid_1d
smooths3 <- vector("list", 3)
for (j in 1:3) {
  pw <- post_pointwise(fit3$f_main[[j]])
  rc <- recenter_smooth(pw$mean, pw$lo, pw$hi, xs,
                        sim3$data[[paste0("X", j)]])
  smooths3[[j]] <- rc
}

pdf(file.path(OUT_DIR, "p3_smooths.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), las = 1)
for (j in 1:3) {
  Xlo <- sim3$X_scale[[j]]$lo; Xhi <- sim3$X_scale[[j]]$hi
  xs_raw <- Xlo + xs * (Xhi - Xlo)
  rc <- smooths3[[j]]
  rng <- range(c(rc$lo, rc$hi))
  amp <- diff(range(rc$mean))
  plot(xs_raw, rc$mean, type = "n", ylim = rng,
       xlab = sprintf("%s (raw units)", covnames[j]),
       ylab = sprintf("centred f(%s)", covnames[j]),
       main = sprintf("%s   amp/sd(y) = %.2f", covnames[j],
                      amp / sd(sim3$data$y)))
  polygon(c(xs_raw, rev(xs_raw)), c(rc$lo, rev(rc$hi)),
          col = grDevices::adjustcolor("steelblue", 0.30), border = NA)
  lines(xs_raw, rc$mean, lwd = 2, col = "steelblue4")
  abline(h = 0, lty = 3, col = "grey50")
  # rug of observed X to show density
  rug(Xlo + sim3$data[[paste0("X", j)]] * (Xhi - Xlo),
      ticksize = 0.02, col = grDevices::adjustcolor("grey40", 0.4))
}
mtext("P3 (n=1500): main-effect smooths, posterior mean + 95% CI",
      side = 3, line = -1.4, outer = TRUE, cex = 1.0)
par(op); dev.off()

# =============================================================
# 2. p3_spatial_re.pdf  -- posterior mean s(loc_i) on raw OSGB coords
# =============================================================
cat("[plot 2/5] p3_spatial_re.pdf\n")
s_mean <- rowMeans(fit3$s_obs)
s_sd   <- apply(fit3$s_obs, 1, sd)

pdf(file.path(OUT_DIR, "p3_spatial_re.pdf"), width = 12, height = 6)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 5), las = 1)

# Posterior mean panel
brks <- pretty(s_mean, n = 11)
cols <- grDevices::hcl.colors(length(brks) - 1L, palette = "RdBu", rev = TRUE)
cuts <- cut(s_mean, breaks = brks, include.lowest = TRUE)
xy <- sim3$coords_raw
plot(xy[, 1] / 1000, xy[, 2] / 1000, pch = 19, cex = 0.55,
     col = cols[as.integer(cuts)],
     xlab = "Easting (km, OSGB)", ylab = "Northing (km, OSGB)",
     main = sprintf("Posterior mean s(loc)  (n=%d)", sim3$n),
     asp = 1)
# Manual horizontal legend strip on the right of the plot
usr <- par("usr")
lx <- usr[2] + 0.02 * diff(usr[1:2]); rx <- usr[2] + 0.05 * diff(usr[1:2])
ny <- length(brks) - 1L
ys <- seq(usr[3], usr[4], length.out = ny + 1)
for (k in seq_len(ny))
  rect(lx, ys[k], rx, ys[k + 1], col = cols[k], border = NA, xpd = NA)
text(rx + 0.01 * diff(usr[1:2]),
     ys[seq(1, ny + 1, length.out = 5)],
     labels = signif(brks[seq(1, ny + 1, length.out = 5)], 2),
     pos = 4, xpd = NA, cex = 0.7)

# Posterior SD panel
brks2 <- pretty(s_sd, n = 11)
cols2 <- grDevices::hcl.colors(length(brks2) - 1L, palette = "Inferno", rev = TRUE)
cuts2 <- cut(s_sd, breaks = brks2, include.lowest = TRUE)
plot(xy[, 1] / 1000, xy[, 2] / 1000, pch = 19, cex = 0.55,
     col = cols2[as.integer(cuts2)],
     xlab = "Easting (km, OSGB)", ylab = "Northing (km, OSGB)",
     main = sprintf("Posterior sd of s(loc)  (n=%d)", sim3$n),
     asp = 1)
usr <- par("usr")
lx <- usr[2] + 0.02 * diff(usr[1:2]); rx <- usr[2] + 0.05 * diff(usr[1:2])
ny <- length(brks2) - 1L
ys <- seq(usr[3], usr[4], length.out = ny + 1)
for (k in seq_len(ny))
  rect(lx, ys[k], rx, ys[k + 1], col = cols2[k], border = NA, xpd = NA)
text(rx + 0.01 * diff(usr[1:2]),
     ys[seq(1, ny + 1, length.out = 5)],
     labels = signif(brks2[seq(1, ny + 1, length.out = 5)], 2),
     pos = 4, xpd = NA, cex = 0.7)

par(op); dev.off()

# =============================================================
# 3. p3_interaction_surfaces.pdf  -- f_{u,v} with significance overlay
# =============================================================
cat("[plot 3/5] p3_interaction_surfaces.pdf\n")
u_grid <- sim3$x_grid_2d$u
v_grid <- sim3$x_grid_2d$v
stopifnot(length(u_grid) * length(v_grid) == nrow(fit3$f_int[[1]]))

pdf(file.path(OUT_DIR, "p3_interaction_surfaces.pdf"),
    width = 13.5, height = 4.5)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 5), las = 1)
for (key in names(fit3$f_int)) {
  ij  <- int_pair_idx[[key]]
  ju  <- ij[1]; jv <- ij[2]
  Mu  <- fit3$f_int[[key]]
  zmean <- rowMeans(Mu)
  zlo   <- apply(Mu, 1, quantile, 0.025)
  zhi   <- apply(Mu, 1, quantile, 0.975)
  # u varies fastest -> rows of matrix(., 30, 30) index u, cols index v
  Z   <- matrix(zmean, nrow = length(u_grid), ncol = length(v_grid))
  Zlo <- matrix(zlo,   nrow = length(u_grid), ncol = length(v_grid))
  Zhi <- matrix(zhi,   nrow = length(u_grid), ncol = length(v_grid))

  u_raw <- sim3$X_scale[[ju]]$lo + u_grid * (sim3$X_scale[[ju]]$hi - sim3$X_scale[[ju]]$lo)
  v_raw <- sim3$X_scale[[jv]]$lo + v_grid * (sim3$X_scale[[jv]]$hi - sim3$X_scale[[jv]]$lo)

  zrange <- range(Z)
  brks   <- pretty(zrange, n = 21)
  cols   <- grDevices::hcl.colors(length(brks) - 1L, palette = "BrBG", rev = FALSE)

  image(u_raw, v_raw, Z, col = cols, breaks = brks,
        xlab = sprintf("%s (raw)", covnames[ju]),
        ylab = sprintf("%s (raw)", covnames[jv]),
        main = sprintf("f_{%s} (%s x %s)", key, covnames[ju], covnames[jv]),
        useRaster = TRUE)
  # Significance mask: hatch cells whose 95% CI includes zero
  sig <- (Zlo > 0) | (Zhi < 0)
  if (any(!sig)) {
    # Overlay translucent grey on non-significant cells via a second image
    NS <- ifelse(sig, NA_real_, 1)
    image(u_raw, v_raw, NS, col = grDevices::adjustcolor("white", 0.55),
          add = TRUE, useRaster = TRUE)
  }
  # Contour lines at zero
  contour(u_raw, v_raw, Z, add = TRUE, levels = 0, lwd = 1.5, col = "black")

  # Color-strip legend
  usr <- par("usr")
  lx <- usr[2] + 0.04 * diff(usr[1:2]); rx <- usr[2] + 0.10 * diff(usr[1:2])
  ny <- length(brks) - 1L
  ys <- seq(usr[3], usr[4], length.out = ny + 1)
  for (k in seq_len(ny))
    rect(lx, ys[k], rx, ys[k + 1], col = cols[k], border = NA, xpd = NA)
  tk <- seq(1, ny + 1, length.out = 5)
  text(rx + 0.005 * diff(usr[1:2]),
       ys[tk], labels = signif(brks[tk], 2), pos = 4, xpd = NA, cex = 0.75)
}
mtext("P3 interaction surfaces: posterior mean, white overlay = 95% CI crosses 0, black contour = 0 level",
      side = 3, line = -1.4, outer = TRUE, cex = 0.9)
par(op); dev.off()

# =============================================================
# 4. p3_traces.pdf  -- variance components: trace + running mean + ESS
# =============================================================
cat("[plot 4/5] p3_traces.pdf\n")
vc <- fit3$var_comp
trace_panels <- list(
  list(x = vc$sigma2,            lab = "sigma^2"),
  list(x = vc$tau2_s,             lab = "tau^2 (spatial GP variance)"),
  list(x = vc$rho,                lab = "rho (Matern range)"),
  list(x = vc$tau2[[1]],          lab = "tau^2_{s,1}  (NO2)"),
  list(x = vc$tau2[[2]],          lab = "tau^2_{s,2}  (NDVI)"),
  list(x = vc$tau2[[3]],          lab = "tau^2_{s,3}  (IMD)"),
  list(x = vc$tau2_int[["1_2"]],  lab = "tau^2_{int,1_2}  (NO2 x NDVI)"),
  list(x = vc$tau2_int[["1_3"]],  lab = "tau^2_{int,1_3}  (NO2 x IMD)"),
  list(x = vc$tau2_int[["2_3"]],  lab = "tau^2_{int,2_3}  (NDVI x IMD)")
)

pdf(file.path(OUT_DIR, "p3_traces.pdf"), width = 12, height = 11)
op <- par(mfrow = c(3, 3), mar = c(3, 4, 3, 1), las = 1)
for (pn in trace_panels) {
  x <- pn$x
  ess <- ess_one(x)
  rm <- cumsum(x) / seq_along(x)
  plot(seq_along(x), x, type = "l", col = grDevices::adjustcolor("grey30", 0.5),
       xlab = "post-burn iteration",
       ylab = pn$lab,
       main = sprintf("%s   ESS=%.1f / %d",
                      pn$lab, ess, length(x)))
  lines(seq_along(x), rm, col = "firebrick3", lwd = 2)
  abline(h = mean(x), col = "steelblue4", lty = 2)
  # Highlight ESS < 100 in red title
  if (ess < 100) {
    mtext("LOW ESS", side = 3, line = -1.2, adj = 0.99, col = "firebrick3",
          cex = 0.8, font = 2)
  }
}
mtext("P3 (n=1500) variance-component traces. Red = running mean; blue dashed = posterior mean",
      side = 3, line = -1.4, outer = TRUE, cex = 0.95)
par(op); dev.off()

# =============================================================
# 5. p3_vs_p1plus  -- side-by-side credible bands and varcomp table
# =============================================================
cat("[plot 5/5] p3_vs_p1plus_smooths.pdf + p3_vs_p1plus_varcomp.txt\n")

# Re-center the n=500 fit the same way
smooths5 <- vector("list", 3)
for (j in 1:3) {
  pw <- post_pointwise(fit5$f_main[[j]])
  rc <- recenter_smooth(pw$mean, pw$lo, pw$hi, sim5$x_grid_1d,
                        sim5$data[[paste0("X", j)]])
  smooths5[[j]] <- rc
}

pdf(file.path(OUT_DIR, "p3_vs_p1plus_smooths.pdf"), width = 12, height = 4)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), las = 1)
for (j in 1:3) {
  # Use the n=1500 sim's X_scale as canonical for the raw axis
  Xlo <- sim3$X_scale[[j]]$lo; Xhi <- sim3$X_scale[[j]]$hi
  xs_raw <- Xlo + sim3$x_grid_1d * (Xhi - Xlo)
  xs5_raw <- sim5$X_scale[[j]]$lo +
             sim5$x_grid_1d * (sim5$X_scale[[j]]$hi - sim5$X_scale[[j]]$lo)
  rc3 <- smooths3[[j]]; rc5 <- smooths5[[j]]
  rng <- range(c(rc3$lo, rc3$hi, rc5$lo, rc5$hi))
  plot(xs_raw, rc3$mean, type = "n", ylim = rng,
       xlab = sprintf("%s (raw)", covnames[j]),
       ylab = sprintf("centred f(%s)", covnames[j]),
       main = covnames[j])
  polygon(c(xs5_raw, rev(xs5_raw)), c(rc5$lo, rev(rc5$hi)),
          col = grDevices::adjustcolor("darkorange", 0.18), border = NA)
  polygon(c(xs_raw,  rev(xs_raw)),  c(rc3$lo, rev(rc3$hi)),
          col = grDevices::adjustcolor("steelblue", 0.30), border = NA)
  lines(xs5_raw, rc5$mean, lwd = 2, col = "darkorange3", lty = 2)
  lines(xs_raw,  rc3$mean, lwd = 2, col = "steelblue4")
  abline(h = 0, lty = 3, col = "grey50")
  legend("topright",
         legend = c(sprintf("p3  (n=%d)", sim3$n),
                    sprintf("p1plus (n=%d)", sim5$n)),
         col = c("steelblue4", "darkorange3"),
         lwd = 2, lty = c(1, 2), bty = "n", cex = 0.8)
}
mtext("P3 vs P1plus: main-effect smooths, posterior mean + 95% CI",
      side = 3, line = -1.4, outer = TRUE, cex = 1.0)
par(op); dev.off()

# Variance-component comparison table
summ_one <- function(x) c(mean = mean(x),
                          q025 = unname(quantile(x, 0.025)),
                          q975 = unname(quantile(x, 0.975)),
                          ess  = ess_one(x))
rows <- list(
  "sigma2 (noise)"       = list(p5 = vc5_summ <- summ_one(fit5$var_comp$sigma2),
                                 p3 = summ_one(fit3$var_comp$sigma2)),
  "tau2_s (spatial)"     = list(p5 = summ_one(fit5$var_comp$tau2_s),
                                 p3 = summ_one(fit3$var_comp$tau2_s)),
  "rho (Matern range)"   = list(p5 = summ_one(fit5$var_comp$rho),
                                 p3 = summ_one(fit3$var_comp$rho)),
  "tau2_s_main[NO2]"     = list(p5 = summ_one(fit5$var_comp$tau2[[1]]),
                                 p3 = summ_one(fit3$var_comp$tau2[[1]])),
  "tau2_s_main[NDVI]"    = list(p5 = summ_one(fit5$var_comp$tau2[[2]]),
                                 p3 = summ_one(fit3$var_comp$tau2[[2]])),
  "tau2_s_main[IMD]"     = list(p5 = summ_one(fit5$var_comp$tau2[[3]]),
                                 p3 = summ_one(fit3$var_comp$tau2[[3]])),
  "tau2_int[NO2xNDVI]"   = list(p5 = summ_one(fit5$var_comp$tau2_int[["1_2"]]),
                                 p3 = summ_one(fit3$var_comp$tau2_int[["1_2"]])),
  "tau2_int[NO2xIMD]"    = list(p5 = summ_one(fit5$var_comp$tau2_int[["1_3"]]),
                                 p3 = summ_one(fit3$var_comp$tau2_int[["1_3"]])),
  "tau2_int[NDVIxIMD]"   = list(p5 = summ_one(fit5$var_comp$tau2_int[["2_3"]]),
                                 p3 = summ_one(fit3$var_comp$tau2_int[["2_3"]]))
)

txt_path <- file.path(OUT_DIR, "p3_vs_p1plus_varcomp.txt")
sink(txt_path)
cat("==============================================================\n")
cat("  P3 (n=1500) vs P1plus (n=500)  -- variance components\n")
cat("  fit3 RDS: ", FIT_P3_RDS,  "\n")
cat("  fit5 RDS: ", FIT_P1P_RDS, "\n")
cat("==============================================================\n\n")
header <- sprintf("%-22s | %10s %10s %10s %6s | %10s %10s %10s %6s",
                  "parameter",
                  "p5_mean","p5_q025","p5_q975","p5_ESS",
                  "p3_mean","p3_q025","p3_q975","p3_ESS")
cat(header, "\n", sep = "")
cat(strrep("-", nchar(header)), "\n", sep = "")
for (nm in names(rows)) {
  r <- rows[[nm]]
  cat(sprintf("%-22s | %10.4g %10.4g %10.4g %6.0f | %10.4g %10.4g %10.4g %6.0f\n",
              nm,
              r$p5["mean"], r$p5["q025"], r$p5["q975"], r$p5["ess"],
              r$p3["mean"], r$p3["q025"], r$p3["q975"], r$p3["ess"]))
}
cat("\n")
cat(sprintf("MH accept p3: sigma2=%.3f tau2=%.3f rho=%.3f  (wall=%.2f h)\n",
            fit3$convergence$mh_accept[["sigma2"]],
            fit3$convergence$mh_accept[["tau2"]],
            fit3$convergence$mh_accept[["rho"]],
            fit3$reproducibility$wall_sec / 3600))
cat(sprintf("MH accept p5: sigma2=%.3f tau2=%.3f rho=%.3f\n",
            fit5$convergence$mh_accept[["sigma2"]],
            fit5$convergence$mh_accept[["tau2"]],
            fit5$convergence$mh_accept[["rho"]]))
cat("\nNotes:\n")
cat(" - p5 = fit_medsat_p1plus_allpairs_n500_seed1_local.rds (n=500, all-pairs)\n")
cat(" - p3 = fit_medsat_n1500_seed1_allpairs_hb.rds          (n=1500, all-pairs)\n")
cat(" - ESS reported out of 2000 post-burn draws in both runs.\n")
cat(" - tau2_int ESS is the headline mixing concern at n=1500\n")
cat("   (8-21 / 2000); we should not lean on these blocks for inference\n")
cat("   from this single run alone.\n")
sink()

# Echo to stdout as well
cat("\n----- p3_vs_p1plus_varcomp.txt -----\n")
cat(readLines(txt_path), sep = "\n")
cat("\n----- end -----\n")

# =============================================================
# Final report
# =============================================================
cat("\n==============================================================\n")
cat(" P3 diagnostics summary\n")
cat("==============================================================\n")
cat(sprintf("Outputs written under %s/:\n", OUT_DIR))
for (f in c("p3_smooths.pdf","p3_spatial_re.pdf","p3_interaction_surfaces.pdf",
            "p3_traces.pdf","p3_vs_p1plus_smooths.pdf","p3_vs_p1plus_varcomp.txt")) {
  fp <- file.path(OUT_DIR, f)
  if (file.exists(fp)) {
    sz <- file.info(fp)$size
    cat(sprintf("  %-36s  (%6.1f KB)\n", f, sz/1024))
  } else {
    cat(sprintf("  %-36s  [MISSING]\n", f))
  }
}

# Print interaction-significance summary inline
cat("\nInteraction-surface 95% CI excludes zero in:\n")
for (key in names(fit3$f_int)) {
  ij <- int_pair_idx[[key]]
  Mu <- fit3$f_int[[key]]
  zlo <- apply(Mu, 1, quantile, 0.025)
  zhi <- apply(Mu, 1, quantile, 0.975)
  sig <- (zlo > 0) | (zhi < 0)
  cat(sprintf("  %s (%s x %s):  %d / %d cells (%.1f%%) significant\n",
              key, covnames[ij[1]], covnames[ij[2]],
              sum(sig), length(sig), 100 * mean(sig)))
}

cat("\nDone.\n")
