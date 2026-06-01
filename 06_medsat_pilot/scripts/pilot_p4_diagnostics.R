## ============================================================
## scripts/pilot_p4_diagnostics.R
##
## Comparative diagnostics across the four n=1500 fits:
##   hl1 -- seed 1, longer chain (n_iter=6000, n_thin=2)   [HEADLINE]
##   s2  -- seed 2, baseline chain (n_iter=3000)
##   s3  -- seed 3, baseline chain
##   s4  -- seed 4, baseline chain
##
## Produces under output/plots/p4/:
##   p4_smooths_all_seeds.pdf
##   p4_interaction_surfaces_all_seeds.pdf
##   p4_spatial_re_all_seeds.pdf
##   p4_variance_table.csv
##   p4_credible_zero_summary.txt
##   p4_trace_grid.pdf
##   p4_ess_comparison.csv
##
## Read-only on the RDS files; does not refit.
## ============================================================

suppressPackageStartupMessages({ library(coda) })

OUT_DIR <- "output/plots/p4"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SEED_TAGS <- c("hl1", "s2", "s3", "s4")
SEED_LABS <- c(hl1 = "seed 1 (headline)", s2 = "seed 2", s3 = "seed 3", s4 = "seed 4")
SEED_COL  <- c(hl1 = "steelblue4", s2 = "darkorange3",
               s3 = "forestgreen", s4 = "purple3")
SEED_LTY  <- c(hl1 = 1, s2 = 2, s3 = 3, s4 = 4)

FIT_PATH <- c(
  hl1 = "output/fit_medsat_n1500_seed1_headline_hb.rds",
  s2  = "output/fit_medsat_n1500_seed2_allpairs_hb.rds",
  s3  = "output/fit_medsat_n1500_seed3_allpairs_hb.rds",
  s4  = "output/fit_medsat_n1500_seed4_allpairs_hb.rds"
)
SIM_PATH <- c(
  hl1 = "output/sim_medsat_london_asthma_n1500_seed1_v2.rds",
  s2  = "output/sim_medsat_london_asthma_n1500_seed2_v2.rds",
  s3  = "output/sim_medsat_london_asthma_n1500_seed3_v2.rds",
  s4  = "output/sim_medsat_london_asthma_n1500_seed4_v2.rds"
)
stopifnot(all(file.exists(FIT_PATH)), all(file.exists(SIM_PATH)))

cat("[load] reading 4 fits + 4 sims...\n")
fits <- lapply(FIT_PATH, readRDS)
sims <- lapply(SIM_PATH, readRDS)
names(fits) <- SEED_TAGS; names(sims) <- SEED_TAGS

# Structural sanity
for (tg in SEED_TAGS) {
  stopifnot(isTRUE(fits[[tg]]$settings$fit_interactions))
  stopifnot(length(fits[[tg]]$f_int) == 3L)
  stopifnot(setequal(names(fits[[tg]]$f_int), c("1_2","1_3","2_3")))
  stopifnot(sims[[tg]]$n == 1500L, sims[[tg]]$p == 3L)
}
n_draws <- fits$hl1$settings$n_draws
stopifnot(all(vapply(fits, function(f) f$settings$n_draws, integer(1)) == n_draws))

covnames <- c("NO2", "NDVI", "IMD")
int_pair_idx <- list("1_2" = c(1L, 2L), "1_3" = c(1L, 3L), "2_3" = c(2L, 3L))

# ============================================================
# Shared helpers
# ============================================================
ess_one <- function(x) as.numeric(coda::effectiveSize(coda::as.mcmc(x)))

post_pointwise <- function(M, probs = c(0.025, 0.975)) {
  list(mean = rowMeans(M),
       lo   = apply(M, 1, quantile, probs[1]),
       hi   = apply(M, 1, quantile, probs[2]))
}

# Re-center a smooth's posterior mean curve against an EXTERNAL reference
# set of observed X_j values (here: seed 1's data). This puts all seeds on
# the same vertical anchor for shape comparison.
recenter_against <- function(curve_mean_or_band, xs, ref_obs_X) {
  curve_at_ref <- stats::approx(x = xs, y = curve_mean_or_band,
                                xout = ref_obs_X, rule = 2)$y
  curve_mean_or_band - mean(curve_at_ref)
}

# ============================================================
# 1. p4_smooths_all_seeds.pdf — 3 panels, 4 curves each, on raw X
# ============================================================
cat("\n[1/7] p4_smooths_all_seeds.pdf\n")
ref_sim <- sims$hl1
xs <- ref_sim$x_grid_1d

# Per-(seed, j) re-centered mean and bands
smooths <- list()
for (tg in SEED_TAGS) {
  smooths[[tg]] <- lapply(seq_len(3), function(j) {
    pw <- post_pointwise(fits[[tg]]$f_main[[j]])
    # Re-center each seed's posterior mean curve against seed 1's data
    anc <- mean(stats::approx(x = xs, y = pw$mean,
                              xout = ref_sim$data[[paste0("X", j)]],
                              rule = 2)$y)
    list(mean = pw$mean - anc,
         lo   = pw$lo   - anc,
         hi   = pw$hi   - anc)
  })
}

pdf(file.path(OUT_DIR, "p4_smooths_all_seeds.pdf"), width = 13, height = 4.4)
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.4, 3.3, 1), las = 1, cex.main = 0.95)
for (j in 1:3) {
  Xlo <- ref_sim$X_scale[[j]]$lo; Xhi <- ref_sim$X_scale[[j]]$hi
  xs_raw <- Xlo + xs * (Xhi - Xlo)
  curves <- sapply(SEED_TAGS, function(tg) smooths[[tg]][[j]]$mean)  # 101 x 4
  rng <- range(unlist(lapply(SEED_TAGS, function(tg)
    c(smooths[[tg]][[j]]$lo, smooths[[tg]][[j]]$hi))))

  # Average pairwise cor of posterior-mean curves
  cor_mat <- cor(curves)
  pair_cors <- cor_mat[upper.tri(cor_mat)]
  cor_summary <- sprintf("mean pairwise cor = %.3f  (range [%.3f, %.3f])",
                         mean(pair_cors), min(pair_cors), max(pair_cors))

  plot(xs_raw, curves[, 1], type = "n", ylim = rng,
       xlab = sprintf("%s (raw units)", covnames[j]),
       ylab = sprintf("re-centred f(%s)", covnames[j]),
       main = sprintf("%s\n%s", covnames[j], cor_summary))
  # Light shading of headline credible band
  polygon(c(xs_raw, rev(xs_raw)),
          c(smooths$hl1[[j]]$lo, rev(smooths$hl1[[j]]$hi)),
          col = grDevices::adjustcolor(SEED_COL["hl1"], 0.18), border = NA)
  abline(h = 0, lty = 3, col = "grey50")
  for (tg in SEED_TAGS) {
    lines(xs_raw, smooths[[tg]][[j]]$mean,
          lwd = 2, col = SEED_COL[tg], lty = SEED_LTY[tg])
  }
  rug(Xlo + ref_sim$data[[paste0("X", j)]] * (Xhi - Xlo),
      ticksize = 0.018, col = grDevices::adjustcolor("grey40", 0.4))
  legend("topright",
         legend = SEED_LABS[SEED_TAGS],
         col    = SEED_COL[SEED_TAGS],
         lty    = SEED_LTY[SEED_TAGS],
         lwd = 2, bty = "n", cex = 0.75)
}
mtext("P4: posterior-mean main-effect smooths across 4 random subsamples (n=1500). Shaded band = headline 95% CI.",
      side = 3, line = -1.4, outer = TRUE, cex = 0.9)
par(op); dev.off()

# ============================================================
# 2. p4_interaction_surfaces_all_seeds.pdf — 3x4 grid
# ============================================================
cat("[2/7] p4_interaction_surfaces_all_seeds.pdf\n")
u_grid <- fits$hl1$x_grid_2d$u
v_grid <- fits$hl1$x_grid_2d$v
n_int_cells <- length(u_grid) * length(v_grid)

# Precompute pointwise summaries for each (seed, key)
int_summ <- list()
sig_count <- matrix(0L, nrow = 4, ncol = 3,
                    dimnames = list(SEED_TAGS, c("1_2","1_3","2_3")))
max_abs_mean <- matrix(0.0, nrow = 4, ncol = 3,
                       dimnames = list(SEED_TAGS, c("1_2","1_3","2_3")))
for (tg in SEED_TAGS) {
  int_summ[[tg]] <- list()
  for (key in c("1_2","1_3","2_3")) {
    Mu <- fits[[tg]]$f_int[[key]]
    z_mean <- rowMeans(Mu)
    z_lo   <- apply(Mu, 1, quantile, 0.025)
    z_hi   <- apply(Mu, 1, quantile, 0.975)
    sig    <- (z_lo > 0) | (z_hi < 0)
    int_summ[[tg]][[key]] <- list(mean = z_mean, lo = z_lo, hi = z_hi, sig = sig)
    sig_count[tg, key]    <- sum(sig)
    max_abs_mean[tg, key] <- max(abs(z_mean))
  }
}

# Shared z range across all 12 surfaces for honest visual comparison
z_all <- unlist(lapply(int_summ, function(L) lapply(L, function(x) x$mean)))
z_lim <- max(abs(z_all)) * 1.02
brks <- pretty(c(-z_lim, z_lim), n = 21)
cols <- grDevices::hcl.colors(length(brks) - 1L, palette = "BrBG", rev = FALSE)

pdf(file.path(OUT_DIR, "p4_interaction_surfaces_all_seeds.pdf"),
    width = 13.5, height = 10.5)
op <- par(mfcol = c(3, 4), mar = c(3.7, 3.7, 2.6, 1.0), las = 1,
          oma = c(0, 0, 3, 4))
for (tg in SEED_TAGS) {
  for (key in c("1_2", "1_3", "2_3")) {
    ij <- int_pair_idx[[key]]
    ju <- ij[1]; jv <- ij[2]
    L  <- int_summ[[tg]][[key]]
    Z  <- matrix(L$mean, nrow = length(u_grid), ncol = length(v_grid))
    # Raw-axis u/v using THIS seed's X_scale for self-consistency
    u_raw <- sims[[tg]]$X_scale[[ju]]$lo + u_grid *
             (sims[[tg]]$X_scale[[ju]]$hi - sims[[tg]]$X_scale[[ju]]$lo)
    v_raw <- sims[[tg]]$X_scale[[jv]]$lo + v_grid *
             (sims[[tg]]$X_scale[[jv]]$hi - sims[[tg]]$X_scale[[jv]]$lo)
    image(u_raw, v_raw, Z, col = cols, breaks = brks, useRaster = TRUE,
          xlab = sprintf("%s", covnames[ju]),
          ylab = sprintf("%s", covnames[jv]),
          main = sprintf("%s  %sx%s   sig=%d/%d  max|mean|=%.3g",
                         SEED_LABS[tg], covnames[ju], covnames[jv],
                         sig_count[tg, key], n_int_cells,
                         max_abs_mean[tg, key]))
    contour(u_raw, v_raw, Z, add = TRUE, levels = 0, lwd = 1.4, col = "black")
  }
}
# Shared color legend on the right margin
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
    new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
ny <- length(brks) - 1L
ys <- seq(0.08, 0.92, length.out = ny + 1)
for (k in seq_len(ny))
  rect(0.965, ys[k], 0.985, ys[k + 1], col = cols[k], border = NA, xpd = NA)
tk <- seq(1, ny + 1, length.out = 5)
text(0.988, ys[tk], labels = signif(brks[tk], 2), pos = 4, xpd = NA, cex = 0.7)
text(0.975, 0.96, labels = "f_{u,v}", pos = NULL, xpd = NA, cex = 0.85)
mtext("P4 interaction surfaces (rows = pair, columns = seed). Shared colour scale; black contour = 0 level.",
      side = 3, line = -1, outer = TRUE, cex = 0.95)
par(op); dev.off()

# ============================================================
# 3. p4_spatial_re_all_seeds.pdf — 4 panels with shared color scale
# ============================================================
cat("[3/7] p4_spatial_re_all_seeds.pdf\n")
s_mean_list <- lapply(SEED_TAGS, function(tg) rowMeans(fits[[tg]]$s_obs))
names(s_mean_list) <- SEED_TAGS
abs_lim <- max(abs(unlist(s_mean_list)))
brks_s <- seq(-abs_lim, abs_lim, length.out = 21)
cols_s <- grDevices::hcl.colors(length(brks_s) - 1L, palette = "RdBu", rev = TRUE)

pdf(file.path(OUT_DIR, "p4_spatial_re_all_seeds.pdf"), width = 13, height = 11)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 5), las = 1)
for (tg in SEED_TAGS) {
  s_mean <- s_mean_list[[tg]]
  cuts   <- cut(s_mean, breaks = brks_s, include.lowest = TRUE)
  xy <- sims[[tg]]$coords_raw
  plot(xy[, 1] / 1000, xy[, 2] / 1000, pch = 19, cex = 0.6,
       col = cols_s[as.integer(cuts)],
       xlab = "Easting (km, OSGB)", ylab = "Northing (km, OSGB)",
       main = sprintf("%s : posterior mean s(loc)", SEED_LABS[tg]),
       asp = 1)
  usr <- par("usr")
  lx <- usr[2] + 0.02 * diff(usr[1:2]); rx <- usr[2] + 0.05 * diff(usr[1:2])
  ny <- length(brks_s) - 1L
  ys <- seq(usr[3], usr[4], length.out = ny + 1)
  for (k in seq_len(ny))
    rect(lx, ys[k], rx, ys[k + 1], col = cols_s[k], border = NA, xpd = NA)
  tk <- seq(1, ny + 1, length.out = 5)
  text(rx + 0.01 * diff(usr[1:2]),
       ys[tk], labels = signif(brks_s[tk], 2),
       pos = 4, xpd = NA, cex = 0.7)
}
mtext("P4 spatial random effect: posterior mean across 4 random subsamples. Shared colour scale.",
      side = 3, line = -1.4, outer = TRUE, cex = 0.95)
par(op); dev.off()

# ============================================================
# 4. p4_variance_table.csv
# ============================================================
cat("[4/7] p4_variance_table.csv\n")
# Define parameter list with extractor closures
param_defs <- list(
  list(label = "sigma2 (noise)",            getter = function(f) f$var_comp$sigma2,
       across = TRUE),
  list(label = "tau2 (spatial GP variance)",getter = function(f) f$var_comp$tau2_s,
       across = TRUE),
  list(label = "rho (Matern range)",        getter = function(f) f$var_comp$rho,
       across = TRUE),
  list(label = "tau2_s_main[NO2]",          getter = function(f) f$var_comp$tau2[[1]],
       across = FALSE),
  list(label = "tau2_s_main[NDVI]",         getter = function(f) f$var_comp$tau2[[2]],
       across = FALSE),
  list(label = "tau2_s_main[IMD]",          getter = function(f) f$var_comp$tau2[[3]],
       across = FALSE),
  list(label = "tau2_int[NO2xNDVI]",        getter = function(f) f$var_comp$tau2_int[["1_2"]],
       across = FALSE),
  list(label = "tau2_int[NO2xIMD]",         getter = function(f) f$var_comp$tau2_int[["1_3"]],
       across = FALSE),
  list(label = "tau2_int[NDVIxIMD]",        getter = function(f) f$var_comp$tau2_int[["2_3"]],
       across = FALSE)
)

fmt_cell <- function(x) {
  q <- quantile(x, c(0.025, 0.975))
  e <- ess_one(x)
  sprintf("%.4g (%.4g, %.4g) [%.0f]", mean(x), q[1], q[2], e)
}

rows_list <- list()
for (pd in param_defs) {
  cells <- vapply(SEED_TAGS, function(tg) fmt_cell(pd$getter(fits[[tg]])),
                  character(1))
  if (pd$across) {
    means <- vapply(SEED_TAGS, function(tg) mean(pd$getter(fits[[tg]])),
                    numeric(1))
    across_sd <- sprintf("%.4g", sd(means))
  } else {
    across_sd <- ""
  }
  rows_list[[length(rows_list) + 1L]] <- c(parameter = pd$label,
                                            hl1 = unname(cells["hl1"]),
                                            s2  = unname(cells["s2"]),
                                            s3  = unname(cells["s3"]),
                                            s4  = unname(cells["s4"]),
                                            across_seed_SD_of_means = across_sd)
}
varcomp_df <- as.data.frame(do.call(rbind, rows_list), stringsAsFactors = FALSE)
write.csv(varcomp_df, file.path(OUT_DIR, "p4_variance_table.csv"),
          row.names = FALSE)

# ============================================================
# 5. p4_credible_zero_summary.txt
# ============================================================
cat("[5/7] p4_credible_zero_summary.txt\n")
txt_path <- file.path(OUT_DIR, "p4_credible_zero_summary.txt")
sink(txt_path)
cat("==============================================================\n")
cat(" P4 credible-band-excludes-zero summary\n")
cat(" 30 x 30 = 900 grid points per interaction surface\n")
cat(" 4 random subsamples (n=1500) x 3 interaction surfaces = 12 cells\n")
cat("==============================================================\n\n")
header <- sprintf("%-22s | %-6s %-6s %-6s %-6s | %s",
                  "interaction", "hl1", "s2", "s3", "s4", "max |mean| across seeds")
cat(header, "\n", sep = "")
cat(strrep("-", nchar(header)), "\n", sep = "")
for (key in c("1_2","1_3","2_3")) {
  ij <- int_pair_idx[[key]]
  label <- sprintf("%s x %s (%s)", covnames[ij[1]], covnames[ij[2]], key)
  cat(sprintf("%-22s | %3d/%-3d %3d/%-3d %3d/%-3d %3d/%-3d | %.4g\n",
              label,
              sig_count["hl1", key], n_int_cells,
              sig_count["s2",  key], n_int_cells,
              sig_count["s3",  key], n_int_cells,
              sig_count["s4",  key], n_int_cells,
              max(max_abs_mean[, key])))
}
cat("\n")
cat(sprintf("TOTAL grid points with 95%% CI excluding zero across all 12 surfaces: %d / %d\n",
            sum(sig_count), 4 * 3 * n_int_cells))
cat("Per surface (mean over seeds):\n")
for (key in c("1_2","1_3","2_3")) {
  ij <- int_pair_idx[[key]]
  cat(sprintf("  %s x %s : mean %.1f / %d cells across seeds\n",
              covnames[ij[1]], covnames[ij[2]],
              mean(sig_count[, key]), n_int_cells))
}
cat("\n")
cat("Max |posterior-mean| on the grid (per seed x pair):\n")
print(round(max_abs_mean, 4))
cat("\nFor reference, sd(y) at each seed:\n")
for (tg in SEED_TAGS)
  cat(sprintf("  %-18s : %.4f\n", SEED_LABS[tg], sd(sims[[tg]]$data$y)))
sink()

# ============================================================
# 6. p4_trace_grid.pdf — 12-panel: {sigma2, tau2_s, rho} x 4 seeds
# ============================================================
cat("[6/7] p4_trace_grid.pdf\n")
trace_specs <- list(
  list(name = "sigma2",  getter = function(f) f$var_comp$sigma2),
  list(name = "tau2_s",  getter = function(f) f$var_comp$tau2_s),
  list(name = "rho",     getter = function(f) f$var_comp$rho)
)

pdf(file.path(OUT_DIR, "p4_trace_grid.pdf"), width = 13, height = 9)
op <- par(mfrow = c(3, 4), mar = c(3, 4, 3, 1), las = 1, cex.main = 0.9)
for (ts in trace_specs) {
  for (tg in SEED_TAGS) {
    x <- ts$getter(fits[[tg]])
    e <- ess_one(x)
    rm <- cumsum(x) / seq_along(x)
    plot(seq_along(x), x, type = "l",
         col = grDevices::adjustcolor("grey30", 0.5),
         xlab = "post-burn iteration",
         ylab = ts$name,
         main = sprintf("%s | %s   ESS=%.0f", SEED_LABS[tg], ts$name, e))
    lines(seq_along(x), rm, col = "firebrick3", lwd = 2)
    abline(h = mean(x), col = "steelblue4", lty = 2)
  }
}
mtext("P4 traces: rows {sigma2, tau2_s (spatial), rho}, columns = seeds. Red = running mean; blue dashed = posterior mean.",
      side = 3, line = -1.2, outer = TRUE, cex = 0.95)
par(op); dev.off()

# ============================================================
# 7. p4_ess_comparison.csv
# ============================================================
cat("[7/7] p4_ess_comparison.csv\n")
ess_rows <- list()
for (pd in param_defs) {
  ess_vals <- vapply(SEED_TAGS, function(tg) ess_one(pd$getter(fits[[tg]])),
                     numeric(1))
  baseline_mean <- mean(ess_vals[c("s2","s3","s4")])
  ratio <- if (baseline_mean > 0) ess_vals["hl1"] / baseline_mean else NA_real_
  ess_rows[[length(ess_rows) + 1L]] <- c(
    parameter = pd$label,
    hl1 = sprintf("%.1f", unname(ess_vals["hl1"])),
    s2  = sprintf("%.1f", unname(ess_vals["s2"])),
    s3  = sprintf("%.1f", unname(ess_vals["s3"])),
    s4  = sprintf("%.1f", unname(ess_vals["s4"])),
    baseline_mean_s2_s4 = sprintf("%.1f", baseline_mean),
    headline_over_baseline = sprintf("%.2f", unname(ratio))
  )
}
ess_df <- as.data.frame(do.call(rbind, ess_rows), stringsAsFactors = FALSE)
write.csv(ess_df, file.path(OUT_DIR, "p4_ess_comparison.csv"),
          row.names = FALSE)

# ============================================================
# Console summary
# ============================================================
cat("\n==============================================================\n")
cat(" P4 diagnostics summary\n")
cat("==============================================================\n")
cat(sprintf("Outputs written to %s/:\n", OUT_DIR))
for (f in c("p4_smooths_all_seeds.pdf",
            "p4_interaction_surfaces_all_seeds.pdf",
            "p4_spatial_re_all_seeds.pdf",
            "p4_variance_table.csv",
            "p4_credible_zero_summary.txt",
            "p4_trace_grid.pdf",
            "p4_ess_comparison.csv")) {
  fp <- file.path(OUT_DIR, f)
  if (file.exists(fp)) {
    cat(sprintf("  %-44s  (%7.1f KB)\n", f, file.info(fp)$size / 1024))
  } else {
    cat(sprintf("  %-44s  [MISSING]\n", f))
  }
}

cat("\n=== Credible-band exclusion of zero (sig / 900 cells) ===\n")
print(sig_count)
cat(sprintf("\nTotal across 12 surfaces: %d / %d (%.3f%%)\n",
            sum(sig_count), 4 * 3 * n_int_cells,
            100 * sum(sig_count) / (4 * 3 * n_int_cells)))

cat("\n=== Across-seed SD of posterior MEANS (vs typical within-seed posterior SD) ===\n")
for (pd in param_defs[1:3]) {
  means   <- vapply(SEED_TAGS, function(tg) mean(pd$getter(fits[[tg]])),  numeric(1))
  sds_in  <- vapply(SEED_TAGS, function(tg) sd(pd$getter(fits[[tg]])),    numeric(1))
  cat(sprintf("  %-30s  across-seed SD = %.4g   |  within-seed SD range = [%.4g, %.4g]\n",
              pd$label, sd(means), min(sds_in), max(sds_in)))
}

cat("\nDone.\n")
