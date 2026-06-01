## ============================================================
## scripts/pilot_p0_mgcv.R
##
## P0 sanity-check fit using mgcv's penalized GAM. Two models:
##   - m_additive : s(NO2) + s(NDVI) + s(IMD) + s(lon,lat)
##   - m_interact : m_additive + ti(NO2, NDVI)
##
## Purpose: confirm a standard penalized-spline fit produces sensible
## smooths on this dataset BEFORE we commit to an expensive Bayesian
## fit. NOT a scientific result.
## ============================================================

suppressPackageStartupMessages({
  library(mgcv)
})

SIM_RDS  <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"
OUT_FIT  <- "output/fit_medsat_p0_mgcv_local.rds"
OUT_DIR  <- "output/plots/p0_mgcv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(SIM_RDS))
sim <- readRDS(SIM_RDS)

# ----------------------------------------------------------------
# Fit
# ----------------------------------------------------------------
cat("====== Fitting m_additive ======\n")
t_add <- system.time(
  m_additive <- gam(
    y ~ s(X1, bs = "cr", k = 15) +
        s(X2, bs = "cr", k = 15) +
        s(X3, bs = "cr", k = 15) +
        s(lon, lat, bs = "tp", k = 50),
    data   = sim$data,
    method = "REML"
  )
)
cat("\n====== Fitting m_interact ======\n")
t_int <- system.time(
  m_interact <- gam(
    y ~ s(X1, bs = "cr", k = 15) +
        s(X2, bs = "cr", k = 15) +
        s(X3, bs = "cr", k = 15) +
        ti(X1, X2, bs = "cr", k = 8) +
        s(lon, lat, bs = "tp", k = 50),
    data   = sim$data,
    method = "REML"
  )
)
cat(sprintf("\nadditive runtime: %.2fs   interact runtime: %.2fs\n",
            t_add[["elapsed"]], t_int[["elapsed"]]))

# ----------------------------------------------------------------
# Save fits
# ----------------------------------------------------------------
saveRDS(list(
  additive = m_additive,
  interact = m_interact,
  sim_file = basename(SIM_RDS),
  fit_date = Sys.time(),
  hostname = Sys.info()[["nodename"]],
  runtime  = list(additive = unname(t_add[["elapsed"]]),
                  interact = unname(t_int[["elapsed"]]))
), OUT_FIT)
cat(sprintf("Saved fits: %s   (%.1f KB)\n",
            OUT_FIT, file.info(OUT_FIT)$size / 1024))

# ----------------------------------------------------------------
# Helpers for original-scale-X smooth plots
# ----------------------------------------------------------------
# Plot a 1-D smooth term against its raw-scale axis, with 95% CI band
# and partial residuals. Works for terms named s(Xj).
plot_smooth_rawx <- function(model, j, sim, ...) {
  Xlo  <- sim$X_scale[[j]]$lo
  Xhi  <- sim$X_scale[[j]]$hi
  nm   <- sim$X_scale[[j]]$name
  term <- sprintf("s(X%d)", j)

  xs_scaled <- seq(0, 1, length.out = 200)
  newdf <- as.data.frame(matrix(
    rep(c(mean(sim$data$X1), mean(sim$data$X2), mean(sim$data$X3),
          mean(sim$data$lon), mean(sim$data$lat)),
        each = length(xs_scaled)),
    ncol = 5))
  names(newdf) <- c("X1", "X2", "X3", "lon", "lat")
  newdf[[paste0("X", j)]] <- xs_scaled

  pred <- predict(model, newdata = newdf, type = "terms", se.fit = TRUE)
  ix <- match(term, colnames(pred$fit))
  if (is.na(ix)) stop("term ", term, " not found in predict.gam")
  fit <- pred$fit[, ix]
  se  <- pred$se.fit[, ix]

  pred_obs <- predict(model, type = "terms")
  pres <- pred_obs[, ix] + residuals(model)

  xs_raw  <- Xlo + xs_scaled * (Xhi - Xlo)
  obs_raw <- Xlo + sim$data[[paste0("X", j)]] * (Xhi - Xlo)

  edf <- summary(model)$s.table[term, "edf"]
  pval <- summary(model)$s.table[term, "p-value"]

  ylim <- range(c(fit + 2 * se, fit - 2 * se, pres))
  plot(xs_raw, fit, type = "n", ylim = ylim,
       xlab = sprintf("%s (raw)", nm),
       ylab = sprintf("s(%s) - centred", nm),
       main = sprintf("s(%s)   edf=%.1f   p=%.2g", nm, edf, pval),
       ...)
  polygon(c(xs_raw, rev(xs_raw)),
          c(fit + 2 * se, rev(fit - 2 * se)),
          col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
  points(obs_raw, pres, pch = 20, col = grDevices::adjustcolor("black", 0.30),
         cex = 0.5)
  lines(xs_raw, fit, lwd = 2, col = "steelblue4")
  abline(h = 0, lty = 3, col = "grey50")
}

# ----------------------------------------------------------------
# (a) p0_smooths_additive.pdf
# ----------------------------------------------------------------
pdf(file.path(OUT_DIR, "p0_smooths_additive.pdf"),
    width = 10, height = 8)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in 1:3) plot_smooth_rawx(m_additive, j, sim)
plot(m_additive, select = 4, scheme = 2, too.far = 0.1,
     main = "s(lon, lat)", asp = 1)
par(op); dev.off()

# ----------------------------------------------------------------
# (b) p0_smooths_interact.pdf  (5 smooths)
# ----------------------------------------------------------------
pdf(file.path(OUT_DIR, "p0_smooths_interact.pdf"),
    width = 12, height = 8)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (j in 1:3) plot_smooth_rawx(m_interact, j, sim)
# ti(X1, X2) — show on original-scale axes too
plot(m_interact, select = 4, scheme = 2, too.far = 0.1,
     main = "ti(NO2, NDVI)",
     xlab = "X1 (scaled NO2)", ylab = "X2 (scaled NDVI)")
plot(m_interact, select = 5, scheme = 2, too.far = 0.1,
     main = "s(lon, lat)", asp = 1)
plot.new()
par(op); dev.off()

# ----------------------------------------------------------------
# (c) p0_residuals.pdf  — gam.check 4-panel + text capture
# ----------------------------------------------------------------
pdf(file.path(OUT_DIR, "p0_residuals.pdf"),
    width = 9, height = 9)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
gam_check_text <- capture.output(gam.check(m_additive))
par(op); dev.off()

# ----------------------------------------------------------------
# (d) p0_interaction_surface.pdf  — vis.gam contour + perspective
# ----------------------------------------------------------------
pdf(file.path(OUT_DIR, "p0_interaction_surface.pdf"),
    width = 10, height = 5)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
vis.gam(m_interact, view = c("X1", "X2"),
        plot.type = "contour", color = "topo", too.far = 0.10,
        main = "Fitted surface vs (NO2, NDVI)\n(other covs at median)",
        xlab = "X1 (scaled NO2)", ylab = "X2 (scaled NDVI)")
vis.gam(m_interact, view = c("X1", "X2"),
        plot.type = "persp", color = "topo", too.far = 0.10,
        theta = 30, phi = 30,
        main = "Fitted surface (perspective)",
        xlab = "X1", ylab = "X2", zlab = "y")
par(op); dev.off()

# ----------------------------------------------------------------
# Report
# ----------------------------------------------------------------
cat("\n==============================================================\n")
cat(" PER-MODEL SUMMARY\n")
cat("==============================================================\n")

report_model <- function(model, name, et) {
  s <- summary(model)
  cat(sprintf("\n--- %s ---\n", name))
  cat(sprintf("  runtime         : %.2fs\n", et))
  cat(sprintf("  R-sq (adj)      : %.4f\n", s$r.sq))
  cat(sprintf("  Deviance expl.  : %.2f%%\n", 100 * s$dev.expl))
  cat(sprintf("  REML score      : %.3f\n", model$gcv.ubre))
  cat(sprintf("  AIC             : %.2f\n", AIC(model)))
  cat("\n  Smooth terms (edf / Ref.df / F / p-value):\n")
  print(round(s$s.table, 4))
}
report_model(m_additive, "m_additive", t_add[["elapsed"]])
report_model(m_interact, "m_interact", t_int[["elapsed"]])

cat("\n==============================================================\n")
cat(" MODEL COMPARISON\n")
cat("==============================================================\n")
cat(sprintf("AIC(m_additive) = %.2f\nAIC(m_interact) = %.2f\nDelta AIC      = %.2f  (negative favours interact)\n",
            AIC(m_additive), AIC(m_interact), AIC(m_interact) - AIC(m_additive)))
cat("\nanova(m_additive, m_interact, test='F'):\n")
print(anova(m_additive, m_interact, test = "F"))

cat("\n==============================================================\n")
cat(" gam.check(m_additive)  (text output)\n")
cat("==============================================================\n")
cat(paste(gam_check_text, collapse = "\n"), "\n")

cat("\nPlots written to ", OUT_DIR, ":\n", sep = "")
print(list.files(OUT_DIR))

cat("\nDone.\n")
