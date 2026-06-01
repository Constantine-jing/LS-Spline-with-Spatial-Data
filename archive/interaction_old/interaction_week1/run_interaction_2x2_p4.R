# ============================================================
# run_interaction_2x2_p4.R
#
# HEADLINE INTERACTION RUN: p = 4 only, both phases, at n = 1000
#
# Purpose:
#   (1) Reproduce the 2x2 factorial experiment for p = 4 at the same
#       sample size (n = 1000) used for the main-effects Bayesian
#       simulations in Chapter 1 Section 5. This is the chapter
#       headline for the interaction extension.
#   (2) Run across multiple data-seed groups (seed testing) so we can
#       check whether results are stable across random data draws.
#   (3) Output clearer 2x2-layout shared-scale heatmaps alongside the
#       existing 3D wireframe plots.
#
# Configuration:
#   N_OBS       = 1000        (bumped from 300 -> 1000)
#   M_KNOTS     = 8           (same as v2)
#   N_ITER      = 3000
#   N_BURN      = 1000
#   C_INT       = 1.5
#   SEED_GROUPS = length of SEED_OFFSETS below
#                 each entry shifts the default seed tuple by that offset,
#                 producing a fully distinct NULL/INT/MCMC seed triple
#   PHASES      = c(1, 2)     (Phase 1 = X1xX2 only, Phase 2 = all 6 pairs)
#
# Runtime estimate:
#   n = 1000 -> Cholesky factor is ~37x slower per iter than n = 300.
#   Empirically ~2.5-3 hours per cell at n = 1000 on Hellbender.
#   Full grid: 4 cells x 2 phases x N seed groups =
#     N_SEEDS = 1 -> ~22 hours  (fits one SBATCH job under 48 hr)
#     N_SEEDS = 2 -> ~44 hours  (needs to be split into two jobs)
#
# Split strategy (recommended):
#   Run this script twice with PHASE_RESTRICT controlling which phase
#   is processed, so each SBATCH job handles only one phase:
#     Job 1: Rscript run_interaction_2x2_p4.R 1   # phase 1 only
#     Job 2: Rscript run_interaction_2x2_p4.R 2   # phase 2 only
#   Each job takes ~11-22 hours depending on N_SEEDS.
#   With no argument, runs both phases sequentially (single job).
#
# Output:
#   interaction_2x2_p4_summary.csv
#       one row per (phase, seed group, cell). With 2 phases x N_SEEDS
#       seed groups x 4 cells = 8*N_SEEDS rows total.
#   interaction_2x2_p4_results.rds
#   For each (phase, seed group):
#       interaction_2x2_p4_phase{ph}_seed{g}_marginals.pdf
#       interaction_2x2_p4_phase{ph}_seed{g}_surfaces_2x2.pdf  <- NEW shared-scale layout
#       interaction_2x2_p4_phase{ph}_seed{g}_surfaces_3d.pdf
#       interaction_2x2_p4_phase{ph}_seed{g}_diagnostics.pdf
#
# Assumes the existing functions from run_interaction_2x2_v2.R are sourced:
#   simulate_unified(), build_design_unified(), fit_one_cell(),
#   summarise_cell_unified(), run_2x2(), print_2x2_table_v2(),
#   plot_marginals_v2(), plot_diagnostics_v2(), plot_surfaces_3d_v2(),
#   results_to_dataframe_v2()
#
# and the companion files these depend on:
#   spatial_utils.R, ls_basis.R, ls_interaction.R, gibbs_interaction.R
# ============================================================


# ============================================================
# CLI arg: phase restriction ("1", "2", or omitted for both)
# ============================================================
PHASE_RESTRICT <- 2L

# ============================================================
# Configuration
# ============================================================
N_OBS   <- 1000
M_KNOTS <- 8
N_ITER  <- 3000
N_BURN  <- 1000
C_INT   <- 1.5
P_FIX   <- 4L

# Seed-testing groups.
# Each offset produces a fully distinct (data_null, data_int, mcmc) seed triple
# by adding the offset to the default seed scheme in run_2x2().
# To add more seed groups, extend this vector; to run just one, leave as c(0).
SEED_OFFSETS <- c(0L)   

PHASES <- if (is.na(PHASE_RESTRICT)) c(1L, 2L) else as.integer(PHASE_RESTRICT)
stopifnot(all(PHASES %in% c(1L, 2L)))


# ============================================================
# Source existing v2 machinery
# ============================================================
# This assumes run_interaction_2x2_v2.R is in the same directory.
# It provides: simulate_unified, build_design_unified, fit_one_cell,
# summarise_cell_unified, run_2x2, print_2x2_table_v2,
# plot_marginals_v2, plot_surfaces_v2, plot_surfaces_3d_v2,
# plot_diagnostics_v2, results_to_dataframe_v2, compute_waic, etc.
#
# We source run_interaction_2x2_v2.R DEFENSIVELY: we only need its
# function definitions, not its own driver block (which runs the full
# p=2,3,4 sweep at n=300). Strategy: source it with interactive() TRUE
# by wrapping in a temporary override so the final `if (!interactive())`
# driver block is skipped.
.restore_interactive <- interactive
assign("interactive", function() TRUE, envir = globalenv())
source("run_interaction_2x2_v2.R", local = FALSE)
assign("interactive", .restore_interactive, envir = globalenv())
rm(.restore_interactive)

# Also source the underlying dependencies (run_2x2 needs these at call time)
source("spatial_utils.R")
source("ls_basis.R")
source("ls_interaction.R")
source("gibbs_interaction.R")


# ============================================================
# NEW: shared-scale 2x2 layout heatmap plotter
#
# Produces a single-page 2x2 PDF with:
#   Top-left:  TRUE f_12
#   Top-right: POSTERIOR MEAN f_12 (Cell D)
#   Bot-left:  POSTERIOR MEAN f_12 (Cell B, truth = 0)
#   Bot-right: RESIDUAL (Cell D: fitted - truth)
#
# All four panels use the SAME symmetric z-scale and palette so the
# intensity of each color means the same numerical value across panels.
# This is the fix for the "hard to compare" feedback.
#
# A second page in the PDF shows the band-width uncertainty map on its
# own (single-panel) scale, since band width is a different quantity
# and would not be comparable on the shared scale anyway.
# ============================================================
plot_surfaces_2x2_shared <- function(res_obj, outfile) {
  s       <- res_obj$settings
  results <- res_obj$results
  truth_int <- res_obj$truth_grid_int$f12

  surfD <- results$int_M1$surface12
  surfB <- results$null_M1$surface12

  pdf(outfile, width = 10, height = 10)

  if (is.null(surfD)) {
    plot.new(); text(0.5, 0.5, "No interaction surface stored")
    dev.off(); return(invisible(NULL))
  }

  x1 <- surfD$x1; x2 <- surfD$x2

  # ---- SHARED symmetric z-scale across ALL four panels ----
  z_all <- c(
    as.vector(truth_int),
    as.vector(surfD$post_mean),
    if (!is.null(surfB)) as.vector(surfB$post_mean) else numeric(0),
    as.vector(surfD$post_mean - truth_int)
  )
  z_lim <- max(abs(z_all), na.rm = TRUE)
  z_breaks <- seq(-z_lim, z_lim, length.out = 51)
  pal <- hcl.colors(50, "Blue-Red 2")

  # --------------- PAGE 1: 2x2 shared-scale comparison ---------------
  par(mfrow = c(2, 2),
      mar = c(4.2, 4.2, 3.2, 5.5),    # extra right margin for color legend
      oma = c(0, 0, 3, 0),
      cex.main = 1.2, cex.lab = 1.1, cex.axis = 1.0)

  # (a) Truth
  image(x1, x2, truth_int, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("(a) TRUE f_12\nrange [%.2f, %.2f]",
                       min(truth_int), max(truth_int)))
  contour(x1, x2, truth_int, add = TRUE, lwd = 1.0, labcex = 0.75)
  .add_colorbar(z_breaks, pal)

  # (b) Cell D posterior mean (recovery)
  image(x1, x2, surfD$post_mean, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("(b) POSTERIOR MEAN (Cell D: INT + M1)\nRMSE = %.3f, cov = %.2f",
                       results$int_M1$rmse_f12, results$int_M1$cov_f12))
  contour(x1, x2, surfD$post_mean, add = TRUE, lwd = 1.0, labcex = 0.75)
  .add_colorbar(z_breaks, pal)

  # (c) Cell B posterior mean (null shrinkage)
  if (!is.null(surfB)) {
    image(x1, x2, surfB$post_mean, col = pal, breaks = z_breaks,
          xlab = "X1", ylab = "X2",
          main = sprintf("(c) POSTERIOR MEAN (Cell B: NULL + M1, truth=0)\nRMSE = %.4f, cov of 0 = %.2f",
                         results$null_M1$rmse_f12, results$null_M1$cov_f12))
    contour(x1, x2, surfB$post_mean, add = TRUE, lwd = 1.0, labcex = 0.75)
  } else {
    plot.new()
    title(main = "(c) Cell B surface not available")
  }
  .add_colorbar(z_breaks, pal)

  # (d) Residual (Cell D fitted - truth)
  resid_surf <- surfD$post_mean - truth_int
  image(x1, x2, resid_surf, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("(d) RESIDUAL: Cell D - truth\nrange [%.2f, %.2f]",
                       min(resid_surf), max(resid_surf)))
  contour(x1, x2, resid_surf, add = TRUE, lwd = 1.0, labcex = 0.75)
  .add_colorbar(z_breaks, pal)

  # Overall title
  mtext(sprintf("Interaction f_12: p = %d, phase %d  (shared z-scale, z_lim = %.2f)",
                s$p, s$phase, z_lim),
        outer = TRUE, cex = 1.3, font = 2, line = 1)

  # --------------- PAGE 2: Uncertainty (own scale) ---------------
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.5, 5.5), oma = c(0, 0, 0, 0),
      cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

  pal_unc <- hcl.colors(50, "YlOrRd", rev = TRUE)
  bw_max <- max(surfD$band_width, na.rm = TRUE)
  bw_breaks <- seq(0, bw_max, length.out = 51)

  image(x1, x2, surfD$band_width, col = pal_unc, breaks = bw_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("UNCERTAINTY: 95%% band width  (Cell D, p=%d, phase %d)\nmean = %.3f, max = %.3f",
                       s$p, s$phase, mean(surfD$band_width), max(surfD$band_width)))
  contour(x1, x2, surfD$band_width, add = TRUE, lwd = 1.0, labcex = 0.75)
  .add_colorbar(bw_breaks, pal_unc)

  dev.off()
  cat(sprintf("  2x2 shared-scale surfaces: %s\n", outfile))
}


# ============================================================
# Helper: minimal in-plot color legend bar
# Adds a narrow vertical colorbar on the right side of the current panel,
# using the current user coordinates. Call immediately after an image().
# ============================================================
.add_colorbar <- function(breaks, palette) {
  # Save current coords, switch to figure-fraction for the bar
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  usr <- par("usr")     # c(xmin, xmax, ymin, ymax) in plot coords
  xr  <- usr[2] - usr[1]
  yr  <- usr[4] - usr[3]

  # Place bar in a strip to the right of the plot region, in plot coords.
  # We rely on par(xpd = NA) so we can draw outside the plot box.
  par(xpd = NA)
  xl <- usr[2] + 0.04 * xr
  xr_ <- usr[2] + 0.09 * xr
  y_bot <- usr[3]
  y_top <- usr[4]
  n <- length(palette)
  ys <- seq(y_bot, y_top, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(xl, ys[i], xr_, ys[i + 1], col = palette[i], border = NA)
  }
  rect(xl, y_bot, xr_, y_top, border = "black", lwd = 0.5)

  # Tick labels: show min, 0, max (or endpoints if 0 not in range)
  br_min <- breaks[1]; br_max <- tail(breaks, 1)
  tick_vals <- if (br_min < 0 && br_max > 0) c(br_min, 0, br_max) else c(br_min, br_max)
  tick_ys <- y_bot + (tick_vals - br_min) / (br_max - br_min) * (y_top - y_bot)
  for (k in seq_along(tick_vals)) {
    lines(c(xr_, xr_ + 0.01 * xr), c(tick_ys[k], tick_ys[k]), lwd = 0.6)
    text(xr_ + 0.015 * xr, tick_ys[k],
         sprintf("%.2f", tick_vals[k]),
         adj = c(0, 0.5), cex = 0.75)
  }
}


# ============================================================
# Driver
# ============================================================
cat("==============================================================\n")
cat(" p = 4 HEADLINE RUN at n =", N_OBS, "\n")
cat(" phases:        ", paste(PHASES, collapse = ", "), "\n")
cat(" seed groups:   ", length(SEED_OFFSETS), " (offsets:",
    paste(SEED_OFFSETS, collapse = ", "), ")\n")
cat(" M =", M_KNOTS, " n_iter =", N_ITER, " n_burn =", N_BURN, " c_int =", C_INT, "\n")
cat("==============================================================\n")

t_total <- proc.time()

all_res <- list()
for (phase_ in PHASES) {
  for (g in seq_along(SEED_OFFSETS)) {
    offset <- SEED_OFFSETS[g]
    tag <- sprintf("phase%d_seed%d", phase_, g)
    cat(sprintf("\n###############  p = 4, phase %d, seed group %d (offset = %d)  ###############\n",
                phase_, g, offset))

    # Build seeds. run_2x2 by default uses 100*p+phase etc.; we add offset
    # to get fully distinct data realizations between seed groups.
    data_seed_null <- 100L * P_FIX + phase_ + offset
    data_seed_int  <- 200L * P_FIX + phase_ + offset
    mcmc_seed      <- 300L * P_FIX + phase_ + offset

    res <- run_2x2(
      p = P_FIX, phase = phase_,
      n = N_OBS, M = M_KNOTS,
      n_iter = N_ITER, n_burn = N_BURN, c_int = C_INT,
      data_seed_null = data_seed_null,
      data_seed_int  = data_seed_int,
      mcmc_seed      = mcmc_seed
    )
    # Record seed metadata on the settings block for downstream reporting
    res$settings$seed_group <- g
    res$settings$seed_offset <- offset
    res$settings$data_seed_null <- data_seed_null
    res$settings$data_seed_int  <- data_seed_int
    res$settings$mcmc_seed      <- mcmc_seed

    all_res[[tag]] <- res

    # Plot outputs for this (phase, seed-group)
    base <- sprintf("interaction_2x2_p4_%s", tag)
    plot_marginals_v2(res,       outfile = paste0(base, "_marginals.pdf"))
    plot_surfaces_v2(res,        outfile = paste0(base, "_surfaces_flat.pdf"))     # old style
    plot_surfaces_2x2_shared(res, outfile = paste0(base, "_surfaces_2x2.pdf"))     # NEW shared-scale layout
    plot_surfaces_3d_v2(res,     outfile = paste0(base, "_surfaces_3d.pdf"))
    plot_diagnostics_v2(res,     outfile = paste0(base, "_diagnostics.pdf"))

    # Print summary table for the console log
    print_2x2_table_v2(res$results, p = P_FIX, phase = phase_)
  }
}

total_elapsed <- as.numeric((proc.time() - t_total)["elapsed"])
cat(sprintf("\n##########  TOTAL: %.1f sec (%.2f hours)  ##########\n",
            total_elapsed, total_elapsed / 3600))

# ----- Save artifacts -----
saveRDS(list(all = all_res, total_sec = total_elapsed,
             config = list(N_OBS = N_OBS, M_KNOTS = M_KNOTS,
                           N_ITER = N_ITER, N_BURN = N_BURN,
                           C_INT = C_INT, P_FIX = P_FIX,
                           SEED_OFFSETS = SEED_OFFSETS,
                           PHASES = PHASES)),
        file = "interaction_2x2_p4_results.rds")

# Build summary CSV. We reuse results_to_dataframe_v2() and add the seed-group
# metadata so each row is uniquely identified by (phase, seed_group, cell).
df_list <- list()
for (tag in names(all_res)) {
  res <- all_res[[tag]]
  df  <- results_to_dataframe_v2(res)
  df$seed_group     <- res$settings$seed_group
  df$seed_offset    <- res$settings$seed_offset
  df$data_seed_null <- res$settings$data_seed_null
  df$data_seed_int  <- res$settings$data_seed_int
  df$mcmc_seed      <- res$settings$mcmc_seed
  df_list[[tag]]    <- df
}
df_all <- do.call(rbind, df_list)
rownames(df_all) <- NULL
write.csv(df_all, "interaction_2x2_p4_summary.csv", row.names = FALSE)

cat("\nCSV saved to interaction_2x2_p4_summary.csv\n")
cat("RDS saved to interaction_2x2_p4_results.rds\n")
