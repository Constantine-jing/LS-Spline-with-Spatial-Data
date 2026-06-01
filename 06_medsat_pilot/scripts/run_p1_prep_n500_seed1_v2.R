## ============================================================
## scripts/run_p1_prep_n500_seed1_v2.R
##
## Step B: run prepare_medsat_london_test() (smoke).
## Step C: rerun real prep with drop_ndvi_ceiling = TRUE; save to
##         output/sim_medsat_london_asthma_n500_seed1_v2.rds.
##         The v1 RDS is left untouched.
## Step D: report diagnostics.
##
## Does NOT fit a model.
## ============================================================

OUT_RDS_V2 <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"
stopifnot(file.exists("output/sim_medsat_london_asthma_n500_seed1.rds"))  # v1 must exist
stopifnot(!file.exists(OUT_RDS_V2))                                       # v2 must not yet

env <- new.env(parent = globalenv())
sys.source("R/prepare_medsat_london.R", chdir = TRUE, envir = env)

skew <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

# ------------------------------------------------------------
# B. Smoke test
# ------------------------------------------------------------
cat("==============================================================\n")
cat(" B. Smoke test: prepare_medsat_london_test()\n")
cat("==============================================================\n")
sim_synth <- env$prepare_medsat_london_test()
cat(sprintf("\nSmoke-test sim$meta$ndvi_ceiling_dropped = %d (expected 99)\n",
            sim_synth$meta$ndvi_ceiling_dropped))
cat(sprintf("Smoke-test sim$meta$na_drops$total      = %d (expected 0)\n",
            sim_synth$meta$na_drops$total))
cat(sprintf("Smoke-test nrow(sim$meta$low_y_rows)    = %d\n",
            nrow(sim_synth$meta$low_y_rows)))

# ------------------------------------------------------------
# C. Real prep v2
# ------------------------------------------------------------
cat("\n==============================================================\n")
cat(" C. Real prep v2 (n=500 seed=1, drop_ndvi_ceiling=TRUE)\n")
cat("==============================================================\n")
sim <- env$prepare_medsat_london(
  medsat_path       = "data/medsat/2019_spatial_raw_master.csv",
  imd_path          = "data/medsat/imd_2019_lsoa21.csv",
  outcome           = "asthma",
  covariates        = c("NO2", "NDVI", "IMD"),
  n_subsample       = 500L,
  seed              = 1L,
  drop_ndvi_ceiling = TRUE
)
saveRDS(sim, OUT_RDS_V2)

# ------------------------------------------------------------
# D. Report
# ------------------------------------------------------------
cat("\n==============================================================\n")
cat(" D. Diagnostics for v2 sim\n")
cat("==============================================================\n")

cat(sprintf("\nSaved file:  %s   (%.1f KB)\n",
            OUT_RDS_V2, file.info(OUT_RDS_V2)$size / 1024))
cat(sprintf("V1 (kept):   output/sim_medsat_london_asthma_n500_seed1.rds   (%.1f KB)\n",
            file.info("output/sim_medsat_london_asthma_n500_seed1.rds")$size / 1024))

cat(sprintf("\nsim$n = %d   sim$p = %d\n", sim$n, sim$p))

# borough coverage
ub <- sort(unique(sim$borough))
cat(sprintf("Unique boroughs in subsample: %d\n", length(ub)))
btab <- sort(table(sim$borough), decreasing = TRUE)
cat("Borough counts in subsample (top 10 / bottom 5):\n")
print(head(btab, 10))
cat("  ...\n")
print(tail(btab, 5))
cat(sprintf("  min/borough = %d, max/borough = %d, mean = %.2f\n",
            min(btab), max(btab), mean(btab)))

# y diagnostics
y <- sim$data$y
cat("\nResponse y (transformed=", sim$meta$response_transform, "):\n", sep = "")
cat(sprintf("  n=%d  mean=%.4f  sd=%.4f  median=%.4f  min=%.4f  max=%.4f  skew=%.4f\n",
            length(y), mean(y), sd(y), median(y), min(y), max(y), skew(y)))

# X2 (NDVI after scaling) sanity check
x2 <- sim$data$X2
cat("\nsim$data$X2 (NDVI after [0,1] min-max scaling):\n")
cat(sprintf("  min = %.6f  max = %.6f\n", min(x2), max(x2)))
cat(sprintf("  # values >= 0.99 (scaled): %d\n", sum(x2 >= 0.99)))
cat(sprintf("  raw NDVI lo = %.4f  hi = %.4f  (X_scale[[2]])\n",
            sim$X_scale[[2]]$lo, sim$X_scale[[2]]$hi))
cat(sprintf("  # raw NDVI rows == 1.000 in subsample: %d (expect 0)\n",
            sum(sim$data$X2 * (sim$X_scale[[2]]$hi - sim$X_scale[[2]]$lo) +
                sim$X_scale[[2]]$lo == 1.0)))

# meta$ndvi_ceiling_dropped
cat(sprintf("\nsim$meta$ndvi_ceiling_dropped = %d  (NDVI=1.0 rows removed pre-subsample)\n",
            sim$meta$ndvi_ceiling_dropped))

# meta$na_drops
nad <- sim$meta$na_drops
cat(sprintf("\nsim$meta$na_drops$total = %d\n", nad$total))
cat("sim$meta$na_drops$by_column (non-zero rows):\n")
print(nad$by_column[nad$by_column > 0L])
cat("sim$meta$na_drops$by_borough (top 5 by % dropped):\n")
print(head(nad$by_borough, 5))

# meta$low_y_rows
ly <- sim$meta$low_y_rows
cat(sprintf("\nsim$meta$low_y_rows: nrow = %d   (subsample rows with y_raw < 1.0)\n",
            nrow(ly)))
if (nrow(ly) > 0L) {
  cat("Breakdown by borough:\n")
  print(sort(table(ly$borough), decreasing = TRUE))
}

# truth_params
cat("\nsim$truth_params:\n"); print(sim$truth_params)
stopifnot(identical(sim$truth_params$nu, 1.0))

cat("\n==============================================================\n")
cat(" Done. No model fitted.\n")
cat("==============================================================\n")
