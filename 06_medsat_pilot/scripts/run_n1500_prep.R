## ============================================================
## scripts/run_n1500_prep.R
##
## Build the n=1500 v2 sim object for the Hellbender scale-up run.
##
##   prepare_medsat_london(
##     n_subsample = 1500, seed = 1, drop_ndvi_ceiling = TRUE, ...
##   )
##   -> output/sim_medsat_london_asthma_n1500_seed1_v2.rds
##
## Does NOT fit a model. Designed to be run once locally to produce
## the sim RDS; the RDS is then transferred to Hellbender alongside
## scripts/run_n1500_allpairs.R.
## ============================================================

OUT_RDS <- "output/sim_medsat_london_asthma_n1500_seed1_v2.rds"

if (file.exists(OUT_RDS)) {
  stop(sprintf(
    "Output already exists: %s\n  Delete it first if you want to regenerate.",
    OUT_RDS))
}

skew <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

env <- new.env(parent = globalenv())
sys.source("R/prepare_medsat_london.R", chdir = TRUE, envir = env)

cat("==============================================================\n")
cat(" Building n=1500 sim (seed=1, drop_ndvi_ceiling=TRUE)\n")
cat("==============================================================\n")
t0 <- Sys.time()
sim <- env$prepare_medsat_london(
  medsat_path       = "data/medsat/2019_spatial_raw_master.csv",
  imd_path          = "data/medsat/imd_2019_lsoa21.csv",
  outcome           = "asthma",
  covariates        = c("NO2", "NDVI", "IMD"),
  n_subsample       = 1500L,
  seed              = 1L,
  drop_ndvi_ceiling = TRUE
)
cat(sprintf("Prep time: %.2f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

saveRDS(sim, OUT_RDS)
sz_kb <- file.info(OUT_RDS)$size / 1024
cat(sprintf("\nSaved: %s  (%.1f KB)\n", OUT_RDS, sz_kb))

cat("\n==============================================================\n")
cat(" Diagnostics\n")
cat("==============================================================\n")

cat(sprintf("\nsim$n = %d   sim$p = %d   nu = %g\n",
            sim$n, sim$p, sim$truth_params$nu))

# ---- borough coverage ----
ub <- sort(unique(sim$borough))
btab <- sort(table(sim$borough), decreasing = TRUE)
cat(sprintf("Unique boroughs: %d / 33\n", length(ub)))
cat(sprintf("Borough counts: min=%d  max=%d  mean=%.1f\n",
            min(btab), max(btab), mean(btab)))
cat("Top 5 boroughs:\n")
print(head(btab, 5))
cat("Bottom 5 boroughs:\n")
print(tail(btab, 5))

# ---- response ----
y <- sim$data$y
cat(sprintf("\ny (%s): n=%d  mean=%.4f  sd=%.4f  median=%.4f\n",
            sim$meta$response_transform,
            length(y), mean(y), sd(y), median(y)))
cat(sprintf("  range=[%.4f, %.4f]  skew=%.4f\n",
            min(y), max(y), skew(y)))

# ---- covariates ----
cat("\nCovariate original ranges:\n")
for (j in seq_len(sim$p)) {
  xs <- sim$X_scale[[j]]
  cat(sprintf("  X%d (%s): [%.4f, %.4f]\n", j, xs$name, xs$lo, xs$hi))
}

# ---- scaled X sanity ----
X <- as.matrix(sim$data[, paste0("X", seq_len(sim$p))])
cat("\nScaled X range check (all must be [0,1]):\n")
for (j in seq_len(sim$p)) {
  cat(sprintf("  X%d: min=%.6f  max=%.6f\n", j, min(X[,j]), max(X[,j])))
}

# ---- coords ----
cat(sprintf("\nCoords (unit square): lon=[%.4f,%.4f]  lat=[%.4f,%.4f]\n",
            min(sim$data$lon), max(sim$data$lon),
            min(sim$data$lat), max(sim$data$lat)))

# ---- NDVI ceiling ----
cat(sprintf("\nNDVI ceiling rows dropped (Sentinel-2 artifact): %d\n",
            sim$meta$ndvi_ceiling_dropped))
cat(sprintf("Raw NDVI hi after drop: %.6f  (X_scale[[2]]$hi, must be < 1.0)\n",
            sim$X_scale[[2]]$hi))

# ---- NA drops ----
nad <- sim$meta$na_drops
cat(sprintf("\nNA drops total: %d\n", nad$total))
nonzero_cols <- nad$by_column[nad$by_column > 0L]
if (length(nonzero_cols) > 0L) {
  cat("  Non-zero NA counts by column:\n")
  print(nonzero_cols)
} else {
  cat("  (no NAs in required columns)\n")
}

# ---- low-y rows ----
ly <- sim$meta$low_y_rows
cat(sprintf("\nLow-y rows (y_raw < 1.0): %d\n", nrow(ly)))
if (nrow(ly) > 0L) {
  cat("  Borough breakdown:\n")
  print(sort(table(ly$borough), decreasing = TRUE))
}

# ---- schema guard ----
stopifnot(sim$n == 1500L, sim$p == 3L)
stopifnot(identical(sim$truth_params$nu, 1.0))
stopifnot(is.na(sim$truth_params$rho), is.na(sim$truth_params$sigma2))
stopifnot(setequal(sim$int_keys, c("1_2", "1_3", "2_3")))
stopifnot(all(X >= 0 - 1e-12), all(X <= 1 + 1e-12))

# ---- compare with n=500 object ----
SIM500 <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"
if (file.exists(SIM500)) {
  sim500 <- readRDS(SIM500)
  cat("\n==============================================================\n")
  cat(" Comparison vs n=500 sim\n")
  cat("==============================================================\n")
  for (j in seq_len(sim$p)) {
    xs1500 <- sim$X_scale[[j]]
    xs500  <- sim500$X_scale[[j]]
    cat(sprintf("  %s: n=1500 range=[%.4f,%.4f]  n=500 range=[%.4f,%.4f]\n",
                xs1500$name, xs1500$lo, xs1500$hi, xs500$lo, xs500$hi))
  }
  cat(sprintf("  n=1500 y: mean=%.4f sd=%.4f\n", mean(sim$data$y), sd(sim$data$y)))
  cat(sprintf("  n=500  y: mean=%.4f sd=%.4f\n", mean(sim500$data$y), sd(sim500$data$y)))
  cat(sprintf("  n=1500 NDVI ceiling dropped: %d\n", sim$meta$ndvi_ceiling_dropped))
  cat(sprintf("  n=500  NDVI ceiling dropped: %d\n", sim500$meta$ndvi_ceiling_dropped))
}

cat("\n==============================================================\n")
cat(" Done. No model fitted.\n")
cat(sprintf(" RDS: %s\n", OUT_RDS))
cat("==============================================================\n")
