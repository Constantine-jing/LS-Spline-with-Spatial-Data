## ============================================================
## scripts/build_sims_seeds_2_3_4.R
##
## Build n=1500 sim RDS files for seeds 2, 3, 4 — for the
## multi-seed robustness runs of the null-interaction finding.
##
## Mirrors scripts/run_n1500_prep.R (which built seed=1) but
## loops over 2:4 and refuses to overwrite an existing RDS.
##
## Outputs:
##   output/sim_medsat_london_asthma_n1500_seed2_v2.rds
##   output/sim_medsat_london_asthma_n1500_seed3_v2.rds
##   output/sim_medsat_london_asthma_n1500_seed4_v2.rds
##
## At the end, prints a side-by-side comparison vs the existing
## seed=1 sim so we can confirm: similar borough coverage, similar
## raw y distribution, all 33 boroughs in each sim.
## ============================================================

SEEDS <- 2:4

skew <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

env <- new.env(parent = globalenv())
sys.source("R/prepare_medsat_london.R", chdir = TRUE, envir = env)

build_one <- function(seed) {
  out_rds <- sprintf("output/sim_medsat_london_asthma_n1500_seed%d_v2.rds", seed)
  if (file.exists(out_rds)) {
    stop(sprintf(
      "Output already exists: %s\n  Delete it first to regenerate.", out_rds))
  }
  cat(sprintf("\n=== Building sim seed %d ===\n", seed))
  t0 <- Sys.time()
  sim <- env$prepare_medsat_london(
    medsat_path       = "data/medsat/2019_spatial_raw_master.csv",
    imd_path          = "data/medsat/imd_2019_lsoa21.csv",
    outcome           = "asthma",
    covariates        = c("NO2", "NDVI", "IMD"),
    n_subsample       = 1500L,
    seed              = seed,
    drop_ndvi_ceiling = TRUE,
    verbose           = FALSE
  )
  cat(sprintf("  prep time: %.2f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  saveRDS(sim, out_rds)
  cat(sprintf("  saved: %s  (%.1f KB)\n", out_rds, file.info(out_rds)$size / 1024))

  # ---- schema guard (same as in run_n1500_prep.R) ----
  stopifnot(sim$n == 1500L, sim$p == 3L)
  stopifnot(identical(sim$truth_params$nu, 1.0))
  stopifnot(is.na(sim$truth_params$rho), is.na(sim$truth_params$sigma2))
  stopifnot(setequal(sim$int_keys, c("1_2", "1_3", "2_3")))
  X <- as.matrix(sim$data[, paste0("X", seq_len(sim$p))])
  stopifnot(all(X >= 0 - 1e-12), all(X <= 1 + 1e-12))

  sim
}

sims <- lapply(SEEDS, build_one)
names(sims) <- as.character(SEEDS)

# ============================================================
# Per-seed summaries
# ============================================================
cat("\n==============================================================\n")
cat(" Per-seed summaries\n")
cat("==============================================================\n")
for (k in seq_along(SEEDS)) {
  s   <- SEEDS[k]
  sim <- sims[[k]]
  ub  <- sort(unique(sim$borough))
  btab <- sort(table(sim$borough), decreasing = TRUE)
  y    <- sim$data$y

  cat(sprintf("\nSeed %d:\n", s))
  cat(sprintf("  n = %d   p = %d   nu = %g\n", sim$n, sim$p, sim$truth_params$nu))
  cat(sprintf("  boroughs: %d / 33 unique   min/borough = %d   max/borough = %d   mean = %.1f\n",
              length(ub), min(btab), max(btab), mean(btab)))
  cat(sprintf("  y (%s): mean=%.4f sd=%.4f median=%.4f range=[%.4f, %.4f] skew=%.4f\n",
              sim$meta$response_transform,
              mean(y), sd(y), median(y), min(y), max(y), skew(y)))
  cat(sprintf("  NDVI ceiling dropped: %d   |   NA drops (total): %d\n",
              sim$meta$ndvi_ceiling_dropped, sim$meta$na_drops$total))
  cat("  X ranges:\n")
  for (j in seq_len(sim$p)) {
    xs <- sim$X_scale[[j]]
    cat(sprintf("    X%d (%s): [%.4g, %.4g]\n", j, xs$name, xs$lo, xs$hi))
  }
  # Confirm all 33 ONS boroughs present
  ons33 <- c(
    "Barking and Dagenham","Barnet","Bexley","Brent","Bromley","Camden",
    "City of London","Croydon","Ealing","Enfield","Greenwich","Hackney",
    "Hammersmith and Fulham","Haringey","Harrow","Havering","Hillingdon",
    "Hounslow","Islington","Kensington and Chelsea","Kingston upon Thames",
    "Lambeth","Lewisham","Merton","Newham","Redbridge","Richmond upon Thames",
    "Southwark","Sutton","Tower Hamlets","Waltham Forest","Wandsworth",
    "Westminster")
  missing_b <- setdiff(ons33, ub)
  if (length(missing_b) > 0L) {
    cat(sprintf("  WARNING: missing boroughs (%d): %s\n",
                length(missing_b), paste(missing_b, collapse = ", ")))
  } else {
    cat("  all 33 ONS-canonical boroughs represented.\n")
  }
}

# ============================================================
# Side-by-side vs existing seed=1 sim
# ============================================================
SIM1 <- "output/sim_medsat_london_asthma_n1500_seed1_v2.rds"
if (file.exists(SIM1)) {
  cat("\n==============================================================\n")
  cat(" Side-by-side comparison (incl. seed=1 from disk)\n")
  cat("==============================================================\n")
  all_sims <- c(list("1" = readRDS(SIM1)), sims)
  rows <- do.call(rbind, lapply(names(all_sims), function(nm) {
    sim <- all_sims[[nm]]
    y   <- sim$data$y
    data.frame(
      seed       = as.integer(nm),
      n          = sim$n,
      boroughs   = length(unique(sim$borough)),
      y_mean     = round(mean(y), 4),
      y_sd       = round(sd(y), 4),
      y_median   = round(median(y), 4),
      y_min      = round(min(y), 4),
      y_max      = round(max(y), 4),
      y_skew     = round(skew(y), 3),
      ndvi_drop  = sim$meta$ndvi_ceiling_dropped,
      stringsAsFactors = FALSE
    )
  }))
  rownames(rows) <- NULL
  print(rows, row.names = FALSE)
}

cat("\n==============================================================\n")
cat(" Done. Three new sim RDS files saved under output/.\n")
cat("==============================================================\n")
