## ============================================================
## scripts/run_p1_prep_n500_seed1.R
##
## Run the smoke test for prepare_medsat_london(), then run the
## real pilot prep (n=500 London LSOAs, seed 1, P1 covariates), save
## the sim object, and print diagnostics needed before fitting.
##
## DOES NOT fit any model.
## ============================================================

OUT_RDS  <- "output/sim_medsat_london_asthma_n500_seed1.rds"
OUT_HIST <- "output/plots/y_histogram_p1.pdf"

dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)

# Source the prep file with chdir=TRUE so its relative `source("spatial_utils.R")`
# resolves inside the R/ directory.
env <- new.env(parent = globalenv())
sys.source("R/prepare_medsat_london.R", chdir = TRUE, envir = env)

# Skewness helper (no extra packages required).
skew <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

# ------------------------------------------------------------
# 1. Smoke test
# ------------------------------------------------------------
cat("==============================================================\n")
cat(" 1. Running prepare_medsat_london_test() ...\n")
cat("==============================================================\n")
sim_synth <- env$prepare_medsat_london_test()
cat(sprintf("Synthetic subsample boroughs represented: %d/33\n",
            length(unique(sim_synth$borough))))
print(sort(unique(sim_synth$borough)))

# ------------------------------------------------------------
# 2. Real prep
# ------------------------------------------------------------
cat("\n==============================================================\n")
cat(" 2. Running real prep on MEDSAT 2019 London (n=500, seed=1) ...\n")
cat("==============================================================\n")
sim <- env$prepare_medsat_london(
  medsat_path = "data/medsat/2019_spatial_raw_master.csv",
  imd_path    = "data/medsat/imd_2019_lsoa21.csv",
  outcome     = "asthma",
  covariates  = c("NO2", "NDVI", "IMD"),
  n_subsample = 500L,
  seed        = 1L
)

saveRDS(sim, OUT_RDS)
cat(sprintf("\nSaved sim object: %s (%.1f KB)\n",
            OUT_RDS, file.info(OUT_RDS)$size / 1024))

# ------------------------------------------------------------
# 3. Report
# ------------------------------------------------------------
cat("\n==============================================================\n")
cat(" 3. Sim diagnostics\n")
cat("==============================================================\n")

cat(sprintf("\nsim$n = %d   sim$p = %d\n", sim$n, sim$p))

# borough coverage
boroughs <- sim$borough
ub <- sort(unique(boroughs))
cat(sprintf("Unique boroughs in subsample: %d / 33\n", length(ub)))
cat("Borough counts in subsample (n=500):\n")
btab <- sort(table(boroughs), decreasing = TRUE)
print(btab)
cat(sprintf("  min/borough = %d, max/borough = %d, mean = %.2f\n",
            min(btab), max(btab), mean(btab)))

# y diagnostics
y <- sim$data$y
cat("\nResponse y (transformed=", sim$meta$response_transform, "):\n", sep = "")
cat(sprintf("  n      = %d\n", length(y)))
cat(sprintf("  mean   = %.4f\n", mean(y)))
cat(sprintf("  sd     = %.4f\n", sd(y)))
cat(sprintf("  median = %.4f\n", median(y)))
cat(sprintf("  min    = %.4f\n", min(y)))
cat(sprintf("  max    = %.4f\n", max(y)))
cat(sprintf("  skewness = %.4f  (Gaussian: 0; |.|>1 ~ strongly skewed)\n", skew(y)))

# also log-transform skewness for comparison
y_log_skew  <- skew(log(y))
y_sqrt_skew <- skew(sqrt(y))
cat(sprintf("  skewness if log-transformed : %.4f\n", y_log_skew))
cat(sprintf("  skewness if sqrt-transformed: %.4f\n", y_sqrt_skew))

# histogram
pdf(OUT_HIST, width = 7, height = 5)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
hist(y, breaks = 30,
     main = sprintf("y (identity), skew=%.2f", skew(y)),
     xlab = "asthma per capita")
hist(log(y), breaks = 30,
     main = sprintf("log(y), skew=%.2f", y_log_skew),
     xlab = "log(asthma per capita)")
par(op); dev.off()
cat(sprintf("\nSaved histogram: %s\n", OUT_HIST))

# transform recommendation
abs_id   <- abs(skew(y))
abs_log  <- abs(y_log_skew)
abs_sqrt <- abs(y_sqrt_skew)
rec <- if (abs_id < 0.5) "identity (near-Gaussian)" else
       if (abs_log  < abs_id && abs_log  < 0.5) "log" else
       if (abs_sqrt < abs_id && abs_sqrt < 0.5) "sqrt" else
       if (abs_log  < abs_id) "log (still some skew, but best of the three)" else
       "identity (no transform clearly improves Gaussianity)"
cat(sprintf("\nResponse-transform recommendation: %s\n", rec))

# truth_params sanity
cat("\nsim$truth_params (real data; only nu is meaningful):\n")
print(sim$truth_params)
stopifnot(identical(sim$truth_params$nu, 1.0))
stopifnot(is.na(sim$truth_params$rho),
          is.na(sim$truth_params$sigma2),
          is.na(sim$truth_params$tau2),
          is.na(sim$truth_params$tau2_s),
          is.na(sim$truth_params$c_int))
cat("  nu == 1.0 OK; rho/sigma2/tau2/tau2_s/c_int are NA OK.\n")

# settings
cat("\nsim$settings (full):\n")
print(sim$settings)

cat("\n==============================================================\n")
cat(" Done. No model fitted.\n")
cat("==============================================================\n")
