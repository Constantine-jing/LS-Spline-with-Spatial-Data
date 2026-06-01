## ============================================================
## scripts/run_n1500_seed3.R
##
## JOB 3 of the four-job Hellbender batch — multi-seed robustness
## for the null-interaction finding. Identical chain length to the
## original n=1500 run (n_iter=3000, n_draws=2000) but on a fresh
## stratified subsample (seed=3).
##
## Settings are identical to run_n1500_seed2.R / run_n1500_seed4.R
## except for the sim file. See run_n1500_headline_seed1.R for the
## polished single-seed companion.
##
## Output: output/fit_medsat_n1500_seed3_allpairs_hb.rds
##
## Usage on Hellbender (submitted by .sub):
##   Rscript scripts/run_n1500_seed3.R
##
## Local dry-run (tiny, discards output):
##   MEDSAT_DRY_RUN=1 Rscript scripts/run_n1500_seed3.R
## ============================================================

IS_DRY_RUN <- identical(Sys.getenv("MEDSAT_DRY_RUN"), "1")
SEED_TAG   <- 3L

cat("============================================================\n")
cat(sprintf(" run_n1500_seed%d.R  |  %s\n", SEED_TAG, Sys.time()))
cat(sprintf(" Hostname:   %s\n", Sys.info()[["nodename"]]))
cat(sprintf(" R version:  %s\n", R.version.string))
cat(sprintf(" Mode:       %s\n", if (IS_DRY_RUN) "DRY-RUN (tiny settings)" else "PRODUCTION"))
cat("============================================================\n\n")

# ---- Hostname guard (production only) ----
host <- Sys.info()[["nodename"]]
if (!IS_DRY_RUN) {
  if (!grepl("hellbender", host, ignore.case = TRUE)) {
    stop(paste0(
      "Hostname guard: expected a Hellbender node, got '", host, "'.\n",
      "  This script is intended for the Hellbender HPC cluster.\n",
      "  For a local dry-run, set env var MEDSAT_DRY_RUN=1."
    ))
  }
}

# ---- Settings ----
if (IS_DRY_RUN) {
  settings <- list(
    M                = 10L,
    nu               = 1.0,
    n_iter           = 100L,
    n_burn           = 50L,
    n_thin           = 1L,
    n_draws          = 50L,
    fit_interactions = TRUE,
    orthogonalize    = TRUE,
    log_rho_mu       = log(0.2),
    log_rho_sd       = 1.0,
    verbose          = TRUE
  )
  OUT_FIT <- file.path(tempdir(),
                       sprintf("fit_medsat_n1500_seed%d_allpairs_DRYRUN.rds", SEED_TAG))
} else {
  settings <- list(
    M                = 20L,
    nu               = 1.0,
    n_iter           = 3000L,
    n_burn           = 1000L,
    n_thin           = 1L,
    n_draws          = 2000L,
    fit_interactions = TRUE,
    orthogonalize    = TRUE,
    log_rho_mu       = log(0.2),
    log_rho_sd       = 1.0,
    mh_sd_log_rho    = 0.15,
    mh_sd_log_tau2   = 0.6,
    mh_sd_log_sigma2 = 0.3,
    verbose          = TRUE
  )
  OUT_FIT <- sprintf("output/fit_medsat_n1500_seed%d_allpairs_hb.rds", SEED_TAG)
}

expected_draws <- floor((settings$n_iter - settings$n_burn) / settings$n_thin)
stopifnot(expected_draws == settings$n_draws)

cat("[settings]\n"); str(settings, give.attr = FALSE)

SIM_RDS <- sprintf("output/sim_medsat_london_asthma_n1500_seed%d_v2.rds", SEED_TAG)

# ---- Guard: don't overwrite an existing production fit ----
if (!IS_DRY_RUN && file.exists(OUT_FIT)) {
  stop(sprintf(
    "Output already exists: %s\n  Rename or delete it first to avoid overwriting.",
    OUT_FIT))
}

# ---- Load sim ----
cat(sprintf("\n[sim] loading %s\n", SIM_RDS))
if (!file.exists(SIM_RDS)) stop(sprintf("Sim RDS not found: %s", SIM_RDS))
sim <- readRDS(SIM_RDS)
stopifnot(sim$n == 1500L, sim$p == 3L)
stopifnot(identical(sim$truth_params$nu, 1.0))
stopifnot(is.na(sim$truth_params$rho))
stopifnot(setequal(sim$int_keys, c("1_2", "1_3", "2_3")))
stopifnot(identical(sim$seed, SEED_TAG))
cat(sprintf("[sim] n=%d  p=%d  nu=%g  seed=%d  int_keys={%s}\n",
            sim$n, sim$p, sim$truth_params$nu, sim$seed,
            paste(sim$int_keys, collapse = ", ")))

# ---- Compile C++ kernels ----
cat("\n[Rcpp] compiling C++ kernels...\n")
t_cpp <- system.time(
  tryCatch(
    Rcpp::sourceCpp("R/ls_interaction_core.cpp", env = globalenv()),
    error = function(e) stop("Rcpp::sourceCpp failed: ", conditionMessage(e))
  )
)
stopifnot(exists("khatri_rao_cpp"))
cat(sprintf("[Rcpp] compiled in %.1f s\n", t_cpp[3]))

# ---- Source wrapper ----
cat("\n[wrapper] sourcing...\n")
env <- new.env(parent = globalenv())
sys.source("R/fit_ours_interaction_wrapper.R", chdir = TRUE, envir = env)

# ---- Fit ----
cat(sprintf("\n[fit] starting at %s\n", Sys.time()))
run_with_chdir <- function(fn, ...) {
  old <- getwd(); setwd("R"); on.exit(setwd(old), add = TRUE)
  fn(...)
}

t_wall_start <- Sys.time()
fit <- run_with_chdir(env$fit_ours_interaction, sim, settings = settings)
t_wall_end   <- Sys.time()
wall_sec <- as.numeric(difftime(t_wall_end, t_wall_start, units = "secs"))

cat(sprintf("\n[fit] completed at %s\n", t_wall_end))
cat(sprintf("[fit] wall time: %.1f s  (%.2f min  |  %.2f h)\n",
            wall_sec, wall_sec / 60, wall_sec / 3600))
cat(sprintf("[fit] Gibbs: %.1f s   post-proc: %.1f s\n",
            fit$timing$fit_sec, fit$timing$post_sec))
cat(sprintf("[fit] MH accepts: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            fit$convergence$mh_accept[["sigma2"]],
            fit$convergence$mh_accept[["tau2"]],
            fit$convergence$mh_accept[["rho"]]))

# ---- Reproducibility metadata ----
fit$reproducibility <- list(
  R_version   = R.version.string,
  hostname    = host,
  wd          = getwd(),
  script_path = sprintf("scripts/run_n1500_seed%d.R", SEED_TAG),
  sim_file    = basename(SIM_RDS),
  sim_n       = sim$n,
  sim_seed    = sim$seed,
  fit_date    = t_wall_end,
  wall_sec    = wall_sec,
  is_dry_run  = IS_DRY_RUN
)

# ---- Posterior summaries ----
v <- fit$var_comp
post_summary <- function(x, name) {
  q <- quantile(x, c(0.025, 0.5, 0.975))
  data.frame(
    name = name, mean = mean(x), sd = sd(x),
    q025 = q[1], q500 = q[2], q975 = q[3],
    ess  = if (requireNamespace("coda", quietly = TRUE)) {
             as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
           } else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
}

post_tbl <- rbind(
  post_summary(v$sigma2,             "sigma2 (noise)"),
  post_summary(v$tau2_s,             "tau2 (spatial)"),
  post_summary(v$rho,                "rho (Matern range)"),
  post_summary(v$tau2[[1]],          "tau2_s[1] (NO2)"),
  post_summary(v$tau2[[2]],          "tau2_s[2] (NDVI)"),
  post_summary(v$tau2[[3]],          "tau2_s[3] (IMD)"),
  post_summary(v$tau2_int[["1_2"]],  "tau2_int[1_2] (NO2xNDVI)"),
  post_summary(v$tau2_int[["1_3"]],  "tau2_int[1_3] (NO2xIMD)"),
  post_summary(v$tau2_int[["2_3"]],  "tau2_int[2_3] (NDVIxIMD)")
)

cat("\n==============================================================\n")
cat(" POSTERIOR SUMMARIES\n")
cat("==============================================================\n")
print(post_tbl, row.names = FALSE, digits = 4)

# ---- Save ----
if (IS_DRY_RUN) {
  cat(sprintf("\n[dry-run] output written to temp: %s\n", OUT_FIT))
  cat("[dry-run] NOT saving to production path. Discarding after session.\n")
  saveRDS(fit, OUT_FIT)
  sz_mb <- file.info(OUT_FIT)$size / 1024^2
  cat(sprintf("[dry-run] temp file size: %.1f MB\n", sz_mb))
} else {
  saveRDS(fit, OUT_FIT)
  sz_mb <- file.info(OUT_FIT)$size / 1024^2
  cat(sprintf("\nSaved: %s  (%.1f MB)\n", OUT_FIT, sz_mb))
}

cat("\n==============================================================\n")
cat(sprintf(" %s complete\n",
            if (IS_DRY_RUN) "DRY-RUN" else "PRODUCTION RUN"))
cat(sprintf(" %s\n", Sys.time()))
cat("==============================================================\n")
