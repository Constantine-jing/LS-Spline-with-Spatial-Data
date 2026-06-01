## ============================================================
## scripts/pilot_p1_fit.R
##
## P1 -- Bayesian ADDITIVE-only fit (NO interaction blocks sampled),
## production MCMC budget. Uses the new fit_interactions = FALSE flag.
##
## Diagnostic plots live in scripts/pilot_p1_diagnostics.R (separate
## task). This script does fit + save + a minimal text report, plus a
## side-by-side comparison against the renamed full-interaction run
## (output/fit_medsat_p1plus_allpairs_n500_seed1_<host>.rds) if it
## exists on disk.
## ============================================================

SIM_RDS <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"

# 1. Confirm v2 sim unchanged.
stopifnot(file.exists(SIM_RDS))
sim <- readRDS(SIM_RDS)
stopifnot(sim$n == 500L, sim$p == 3L)
stopifnot(sim$truth_params$nu == 1.0)
cat(sprintf("[sim] n=%d p=%d nu=%g  (v2 unchanged)\n",
            sim$n, sim$p, sim$truth_params$nu))

# 2. Compile C++ kernels into globalenv so gibbs_interaction picks them up.
if (requireNamespace("Rcpp", quietly = TRUE)) {
  tryCatch(
    Rcpp::sourceCpp("R/ls_interaction_core.cpp", env = globalenv()),
    error = function(e) message("Rcpp::sourceCpp failed; falling back to pure R: ",
                                conditionMessage(e))
  )
}

# 3. Source the wrapper.
env <- new.env(parent = globalenv())
sys.source("R/fit_ours_interaction_wrapper.R", chdir = TRUE, envir = env)

# 4. Fit. The wrapper does source("fit_spatial_reml.R") inside its body, so
# briefly chdir to R/ around the call so the relative path resolves.
settings <- list(
  M             = 20L,
  nu            = 1.0,
  n_iter        = 3000L,
  n_burn        = 1000L,
  n_thin        = 1L,
  n_draws       = 2000L,
  orthogonalize    = TRUE,
  fit_interactions = FALSE,
  verbose          = TRUE
)
cat("[settings]\n"); str(settings, give.attr = FALSE)

old_wd <- getwd()
setwd("R")
on.exit(setwd(old_wd), add = TRUE)
t_start <- Sys.time()
fit <- env$fit_ours_interaction(sim, settings = settings)
t_end <- Sys.time()
setwd(old_wd)
wall_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
per_iter <- wall_sec / settings$n_iter

# 5. Reproducibility metadata.
fit$reproducibility <- list(
  R_version    = R.version.string,
  hostname     = Sys.info()[["nodename"]],
  wd           = getwd(),
  script_path  = "scripts/pilot_p1_fit.R",
  sim_file     = basename(SIM_RDS),
  fit_date     = Sys.time(),
  wall_sec     = wall_sec
)

# 6. Save.
host_tag <- if (Sys.info()[["nodename"]] == "hellbender") "hb" else "local"
OUT_FIT <- sprintf("output/fit_medsat_p1_additive_n500_seed1_%s.rds", host_tag)
saveRDS(fit, OUT_FIT)
sz_mb <- file.info(OUT_FIT)$size / 1024^2
cat(sprintf("\nSaved fit: %s (%.1f MB)\n", OUT_FIT, sz_mb))

# ----------------------------------------------------------------
# Aliases for canonical-schema reading (wrapper has a naming swap).
# ----------------------------------------------------------------
v <- fit$var_comp
sigma2_draws  <- v$sigma2      # noise sigma^2
tau2_draws    <- v$tau2_s      # spatial GP variance tau^2
rho_draws     <- v$rho         # Matern range
tau2_s_draws  <- v$tau2        # list of per-smooth smoothing variance

post_summary <- function(x, name) {
  q <- quantile(x, c(0.025, 0.5, 0.975))
  data.frame(
    name = name,
    mean = mean(x),
    sd   = sd(x),
    q025 = q[1], q500 = q[2], q975 = q[3],
    ess  = if (requireNamespace("coda", quietly = TRUE)) {
             as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
           } else NA_real_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
post_tbl <- rbind(
  post_summary(sigma2_draws,    "sigma2 (noise)"),
  post_summary(tau2_draws,      "tau2 (spatial)"),
  post_summary(tau2_s_draws[[1]], "tau2_s[1] (NO2)"),
  post_summary(tau2_s_draws[[2]], "tau2_s[2] (NDVI)"),
  post_summary(tau2_s_draws[[3]], "tau2_s[3] (IMD)"),
  post_summary(rho_draws,       "rho (Matern range)")
)

# ----------------------------------------------------------------
# Report
# ----------------------------------------------------------------
cat("\n==============================================================\n")
cat(" P1 ADDITIVE REPORT\n")
cat("==============================================================\n")
cat(sprintf("Wall time         : %.1f s   (~ %.2f min)\n", wall_sec, wall_sec / 60))
cat(sprintf("Per-iteration time: %.4f s   (3000 iters, M=%d)\n",
            per_iter, settings$M))
cat(sprintf("Gibbs only        : %.1f s     post-proc: %.1f s\n",
            fit$timing$fit_sec, fit$timing$post_sec))
cat(sprintf("fit$f_int length  : %d   (expect 0)\n", length(fit$f_int)))
cat(sprintf("fit$settings$fit_interactions: %s   (expect FALSE)\n",
            fit$settings$fit_interactions))
cat(sprintf("fit$settings$int_keys: c(%s)\n",
            paste(shQuote(fit$settings$int_keys), collapse = ",")))
cat(sprintf("RDS file size     : %.1f MB\n", sz_mb))

cat("\nPosterior summaries:\n")
print(post_tbl, row.names = FALSE, digits = 4)

cat(sprintf("\nMH acceptance rates: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            fit$convergence$mh_accept[["sigma2"]],
            fit$convergence$mh_accept[["tau2"]],
            fit$convergence$mh_accept[["rho"]]))

# ----------------------------------------------------------------
# Side-by-side vs full-interaction (early P2-equivalent)
# ----------------------------------------------------------------
ALL_RDS <- sprintf("output/fit_medsat_p1plus_allpairs_n500_seed1_%s.rds", host_tag)
if (file.exists(ALL_RDS)) {
  cat("\n==============================================================\n")
  cat(" SIDE-BY-SIDE: additive  vs  full-interaction (renamed)\n")
  cat("==============================================================\n")
  fit_all <- readRDS(ALL_RDS)
  va <- fit_all$var_comp

  summ_one <- function(x) {
    q <- quantile(x, c(0.025, 0.975))
    c(mean = mean(x), q025 = q[1], q975 = q[2])
  }
  rows <- rbind(
    cbind(parameter = "sigma2 (noise)",
          additive  = round(summ_one(v$sigma2),  4),
          full_int  = round(summ_one(va$sigma2), 4)),
    cbind(parameter = "tau2 (spatial)",
          additive  = round(summ_one(v$tau2_s),  4),
          full_int  = round(summ_one(va$tau2_s), 4)),
    cbind(parameter = "rho (range)",
          additive  = round(summ_one(v$rho),     4),
          full_int  = round(summ_one(va$rho),    4))
  )
  print(as.data.frame(rows))

  add_walltime <- fit$reproducibility$wall_sec
  full_walltime <- if (!is.null(fit_all$reproducibility$wall_sec)) {
    fit_all$reproducibility$wall_sec
  } else {
    fit_all$timing$total_sec
  }
  cat(sprintf("\nWall-time speedup: additive=%.1fs   full=%.1fs   speedup=%.2fx\n",
              add_walltime, full_walltime, full_walltime / add_walltime))

  # One-line interpretation
  shift_sigma <- mean(v$sigma2) - mean(va$sigma2)
  shift_tau   <- mean(v$tau2_s) - mean(va$tau2_s)
  shift_rho   <- mean(v$rho)    - mean(va$rho)
  ratio_sigma <- mean(v$sigma2) / mean(va$sigma2)
  ratio_tau   <- mean(v$tau2_s) / mean(va$tau2_s)
  ratio_rho   <- mean(v$rho)    / mean(va$rho)
  cat(sprintf("\nMean shifts (additive minus full):  sigma2 %+.3f (%.2fx)  tau2 %+.3f (%.2fx)  rho %+.4f (%.2fx)\n",
              shift_sigma, ratio_sigma, shift_tau, ratio_tau,
              shift_rho,   ratio_rho))
} else {
  cat(sprintf("\nFull-interaction RDS not found at %s; skipping side-by-side.\n",
              ALL_RDS))
}

cat("\nDone.\n")
