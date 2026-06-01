## ============================================================
## scripts/pilot_p1_rho_prior_check.R
##
## ρ-prior sensitivity check. Re-runs the P1 additive fit with a
## broader log-normal prior on ρ:
##   log_rho_mu = log(0.05),  log_rho_sd = 2.0
##   ==>  ~95% prior central mass on ρ in [0.001, 2.7]
## Compared to original
##   log_rho_mu = log(0.20),  log_rho_sd = 1.0
##   ==>  ~95% prior central mass on ρ in [0.028, 1.42]
##
## Original P1 additive posterior 95% CI for ρ was [0.0145, 0.0194],
## entirely below the original prior central mass. This run asks
## whether ρ moves under a prior that puts plenty of mass at small ρ.
## ============================================================

SIM_RDS <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"
ORIG_RDS <- "output/fit_medsat_p1_additive_n500_seed1_local.rds"

stopifnot(file.exists(SIM_RDS))
sim <- readRDS(SIM_RDS)
stopifnot(sim$n == 500L, sim$p == 3L, sim$truth_params$nu == 1.0)
cat(sprintf("[sim] n=%d p=%d nu=%g  (v2 unchanged)\n",
            sim$n, sim$p, sim$truth_params$nu))

# Compile C++ kernels.
if (requireNamespace("Rcpp", quietly = TRUE)) {
  tryCatch(
    Rcpp::sourceCpp("R/ls_interaction_core.cpp", env = globalenv()),
    error = function(e) message("Rcpp::sourceCpp failed: ",
                                conditionMessage(e))
  )
}

env <- new.env(parent = globalenv())
sys.source("R/fit_ours_interaction_wrapper.R", chdir = TRUE, envir = env)

settings <- list(
  M             = 20L,
  nu            = 1.0,
  n_iter        = 3000L,
  n_burn        = 1000L,
  n_thin        = 1L,
  n_draws       = 2000L,
  orthogonalize    = TRUE,
  fit_interactions = FALSE,
  log_rho_mu       = log(0.05),
  log_rho_sd       = 2.0,
  verbose          = TRUE
)
cat("[settings]\n"); str(settings, give.attr = FALSE)
cat(sprintf("\nrho prior: log_rho_mu=log(0.05)=%.4f, log_rho_sd=%.2f\n",
            log(0.05), 2.0))
cat(sprintf("approx 95%% prior interval on rho: [%.4f, %.4f]\n",
            exp(log(0.05) - 1.96 * 2.0),
            exp(log(0.05) + 1.96 * 2.0)))

old_wd <- getwd()
setwd("R")
on.exit(setwd(old_wd), add = TRUE)
t_start <- Sys.time()
fit <- env$fit_ours_interaction(sim, settings = settings)
t_end <- Sys.time()
setwd(old_wd)
wall_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))

fit$reproducibility <- list(
  R_version    = R.version.string,
  hostname     = Sys.info()[["nodename"]],
  wd           = getwd(),
  script_path  = "scripts/pilot_p1_rho_prior_check.R",
  sim_file     = basename(SIM_RDS),
  fit_date     = Sys.time(),
  wall_sec     = wall_sec
)

host_tag <- if (Sys.info()[["nodename"]] == "hellbender") "hb" else "local"
OUT_FIT <- sprintf("output/fit_medsat_p1_additive_n500_seed1_diffuse_rho_%s.rds",
                   host_tag)
saveRDS(fit, OUT_FIT)
cat(sprintf("\nSaved fit: %s (%.1f MB)\n",
            OUT_FIT, file.info(OUT_FIT)$size / 1024^2))
cat(sprintf("Wall time: %.1f s   (%.2f min)\n", wall_sec, wall_sec / 60))

# Aliases (wrapper has the naming swap)
v <- fit$var_comp
sigma2_draws <- v$sigma2
tau2_draws   <- v$tau2_s
rho_draws    <- v$rho
tau2_s_draws <- v$tau2

ess_one <- function(x) {
  if (requireNamespace("coda", quietly = TRUE))
    as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
  else NA_real_
}
post_row <- function(x, name) {
  q <- quantile(x, c(0.025, 0.975))
  data.frame(name = name,
             mean = mean(x), q025 = q[1], q975 = q[2],
             ess = ess_one(x),
             row.names = NULL, stringsAsFactors = FALSE)
}

diffuse_tbl <- rbind(
  post_row(sigma2_draws,    "sigma2"),
  post_row(tau2_draws,      "tau2"),
  post_row(rho_draws,       "rho"),
  post_row(tau2_s_draws[[1]], "tau2_s_1"),
  post_row(tau2_s_draws[[2]], "tau2_s_2"),
  post_row(tau2_s_draws[[3]], "tau2_s_3")
)

cat("\n=========== DIFFUSE-PRIOR POSTERIOR SUMMARIES ===========\n")
print(diffuse_tbl, row.names = FALSE, digits = 4)
cat(sprintf("\nMH acceptance rates: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            fit$convergence$mh_accept[["sigma2"]],
            fit$convergence$mh_accept[["tau2"]],
            fit$convergence$mh_accept[["rho"]]))

# ----------------------------------------------------------------
# Side-by-side comparison
# ----------------------------------------------------------------
stopifnot(file.exists(ORIG_RDS))
fit_orig <- readRDS(ORIG_RDS)
vo <- fit_orig$var_comp

orig_tbl <- rbind(
  post_row(vo$sigma2,      "sigma2"),
  post_row(vo$tau2_s,      "tau2"),
  post_row(vo$rho,         "rho"),
  post_row(vo$tau2[[1]],   "tau2_s_1"),
  post_row(vo$tau2[[2]],   "tau2_s_2"),
  post_row(vo$tau2[[3]],   "tau2_s_3")
)

cat("\n=========== SIDE-BY-SIDE: original vs diffuse prior ===========\n")
sxs <- data.frame(
  param          = orig_tbl$name,
  orig_mean      = round(orig_tbl$mean,    4),
  orig_q025      = round(orig_tbl$q025,    4),
  orig_q975      = round(orig_tbl$q975,    4),
  diffuse_mean   = round(diffuse_tbl$mean, 4),
  diffuse_q025   = round(diffuse_tbl$q025, 4),
  diffuse_q975   = round(diffuse_tbl$q975, 4),
  pct_mean_change = round(100 * (diffuse_tbl$mean - orig_tbl$mean) / orig_tbl$mean, 1),
  stringsAsFactors = FALSE
)
print(sxs, row.names = FALSE)

cat(sprintf("\nESS for rho:  original = %.0f   diffuse = %.0f\n",
            orig_tbl[orig_tbl$name == "rho", "ess"],
            diffuse_tbl[diffuse_tbl$name == "rho", "ess"]))

# Verdict
rho_o_lo <- orig_tbl[orig_tbl$name == "rho", "q025"]
rho_o_hi <- orig_tbl[orig_tbl$name == "rho", "q975"]
rho_d_lo <- diffuse_tbl[diffuse_tbl$name == "rho", "q025"]
rho_d_hi <- diffuse_tbl[diffuse_tbl$name == "rho", "q975"]
rho_o_mean <- orig_tbl[orig_tbl$name == "rho", "mean"]
rho_d_mean <- diffuse_tbl[diffuse_tbl$name == "rho", "mean"]
mean_change_pct <- 100 * (rho_d_mean - rho_o_mean) / rho_o_mean

# CI overlap: Jaccard-style
overlap_lo <- max(rho_o_lo, rho_d_lo)
overlap_hi <- min(rho_o_hi, rho_d_hi)
overlap_w  <- max(0, overlap_hi - overlap_lo)
union_lo   <- min(rho_o_lo, rho_d_lo)
union_hi   <- max(rho_o_hi, rho_d_hi)
union_w    <- union_hi - union_lo
overlap_frac <- if (union_w > 0) overlap_w / union_w else 1

cat("\n=========== VERDICT ===========\n")
cat(sprintf("rho mean : original = %.4f   diffuse = %.4f   shift = %+.1f%%\n",
            rho_o_mean, rho_d_mean, mean_change_pct))
cat(sprintf("rho 95%%  : original = [%.4f, %.4f]   diffuse = [%.4f, %.4f]\n",
            rho_o_lo, rho_o_hi, rho_d_lo, rho_d_hi))
cat(sprintf("CI overlap fraction (Jaccard of intervals): %.3f\n", overlap_frac))

verdict <- if (abs(mean_change_pct) <= 20 && overlap_frac >= 0.4) {
  "DATA-IDENTIFIED. rho posterior is stable under prior changes."
} else {
  "PRIOR-PINNED. rho posterior shifted substantially under prior changes."
}
cat("\n", verdict, "\n", sep = "")

cat("\nDone.\n")
