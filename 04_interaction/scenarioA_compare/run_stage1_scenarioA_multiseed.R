# ============================================================
# run_stage1_scenarioA_multiseed.R
# 
# Runs Scenario A across 10 seeds for all three methods (ours,
# mgcv, inla). Saves each fit individually so a crash does not
# lose finished work; resumes by skipping any (method, seed)
# combo whose RDS already exists.
# 
# Outputs:
#   ../output/fits/fit_<method>_scenarioA_seed<N>.rds   per fit
#   ../output/metrics/scenarioA_per_seed.csv            tidy metrics
#   ../output/metrics/scenarioA_summary.csv             aggregate
# 
# Expected runtime, n_seeds = 10:
#   ours:  ~6 hours total (35 min/seed)
#   inla:  ~3 minutes total
#   mgcv:  ~15 seconds total
# 
# To run only a subset of seeds or methods, edit the SEEDS and
# METHODS variables at the top.  Resume just by re-running.
# 
# To run unattended (e.g. overnight): in R/RStudio,
#   sink("stage1_log.txt", split = TRUE)
#   source("run_stage1_scenarioA_multiseed.R")
#   sink()
# ============================================================

# ---- configuration ----
SEEDS    <- 1:10
METHODS  <- c("ours", "mgcv", "inla")
SCENARIO <- "A"

OUT_FITS    <- file.path("..", "output", "fits")
OUT_METRICS <- file.path("..", "output", "metrics")
dir.create(OUT_FITS,    recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_METRICS, recursive = TRUE, showWarnings = FALSE)

# ---- sources ----
source("simulate_scenario_A.R")
source("compute_recovery_metrics.R")
# wrappers are sourced lazily inside the loop so a missing optional
# package (e.g. INLA) only matters if you actually selected that method.

# ---- helpers ----
fit_path <- function(method, seed) {
  file.path(OUT_FITS, sprintf("fit_%s_scenario%s_seed%d.rds",
                              method, SCENARIO, seed))
}

# Clean status print
banner <- function(msg) {
  cat("\n", strrep("-", 60), "\n", msg, "\n",
      strrep("-", 60), "\n", sep = "")
}

# ---- main loop: methods nested inside seeds ----
# Order: for each seed, do mgcv first (fast), then inla, then ours
# (slow). This way we get partial results across seeds quickly.
order_within_seed <- intersect(c("mgcv", "inla", "ours"), METHODS)

t_global <- proc.time()

# Track per-seed timings for ETA
seed_times <- numeric(0)

for (s in seq_along(SEEDS)) {
  seed <- SEEDS[s]
  banner(sprintf("seed %d / %d", seed, length(SEEDS)))
  t_seed_start <- proc.time()
  
  # simulate ONCE per seed -> shared across methods (paired comparison)
  sim <- simulate_scenario_A(seed = seed)
  
  for (m in order_within_seed) {
    fp <- fit_path(m, seed)
    if (file.exists(fp)) {
      cat(sprintf("[skip]  %s  (already done)\n", basename(fp)))
      next
    }
    
    # lazy-load wrapper
    if (m == "ours") source("fit_ours_wrapper.R")
    if (m == "mgcv") source("fit_mgcv_wrapper.R")
    if (m == "inla") source("fit_inla_wrapper.R")
    
    cat(sprintf("[fit]   %s  seed=%d ... ", m, seed))
    flush.console()
    t_m <- proc.time()
    
    fit <- tryCatch({
      switch(m,
        ours = fit_ours(sim, settings = list(verbose = FALSE)),
        mgcv = fit_mgcv(sim),
        inla = fit_inla(sim)
      )
    }, error = function(e) {
      cat(sprintf("FAILED: %s\n", conditionMessage(e)))
      NULL
    })
    
    if (is.null(fit)) next
    
    saveRDS(fit, fp)
    cat(sprintf("done in %.1fs\n", (proc.time() - t_m)[3]))
  }
  
  seed_sec <- (proc.time() - t_seed_start)[3]
  seed_times <- c(seed_times, seed_sec)
  
  # ETA based on average seed time so far
  if (s < length(SEEDS)) {
    avg_seed <- mean(seed_times)
    eta_min  <- (length(SEEDS) - s) * avg_seed / 60
    cat(sprintf("[eta]   %.0f seconds for this seed; ~%.1f min remaining\n",
                seed_sec, eta_min))
  }
}

banner("aggregating metrics")

# ---- collect per-seed metrics across all completed fits ----
all_rows <- list()
for (seed in SEEDS) {
  sim <- simulate_scenario_A(seed = seed)   # same DGP, regen for truth
  for (m in METHODS) {
    fp <- fit_path(m, seed)
    if (!file.exists(fp)) {
      cat(sprintf("[miss]  %s seed=%d  (no fit RDS)\n", m, seed))
      next
    }
    fit <- readRDS(fp)
    row <- compute_recovery_metrics(fit, sim)
    all_rows[[paste(m, seed, sep = "_")]] <- row
  }
}

per_seed <- do.call(rbind, all_rows)
rownames(per_seed) <- NULL

per_seed_path <- file.path(OUT_METRICS, sprintf("scenario%s_per_seed.csv", SCENARIO))
write.csv(per_seed, per_seed_path, row.names = FALSE)
cat(sprintf("\n[save]  %s  (%d rows)\n", per_seed_path, nrow(per_seed)))

# ---- summarize ----
if (nrow(per_seed) > 0L) {
  summ <- summarize_metrics(per_seed)
  summ_path <- file.path(OUT_METRICS, sprintf("scenario%s_summary.csv", SCENARIO))
  write.csv(summ, summ_path, row.names = FALSE)
  cat(sprintf("[save]  %s\n", summ_path))
  
  banner("compact summary (mean +/- MCSE across seeds)")
  
  # print a compact human-readable table
  fmt <- function(m, se) {
    if (is.na(m)) return("    NA   ")
    if (is.na(se) || se == 0) return(sprintf("%6.3f      ", m))
    sprintf("%6.3f \u00b1 %.3f", m, se)
  }
  for (m in METHODS) {
    rows <- summ[summ$method == m, , drop = FALSE]
    if (nrow(rows) == 0L) next
    r <- rows[1, ]
    cat(sprintf("\n[%s]  n_seeds = %d\n", m, r$n_seeds))
    cat(sprintf("  rmse_f1     : %s\n",
                fmt(r$mean_rmse_f1, r$mcse_rmse_f1)))
    cat(sprintf("  rmse_f2     : %s\n",
                fmt(r$mean_rmse_f2, r$mcse_rmse_f2)))
    cat(sprintf("  rmse_f3     : %s\n",
                fmt(r$mean_rmse_f3, r$mcse_rmse_f3)))
    cat(sprintf("  rmse_f_mean : %s\n",
                fmt(r$mean_rmse_f_mean, r$mcse_rmse_f_mean)))
    cat(sprintf("  cov95_f1    : %s\n",
                fmt(r$mean_cov95_f1, r$mcse_cov95_f1)))
    cat(sprintf("  cov95_f2    : %s\n",
                fmt(r$mean_cov95_f2, r$mcse_cov95_f2)))
    cat(sprintf("  cov95_f3    : %s\n",
                fmt(r$mean_cov95_f3, r$mcse_cov95_f3)))
    cat(sprintf("  cov95_f_mean: %s\n",
                fmt(r$mean_cov95_f_mean, r$mcse_cov95_f_mean)))
    cat(sprintf("  rmse_s      : %s\n",
                fmt(r$mean_rmse_s, r$mcse_rmse_s)))
    cat(sprintf("  cor_s       : %s\n",
                fmt(r$mean_cor_s,  r$mcse_cor_s)))
    cat(sprintf("  total_sec   : %s\n",
                fmt(r$mean_total_sec, r$mcse_total_sec)))
  }
}

cat(sprintf("\n[total wall time: %.1f minutes]\n",
            (proc.time() - t_global)[3] / 60))
cat("\nOK.\n")
