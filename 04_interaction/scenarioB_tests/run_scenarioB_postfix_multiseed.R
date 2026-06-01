# ============================================================
# run_scenarioB_postfix_multiseed.R
#
# Multi-seed Scenario B with the patched wrapper.
# Parallel-friendly: each invocation runs ONE seed at a time
# (controlled by the seeds argument). Run multiple R sessions
# with different seeds in parallel for ~10x wall-clock speedup.
#
# USAGE:
#   # Run all 10 seeds sequentially in one session (~20-30 hours):
#   Rscript run_scenarioB_postfix_multiseed.R
#
#   # Or run subset of seeds in one session, rest in another:
#   Rscript run_scenarioB_postfix_multiseed.R 1 2 3
#   Rscript run_scenarioB_postfix_multiseed.R 4 5 6
#   Rscript run_scenarioB_postfix_multiseed.R 7 8 9 10
#
# Each seed produces:
#   comparison/output/scenarioB_postfix_multiseed/
#     fit_postfix_seed{N}.rds
#     log_postfix_seed{N}.txt
#
# After all seeds finish, run compare_postfix_multiseed.R for the
# summary table.
# ============================================================

# ---- parse seed argument(s) ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  seeds <- 1:10
  cat("No seeds specified; running all 10.\n")
} else {
  seeds <- as.integer(args)
  if (any(is.na(seeds))) stop("Invalid seed argument(s): must be integers")
  cat(sprintf("Running seeds: %s\n", paste(seeds, collapse = ", ")))
}

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

out_dir <- file.path("comparison", "output", "scenarioB_postfix_multiseed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- seed loop ----
for (sd in seeds) {

  out_rds <- file.path(out_dir, sprintf("fit_postfix_seed%d.rds", sd))
  log_txt <- file.path(out_dir, sprintf("log_postfix_seed%d.txt",  sd))

  if (file.exists(out_rds)) {
    cat(sprintf("\n[seed %d] already done (%s exists), skipping.\n",
                sd, basename(out_rds)))
    next
  }

  cat(sprintf("\n================================================================\n"))
  cat(sprintf(" Seed %d / Scenario B postfix\n", sd))
  cat(sprintf("================================================================\n"))

  # capture log
  con <- file(log_txt, open = "wt")
  sink(con, type = "message")
  sink(con, type = "output", split = TRUE)

  cat(sprintf("Seed %d started at %s\n", sd, format(Sys.time())))

  sim <- simulate_scenario_B(seed = sd)
  cat(sprintf("Simulated Scenario B at seed=%d: n=%d, p=%d\n",
              sd, sim$n, sim$p))

  t0 <- proc.time()
  fit_sd <- fit_ours_interaction(
    sim,
    settings = list(
      orthogonalize = TRUE,
      verbose       = TRUE
    )
  )
  elapsed_min <- (proc.time() - t0)[3] / 60
  cat(sprintf("\nFit complete in %.1f minutes.\n", elapsed_min))

  saveRDS(fit_sd, out_rds)
  cat(sprintf("Saved: %s\n", out_rds))

  # quick verdict for this seed
  rmse_curve <- function(M, truth) {
    pm <- rowMeans(M); pm <- pm - mean(pm)
    tt <- truth - mean(truth); sqrt(mean((pm - tt)^2))
  }
  cov_curve <- function(M, truth, level = 0.95) {
    a <- (1 - level) / 2
    Mc <- sweep(M, 2, colMeans(M), FUN = "-")
    lo <- apply(Mc, 1, quantile, probs = a)
    hi <- apply(Mc, 1, quantile, probs = 1 - a)
    tt <- truth - mean(truth); mean(tt >= lo & tt <= hi)
  }

  cat("\n---- Quick verdict for this seed ----\n")
  for (j in 1:3) {
    cat(sprintf("  RMSE_f%d = %.4f   cov_f%d = %.3f\n",
                j,
                rmse_curve(fit_sd$f_main[[j]], sim$truth_f_grid[[j]]),
                j,
                cov_curve(fit_sd$f_main[[j]], sim$truth_f_grid[[j]])))
  }
  for (key in names(fit_sd$f_int)) {
    cat(sprintf("  RMSE_f%-3s = %.4f\n", key,
                rmse_curve(fit_sd$f_int[[key]], sim$truth_f_int[[key]])))
  }
  cat(sprintf("Seed %d ended at %s\n", sd, format(Sys.time())))

  sink(type = "output")
  sink(type = "message")
  close(con)
}

cat("\nAll requested seeds complete.\n")
cat("Run compare_postfix_multiseed.R for the summary table.\n")
