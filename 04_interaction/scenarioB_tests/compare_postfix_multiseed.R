# ============================================================
# compare_postfix_multiseed.R
#
# Aggregate the multi-seed Scenario B postfix results.
# Produces:
#   - per-seed RMSE and coverage table
#   - mean +/- SE summary across seeds
#   - CSVs for the dissertation table
#   - quick mgcv comparison if mgcv multi-seed RDS available
# ============================================================

source("simulate_scenario_B.R")

in_dir  <- file.path("comparison", "output", "scenarioB_postfix_multiseed")
out_dir <- file.path("comparison", "output", "scenarioB_postfix_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----
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
rmse_spatial <- function(s_draws, truth_b) {
  sqrt(mean((rowMeans(s_draws) - truth_b)^2))
}

# ---- loop over seeds ----
seeds <- 1:10
rows <- list()

for (sd in seeds) {
  rds <- file.path(in_dir, sprintf("fit_postfix_seed%d.rds", sd))
  if (!file.exists(rds)) {
    cat(sprintf("Seed %d: missing (%s); skipping\n", sd, basename(rds)))
    next
  }

  fit <- readRDS(rds)
  sim <- simulate_scenario_B(seed = sd)

  row <- list(seed = sd)
  for (j in 1:3) {
    row[[sprintf("RMSE_f%d", j)]] <- rmse_curve(fit$f_main[[j]], sim$truth_f_grid[[j]])
    row[[sprintf("cov_f%d",  j)]] <- cov_curve (fit$f_main[[j]], sim$truth_f_grid[[j]])
  }
  for (key in names(fit$f_int)) {
    row[[sprintf("RMSE_f%s", key)]] <- rmse_curve(fit$f_int[[key]], sim$truth_f_int[[key]])
    row[[sprintf("cov_f%s",  key)]] <- cov_curve (fit$f_int[[key]], sim$truth_f_int[[key]])
  }
  row$RMSE_s <- rmse_spatial(fit$s_obs, sim$truth_s_obs)

  # variance components (note: wrapper still has labels swapped)
  # gs$tau2_samples is stored as var_comp$sigma2 (so this is nugget tau2)
  # gs$sigma2_samples is stored as var_comp$tau2_s  (so this is spatial sigma2)
  row$tau2_nugget_pm  <- mean(fit$var_comp$sigma2)
  row$sigma2_spatial_pm <- mean(fit$var_comp$tau2_s)
  row$rho_pm <- mean(fit$var_comp$rho)

  rows[[length(rows) + 1]] <- row
}

if (length(rows) == 0) stop("No seed RDS files found in ", in_dir)

# ---- per-seed table ----
df <- do.call(rbind, lapply(rows, as.data.frame))
cat("\n================================================================\n")
cat(" Per-seed metrics (Scenario B postfix)\n")
cat("================================================================\n\n")
print(round(df, 4))

write.csv(df, file.path(out_dir, "scenarioB_postfix_per_seed.csv"),
          row.names = FALSE)

# ---- summary across seeds ----
metric_cols <- setdiff(names(df), "seed")
n_seeds <- nrow(df)

summary_df <- data.frame(
  metric = metric_cols,
  mean   = sapply(metric_cols, function(k) mean(df[[k]], na.rm = TRUE)),
  sd     = sapply(metric_cols, function(k) sd  (df[[k]], na.rm = TRUE)),
  se     = sapply(metric_cols, function(k) sd  (df[[k]], na.rm = TRUE) / sqrt(n_seeds)),
  median = sapply(metric_cols, function(k) median(df[[k]], na.rm = TRUE)),
  min    = sapply(metric_cols, function(k) min   (df[[k]], na.rm = TRUE)),
  max    = sapply(metric_cols, function(k) max   (df[[k]], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

cat("\n================================================================\n")
cat(sprintf(" Summary across %d seeds\n", n_seeds))
cat("================================================================\n\n")
print(round(summary_df, 4), row.names = FALSE)

write.csv(summary_df, file.path(out_dir, "scenarioB_postfix_summary.csv"),
          row.names = FALSE)

# ---- headline ----
cat("\n================================================================\n")
cat(" HEADLINE\n")
cat("================================================================\n\n")
m <- summary_df$mean
s <- summary_df$se
nm <- summary_df$metric
hd <- function(metric) {
  i <- which(nm == metric)
  if (length(i) == 0) return(NA)
  sprintf("%.4f (SE %.4f)", m[i], s[i])
}
cat(sprintf("Main effect RMSE:\n"))
cat(sprintf("  RMSE_f1 = %s\n", hd("RMSE_f1")))
cat(sprintf("  RMSE_f2 = %s\n", hd("RMSE_f2")))
cat(sprintf("  RMSE_f3 = %s\n", hd("RMSE_f3")))
cat(sprintf("\nMain effect 95%% pointwise coverage:\n"))
cat(sprintf("  cov_f1 = %s\n", hd("cov_f1")))
cat(sprintf("  cov_f2 = %s\n", hd("cov_f2")))
cat(sprintf("  cov_f3 = %s\n", hd("cov_f3")))
cat(sprintf("\nInteraction surface RMSE:\n"))
cat(sprintf("  RMSE_f1_2 = %s\n", hd("RMSE_f1_2")))
cat(sprintf("  RMSE_f1_3 = %s\n", hd("RMSE_f1_3")))
cat(sprintf("  RMSE_f2_3 = %s\n", hd("RMSE_f2_3")))
cat(sprintf("\nSpatial RE RMSE:\n"))
cat(sprintf("  RMSE_s = %s\n", hd("RMSE_s")))

cat(sprintf("\nVariance component posterior means (truth: sigma2_spatial=1, tau2_nugget=1, rho=0.06):\n"))
cat(sprintf("  sigma2_spatial = %s\n", hd("sigma2_spatial_pm")))
cat(sprintf("  tau2_nugget    = %s\n", hd("tau2_nugget_pm")))
cat(sprintf("  rho            = %s\n", hd("rho_pm")))

cat(sprintf("\nReference (mgcv on seed 1):  RMSE_f1 = 0.057\n"))

cat("\nSaved:\n")
cat(sprintf("  %s\n", file.path(out_dir, "scenarioB_postfix_per_seed.csv")))
cat(sprintf("  %s\n", file.path(out_dir, "scenarioB_postfix_summary.csv")))
cat("\nDone.\n")
