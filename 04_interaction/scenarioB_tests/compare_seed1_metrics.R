# ============================================================
# compare_seed1_metrics.R
#
# Compute the headline metrics (RMSE + 95% pointwise coverage) for the two
# saved Scenario B seed-1 fits and print a side-by-side comparison.
# ============================================================

source("simulate_scenario_B.R")

sim    <- simulate_scenario_B(seed = 1L)
fit_y  <- readRDS("fit_scenarioB_seed1_ortho.rds")
fit_n  <- readRDS("fit_scenarioB_seed1_noortho.rds")

# ---------- helpers ----------
rmse_curve <- function(f_draws_matrix, truth) {
  # f_draws_matrix: rows = grid points, cols = draws
  pm <- rowMeans(f_draws_matrix)
  pm <- pm - mean(pm)
  tt <- truth - mean(truth)
  sqrt(mean((pm - tt)^2))
}

cov_curve <- function(f_draws_matrix, truth, level = 0.95) {
  alpha <- (1 - level) / 2
  lo <- apply(f_draws_matrix, 1, quantile, probs = alpha)
  hi <- apply(f_draws_matrix, 1, quantile, probs = 1 - alpha)
  lo_c <- lo - rowMeans(f_draws_matrix) + (rowMeans(f_draws_matrix) - mean(rowMeans(f_draws_matrix)))
  hi_c <- hi - rowMeans(f_draws_matrix) + (rowMeans(f_draws_matrix) - mean(rowMeans(f_draws_matrix)))
  # simpler: center each draw, compare to centered truth
  pm_c <- sweep(f_draws_matrix, 2, colMeans(f_draws_matrix), FUN = "-")
  lo2 <- apply(pm_c, 1, quantile, probs = alpha)
  hi2 <- apply(pm_c, 1, quantile, probs = 1 - alpha)
  tt <- truth - mean(truth)
  mean(tt >= lo2 & tt <= hi2)
}

rmse_spatial <- function(s_draws_matrix, truth_b) {
  pm <- rowMeans(s_draws_matrix)
  sqrt(mean((pm - truth_b)^2))
}

cor_spatial <- function(s_draws_matrix, truth_b) {
  pm <- rowMeans(s_draws_matrix)
  cor(pm, truth_b)
}

# ---------- compute ----------
metrics_one <- function(fit, label) {
  out <- list(label = label)
  for (j in seq_along(fit$f_main)) {
    out[[paste0("RMSE_f", j)]] <- rmse_curve(fit$f_main[[j]], sim$truth_f_grid[[j]])
    out[[paste0("cov_f",  j)]] <- cov_curve (fit$f_main[[j]], sim$truth_f_grid[[j]])
  }
  for (k in names(fit$f_int)) {
    out[[paste0("RMSE_f", k)]] <- rmse_curve(fit$f_int[[k]], sim$truth_f_int[[k]])
    out[[paste0("cov_f",  k)]] <- cov_curve (fit$f_int[[k]], sim$truth_f_int[[k]])
  }
  out$RMSE_s   <- rmse_spatial(fit$s_obs, sim$truth_s_obs)
  out$cor_s    <- cor_spatial(fit$s_obs, sim$truth_s_obs)
  out$sigma2_postmean <- mean(fit$var_comp$sigma2)
  out$tau2_s_postmean <- mean(fit$var_comp$tau2_s)
  out$rho_postmean    <- mean(fit$var_comp$rho)
  out$time_min        <- fit$timing$total_sec / 60
  out
}

m_y <- metrics_one(fit_y, "ortho=TRUE")
m_n <- metrics_one(fit_n, "ortho=FALSE")

# ---------- print ----------
keys <- setdiff(names(m_y), "label")

cat("\n=================================================================\n")
cat("  Scenario B, seed=1: paired comparison\n")
cat("=================================================================\n")
cat(sprintf("  %-20s  %12s  %12s  %12s\n",
            "metric", "ortho=FALSE", "ortho=TRUE", "ratio (Y/N)"))
cat(paste0(rep("-", 65), collapse = ""), "\n")
for (k in keys) {
  vy <- m_y[[k]]; vn <- m_n[[k]]
  if (is.null(vn) || is.null(vy)) next
  ratio <- if (abs(vn) > 1e-12) vy / vn else NA_real_
  cat(sprintf("  %-20s  %12.4f  %12.4f  %12.3f\n", k, vn, vy, ratio))
}

cat("\n  mgcv reference (Section 5.9, seed=1):\n")
cat("    RMSE_f1=0.076  cov_f1=1.00   RMSE_f2=0.079  cov_f2=0.66\n")
cat("    RMSE_f3=0.018  RMSE_f1_2=0.248  cov_f1_2=0.86\n")
cat("    RMSE_f1_3=0.097  RMSE_f2_3=0.149  RMSE_s=0.671\n")
