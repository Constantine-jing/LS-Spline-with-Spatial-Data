# ============================================================
# compare_h3_metrics.R
#
# Side-by-side metrics for the four Scenario B seed-1 fits:
#   1. ortho=FALSE              (original, basis-correlation present)
#   2. ortho=TRUE               (Hypothesis 1 fix: orthogonalized W_uv)
#   3. H3 confirmation          (ortho=TRUE + tau2_s_1 clamped at 0.005)
#   4. H3 prior fix tighten b_s (ortho=TRUE + b_smooth = 0.0005)
#
# Mirrors compare_seed1_metrics.R / meeting_comparison.R conventions
# (centered RMSE, centered-draws coverage, spatial RMSE on raw b).
#
# Reads:
#   fit_scenarioB_seed1_noortho.rds
#   fit_scenarioB_seed1_ortho.rds
#   comparison/output/scenarioB_h3_confirm/fit_h3_confirm_seed1.rds
#   comparison/output/scenarioB_h3_tightbs/fit_h3_tightbs_seed1.rds
#
# Outputs:
#   comparison/output/scenarioB_h3_summary.csv
#   comparison/output/scenarioB_h3_tau2s_table.csv
# ============================================================

source("simulate_scenario_B.R")

# ---- helpers (match compare_seed1_metrics.R) ----
rmse_curve <- function(M, truth) {
  pm <- rowMeans(M); pm <- pm - mean(pm)
  tt <- truth - mean(truth)
  sqrt(mean((pm - tt)^2))
}
cov_curve <- function(M, truth, level = 0.95) {
  a <- (1 - level) / 2
  Mc <- sweep(M, 2, colMeans(M), FUN = "-")
  lo <- apply(Mc, 1, quantile, probs = a)
  hi <- apply(Mc, 1, quantile, probs = 1 - a)
  tt <- truth - mean(truth)
  mean(tt >= lo & tt <= hi)
}
rmse_spatial <- function(s_draws, truth_b) {
  sqrt(mean((rowMeans(s_draws) - truth_b)^2))
}
cor_spatial <- function(s_draws, truth_b) {
  cor(rowMeans(s_draws), truth_b)
}

# ---- load ----
sim <- simulate_scenario_B(seed = 1L)

paths <- list(
  noortho = "fit_scenarioB_seed1_noortho.rds",
  ortho   = "fit_scenarioB_seed1_ortho.rds",
  h3conf  = file.path("comparison", "output", "scenarioB_h3_confirm",
                      "fit_h3_confirm_seed1.rds"),
  h3tbs   = file.path("comparison", "output", "scenarioB_h3_tightbs",
                      "fit_h3_tightbs_seed1.rds")
)

fits <- list()
for (lab in names(paths)) {
  if (file.exists(paths[[lab]])) {
    fits[[lab]] <- readRDS(paths[[lab]])
  } else {
    cat(sprintf("MISSING: %-8s -> %s\n", lab, paths[[lab]]))
  }
}
present <- names(fits)
if (length(present) < 2) {
  stop("Need at least two fits to compare. Run the experiment scripts first.")
}

# ---- compute metrics for each fit ----
metrics_one <- function(fit) {
  out <- list()
  for (j in seq_along(fit$f_main)) {
    out[[paste0("RMSE_f", j)]] <- rmse_curve(fit$f_main[[j]], sim$truth_f_grid[[j]])
    out[[paste0("cov_f",  j)]] <- cov_curve (fit$f_main[[j]], sim$truth_f_grid[[j]])
  }
  for (k in names(fit$f_int)) {
    out[[paste0("RMSE_f", k)]] <- rmse_curve(fit$f_int[[k]], sim$truth_f_int[[k]])
    out[[paste0("cov_f",  k)]] <- cov_curve (fit$f_int[[k]], sim$truth_f_int[[k]])
  }
  out$RMSE_s          <- rmse_spatial(fit$s_obs, sim$truth_s_obs)
  out$cor_s           <- cor_spatial (fit$s_obs, sim$truth_s_obs)
  out$sigma2_postmean <- mean(fit$var_comp$sigma2)
  out$tau2_s_postmean <- mean(fit$var_comp$tau2_s)
  out$rho_postmean    <- mean(fit$var_comp$rho)
  for (j in seq_along(fit$var_comp$tau2)) {
    out[[paste0("tau2_s_main_", j)]] <- mean(fit$var_comp$tau2[[j]])
  }
  for (k in names(fit$var_comp$tau2_int)) {
    out[[paste0("tau2_s_int_", k)]] <- mean(fit$var_comp$tau2_int[[k]])
  }
  out$time_min <- fit$timing$total_sec / 60
  out
}

mets <- lapply(fits, metrics_one)

# ---- print headline table ----
keys <- unique(unlist(lapply(mets, names)))
col_labels <- c(
  noortho = "ortho=FALSE",
  ortho   = "ortho=TRUE",
  h3conf  = "H3 clamp t1",
  h3tbs   = "tight b_s"
)

cat("\n================================================================\n")
cat(" Scenario B, seed = 1: four-way comparison\n")
cat(" ortho=FALSE | ortho=TRUE | H3 clamp tau2_s_1=0.005 | tight b_s=5e-4\n")
cat("================================================================\n\n")

hdr <- sprintf("  %-22s", "metric")
for (lab in present) hdr <- paste0(hdr, sprintf("  %12s", col_labels[lab]))
cat(hdr, "\n")
cat(paste0(rep("-", nchar(hdr)), collapse = ""), "\n")

# Order: main effects first, then interactions, then spatial, then var comps,
# then per-smooth tau2_s, then timing.
section <- function(label) {
  cat(sprintf("\n  -- %s --\n", label))
}

section("Main effect RMSE")
for (j in 1:3) {
  k <- sprintf("RMSE_f%d", j)
  row <- sprintf("  %-22s", k)
  for (lab in present) row <- paste0(row, sprintf("  %12.4f", mets[[lab]][[k]]))
  cat(row, "\n")
}

section("Main effect 95% pointwise coverage")
for (j in 1:3) {
  k <- sprintf("cov_f%d", j)
  row <- sprintf("  %-22s", k)
  for (lab in present) row <- paste0(row, sprintf("  %12.4f", mets[[lab]][[k]]))
  cat(row, "\n")
}

section("Interaction surface RMSE / coverage")
int_keys <- intersect(c("1_2", "1_3", "2_3"),
                      gsub("^RMSE_f", "", grep("^RMSE_f[0-9]_", keys, value = TRUE)))
for (k in int_keys) {
  for (pref in c("RMSE_f", "cov_f")) {
    metric <- paste0(pref, k)
    row <- sprintf("  %-22s", metric)
    for (lab in present) {
      v <- mets[[lab]][[metric]]
      row <- paste0(row, sprintf("  %12s",
                                 if (is.null(v)) "NA" else sprintf("%12.4f", v)))
    }
    cat(row, "\n")
  }
}

section("Spatial random effect")
for (k in c("RMSE_s", "cor_s")) {
  row <- sprintf("  %-22s", k)
  for (lab in present) row <- paste0(row, sprintf("  %12.4f", mets[[lab]][[k]]))
  cat(row, "\n")
}

section("Variance components (truth: sigma2=1, tau2_s=1, rho=0.06)")
for (k in c("sigma2_postmean", "tau2_s_postmean", "rho_postmean")) {
  row <- sprintf("  %-22s", k)
  for (lab in present) row <- paste0(row, sprintf("  %12.4f", mets[[lab]][[k]]))
  cat(row, "\n")
}

section("Per-smooth tau2_s (smaller = more shrinkage = smoother fit)")
for (j in 1:3) {
  k <- sprintf("tau2_s_main_%d", j)
  row <- sprintf("  %-22s", k)
  for (lab in present) row <- paste0(row, sprintf("  %12.5f", mets[[lab]][[k]]))
  cat(row, "\n")
}
for (k in int_keys) {
  metric <- sprintf("tau2_s_int_%s", k)
  row <- sprintf("  %-22s", metric)
  for (lab in present) {
    v <- mets[[lab]][[metric]]
    row <- paste0(row, sprintf("  %12s",
                               if (is.null(v)) "NA" else sprintf("%12.5f", v)))
  }
  cat(row, "\n")
}

section("Timing (minutes)")
row <- sprintf("  %-22s", "total_min")
for (lab in present) row <- paste0(row, sprintf("  %12.2f", mets[[lab]]$time_min))
cat(row, "\n")

# ---- save tidy CSVs ----
out_csv_dir <- file.path("comparison", "output")
dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)

# wide table: rows = metrics, cols = fits
tab_wide <- data.frame(metric = keys, stringsAsFactors = FALSE)
for (lab in present) {
  tab_wide[[col_labels[lab]]] <- vapply(keys,
    function(k) {
      v <- mets[[lab]][[k]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }, numeric(1))
}
write.csv(tab_wide,
          file.path(out_csv_dir, "scenarioB_h3_summary.csv"),
          row.names = FALSE)

# Per-smooth tau2_s table only (the headline diagnostic for H3)
tau_keys <- grep("^tau2_s_main_", keys, value = TRUE)
tau_tab  <- tab_wide[tab_wide$metric %in% tau_keys, , drop = FALSE]
write.csv(tau_tab,
          file.path(out_csv_dir, "scenarioB_h3_tau2s_table.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  comparison/output/scenarioB_h3_summary.csv\n")
cat("  comparison/output/scenarioB_h3_tau2s_table.csv\n\n")

# ---- decision logic ----
if (all(c("ortho", "h3conf") %in% present)) {
  rmse_f1_y <- mets$ortho$RMSE_f1
  rmse_f1_h <- mets$h3conf$RMSE_f1
  cov_f1_y  <- mets$ortho$cov_f1
  cov_f1_h  <- mets$h3conf$cov_f1
  cat("---- H3 CONFIRMATION VERDICT ----\n")
  cat(sprintf("  RMSE_f1: %.4f -> %.4f  (ratio %.2f)\n",
              rmse_f1_y, rmse_f1_h, rmse_f1_h / rmse_f1_y))
  cat(sprintf("  cov_f1 : %.4f -> %.4f\n", cov_f1_y, cov_f1_h))
  pass_rmse <- (rmse_f1_h <= 0.5 * rmse_f1_y)
  pass_cov  <- (cov_f1_h  >  cov_f1_y) && (cov_f1_h >= 0.85)
  if (pass_rmse && pass_cov) {
    cat("  -> H3 CONFIRMED. Undersmoothing on f1 is the cause.\n")
  } else if (pass_rmse && !pass_cov) {
    cat("  -> RMSE drops but coverage does not recover.\n")
    cat("     Possible: bands too narrow because clamp is too tight.\n")
  } else {
    cat("  -> H3 NOT CONFIRMED. Look for a different mechanism.\n")
  }
}

if (all(c("ortho", "h3tbs") %in% present)) {
  rmse_f1_y <- mets$ortho$RMSE_f1
  rmse_f1_t <- mets$h3tbs$RMSE_f1
  rmse_f2_y <- mets$ortho$RMSE_f2
  rmse_f2_t <- mets$h3tbs$RMSE_f2
  rmse_f3_y <- mets$ortho$RMSE_f3
  rmse_f3_t <- mets$h3tbs$RMSE_f3
  cat("\n---- H3 PRIOR-FIX VERDICT (tighten b_s) ----\n")
  cat(sprintf("  RMSE_f1: %.4f -> %.4f  (ratio %.2f)\n",
              rmse_f1_y, rmse_f1_t, rmse_f1_t / rmse_f1_y))
  cat(sprintf("  RMSE_f2: %.4f -> %.4f  (ratio %.2f)\n",
              rmse_f2_y, rmse_f2_t, rmse_f2_t / rmse_f2_y))
  cat(sprintf("  RMSE_f3: %.4f -> %.4f  (ratio %.2f)\n",
              rmse_f3_y, rmse_f3_t, rmse_f3_t / rmse_f3_y))
  drop_f1 <- (rmse_f1_t < rmse_f1_y)
  noregr  <- (rmse_f2_t <= 1.10 * rmse_f2_y) && (rmse_f3_t <= 1.10 * rmse_f3_y)
  if (drop_f1 && noregr) {
    cat("  -> SUCCESS: tighten-b_s fixes f1 without regressing f2/f3.\n")
  } else if (drop_f1 && !noregr) {
    cat("  -> PARTIAL: f1 improves but f2/f3 over-shrunk.\n")
    cat("     Recommendation: half-Cauchy / PC prior on sigma_s.\n")
  } else {
    cat("  -> FAILURE: tighten-b_s did not fix f1.\n")
  }
}

cat("\nDone.\n")
