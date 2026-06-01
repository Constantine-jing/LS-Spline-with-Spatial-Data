# ============================================================
# compute_recovery_metrics.R
# 
# Method-agnostic recovery metrics for cross-method comparison.
# Consumes ONLY the canonical fit_result schema and the truth
# (from simulate_scenario_*).
# 
# Per-fit metrics produced:
#   - rmse_f_j      : RMSE of posterior-mean f_j on the eval grid
#   - cov95_f_j     : pointwise 95%-coverage of f_j on the grid
#   - rmse_s        : RMSE of posterior-mean s(loc) at observed locs
#   - cor_s         : correlation of posterior-mean s with truth
#   - bias_sigma2   : posterior mean of sigma2 minus truth
#   - rmse_sigma2   : sqrt(MSE) for sigma2 across draws around truth
#   - cov95_sigma2  : 1 if true sigma2 in 95% CI, else 0
#   - same for tau2_s, rho (where defined)
#   - timing total/fit/post seconds
# 
# Output: a one-row data.frame, suitable for rbind into a tidy
# results table.
# ============================================================


# ------------------------------------------------------------
# helper: pointwise 95% credible/confidence interval from draws
# ------------------------------------------------------------
.pointwise_ci <- function(draws, probs = c(0.025, 0.975)) {
  apply(draws, 1, quantile, probs = probs, na.rm = TRUE)
}


# ------------------------------------------------------------
# compute_recovery_metrics(fit, sim)
# 
# Args:
#   fit : canonical fit_result (already validated)
#   sim : object from simulate_scenario_*; must contain
#         truth_f_grid, truth_s_obs, truth_params, x_grid_1d
# 
# Returns:
#   one-row data.frame with method, scenario, seed, and all metrics
# ------------------------------------------------------------
compute_recovery_metrics <- function(fit, sim) {
  
  stopifnot(identical(fit$scenario, sim$scenario))
  stopifnot(identical(fit$seed,     sim$seed))
  
  p <- length(sim$truth_f_grid)
  
  out <- data.frame(
    method   = fit$method,
    scenario = fit$scenario,
    seed     = fit$seed,
    n_draws  = ncol(fit$f_main[[1]]),
    stringsAsFactors = FALSE
  )
  
  # ---- main-effect smooths: per-smooth RMSE and coverage ----
  for (j in seq_len(p)) {
    truth_j <- sim$truth_f_grid[[j]]
    draws_j <- fit$f_main[[j]]
    pmean_j <- rowMeans(draws_j)
    qq      <- .pointwise_ci(draws_j)
    rmse_j  <- sqrt(mean((pmean_j - truth_j)^2))
    cov_j   <- mean(truth_j >= qq[1, ] & truth_j <= qq[2, ])
    out[[paste0("rmse_f",  j)]] <- rmse_j
    out[[paste0("cov95_f", j)]] <- cov_j
  }
  # also report mean RMSE / mean coverage across smooths
  out$rmse_f_mean  <- mean(unlist(out[1, paste0("rmse_f",  seq_len(p))]))
  out$cov95_f_mean <- mean(unlist(out[1, paste0("cov95_f", seq_len(p))]))
  
  # ---- spatial RE recovery ----
  if (!is.null(fit$s_obs) && !is.null(sim$truth_s_obs)) {
    s_pmean <- rowMeans(fit$s_obs)
    out$rmse_s <- sqrt(mean((s_pmean - sim$truth_s_obs)^2))
    out$cor_s  <- cor(s_pmean, sim$truth_s_obs)
  } else {
    out$rmse_s <- NA_real_
    out$cor_s  <- NA_real_
  }
  
  # ---- variance-component recovery ----
  # bias = E[posterior] - truth
  # cov95 = 1 if truth in 95% CI of posterior draws, else 0
  vc_score <- function(x_draws, truth_val) {
    if (is.null(x_draws) || all(is.na(x_draws))) {
      return(c(post_mean = NA_real_, bias = NA_real_, cov95 = NA_real_))
    }
    pm <- mean(x_draws)
    ci <- quantile(x_draws, c(0.025, 0.975), na.rm = TRUE)
    cov_in <- as.integer(truth_val >= ci[1] && truth_val <= ci[2])
    c(post_mean = pm, bias = pm - truth_val, cov95 = cov_in)
  }
  
  scores_sigma2 <- vc_score(fit$var_comp$sigma2, sim$truth_params$sigma2)
  out$post_mean_sigma2 <- scores_sigma2["post_mean"]
  out$bias_sigma2      <- scores_sigma2["bias"]
  out$cov95_sigma2     <- scores_sigma2["cov95"]
  
  scores_tau2_s <- vc_score(fit$var_comp$tau2_s, sim$truth_params$tau2_s)
  out$post_mean_tau2_s <- scores_tau2_s["post_mean"]
  out$bias_tau2_s      <- scores_tau2_s["bias"]
  out$cov95_tau2_s     <- scores_tau2_s["cov95"]
  
  # rho is reported in each method's NATIVE parameterization, so
  # recovery against the truth parameter is only meaningful for
  # ours and (with conversion) for INLA. We report posterior mean
  # only; the per-method comparison in the paper uses effective
  # range, computed downstream.
  out$post_mean_rho <- if (!is.null(fit$var_comp$rho))
    mean(fit$var_comp$rho) else NA_real_
  
  # ---- timing ----
  out$total_sec <- fit$timing$total_sec
  out$fit_sec   <- fit$timing$fit_sec
  out$post_sec  <- fit$timing$post_sec
  
  out
}


# ------------------------------------------------------------
# summarize_metrics(metric_df)
# 
# Aggregate per-seed metrics into mean (+/- Monte Carlo SE)
# per (method, scenario). For coverage, MC SE of a Bernoulli mean.
# 
# Returns a long-form summary data.frame.
# ------------------------------------------------------------
summarize_metrics <- function(metric_df) {
  
  stopifnot(all(c("method", "scenario", "seed") %in% colnames(metric_df)))
  
  # determine which columns are metrics (numeric, not seed/n_draws)
  excluded <- c("method", "scenario", "seed", "n_draws")
  metric_cols <- setdiff(colnames(metric_df), excluded)
  metric_cols <- metric_cols[sapply(metric_df[metric_cols], is.numeric)]
  
  # split
  splits <- split(metric_df, list(metric_df$method, metric_df$scenario), drop = TRUE)
  
  summary_rows <- list()
  for (key in names(splits)) {
    df_g <- splits[[key]]
    method_g   <- df_g$method[1]
    scenario_g <- df_g$scenario[1]
    n_seeds <- nrow(df_g)
    
    row <- data.frame(method = method_g, scenario = scenario_g,
                      n_seeds = n_seeds,
                      stringsAsFactors = FALSE)
    for (col in metric_cols) {
      v   <- df_g[[col]]
      v   <- v[is.finite(v)]
      m   <- mean(v)
      mse <- if (length(v) > 1L) sd(v) / sqrt(length(v)) else NA_real_
      row[[paste0("mean_",  col)]] <- m
      row[[paste0("mcse_",  col)]] <- mse
    }
    summary_rows[[key]] <- row
  }
  
  do.call(rbind, summary_rows)
}
