# ============================================================
# run_michigan_lung_poisson.R
# Real data application: Poisson spatial additive model for
# Michigan county-level lung cancer mortality.
#
# Model:
#   y_i ~ Poisson(E_i * theta_i)
#   log(theta_i) = mu + sum_j f_j(x_{ij}) + b_i
#   b_i ~ GP(0, sigma2 * R(rho, nu))
#
# Equivalently:
#   log(lambda_i) = log(E_i) + mu + sum_j f_j(x_{ij}) + b_i
# where log(E_i) is the offset.
#
# This parallels the Munich rent Gaussian application but for
# Poisson count outcomes.  Comparison with Nandy et al. (2017)
# who used Gaussian model on Tukey-transformed rates.
#
# Depends on: ls_basis.R, spatial_utils.R, gibbs_stage_c_poisson.R
# Data:       michigan_lung_data.rda (from prep_michigan_lung_data.R)
#
# Run: Rscript run_michigan_lung_poisson.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("gibbs_stage_c_poisson.R")

# ============================================================
# Settings
# ============================================================
M_C    <- 10          # knots per covariate (fewer counties => fewer knots)
n_iter <- 15000
n_burn <- 5000
n_thin <- 1
nu     <- 1.5

x_grid <- seq(0, 1, length.out = 101)

dir.create("michigan_lung", showWarnings = FALSE)


# ============================================================
# Load data
# ============================================================
cat("=== Michigan Lung Cancer: Poisson Spatial Additive Model ===\n")

if (!file.exists("michigan_lung_data.rda")) {
  stop("michigan_lung_data.rda not found!\n",
       "  Run prep_michigan_lung_data.R first to assemble the data.\n",
       "  See that script for instructions on downloading from CDC WONDER.")
}

load("michigan_lung_data.rda")  # loads `lung_data`

cat(sprintf("  Loaded: %d counties\n", nrow(lung_data)))
cat(sprintf("  Deaths: mean=%.1f, range=[%d, %d]\n",
            mean(lung_data$Deaths), min(lung_data$Deaths),
            max(lung_data$Deaths)))


# ============================================================
# Identify covariates
# ============================================================
# Determine which covariates are available in the dataset.
# We look for the Census variables that Nandy et al. used.
# Start with whatever we have; the model adapts to the
# number of covariates.

potential_covs <- c("poverty", "unemployment", "pct_urban",
                    "pct_65plus", "median_income", "pct_nonwhite",
                    "pct_never_married", "pct_agriculture",
                    "pct_white_collar", "pct_highschool",
                    "pct_under18", "pct_crowding",
                    "pct_foreign_born", "pct_lang_isolation",
                    "pct_same_house", "pct_move_same_county",
                    "pct_move_same_state", "pct_move_diff_state")

available_covs <- intersect(potential_covs, names(lung_data))

if (length(available_covs) == 0) {
  cat("\n  WARNING: No Census covariates found in dataset.\n")
  cat("  Running intercept-only spatial model (no additive terms).\n")
  cat("  This is equivalent to a BYM-type disease mapping model.\n")
  cat("  Add Census covariates to michigan_lung_data.rda for the\n")
  cat("  full geoadditive model.\n\n")
  
  run_additive <- FALSE
} else {
  cat(sprintf("\n  Found %d covariates: %s\n",
              length(available_covs),
              paste(available_covs, collapse = ", ")))
  run_additive <- TRUE
}


# ============================================================
# Prepare spatial structures
# ============================================================
coords <- as.matrix(lung_data[, c("lon", "lat")])

# Rescale coordinates to [0,1] x [0,1]
coord_range <- apply(coords, 2, range)
coords_scaled <- coords
coords_scaled[, 1] <- (coords[, 1] - coord_range[1, 1]) /
  (coord_range[2, 1] - coord_range[1, 1])
coords_scaled[, 2] <- (coords[, 2] - coord_range[1, 2]) /
  (coord_range[2, 2] - coord_range[1, 2])

# Add small jitter for any exact coordinate ties
# (unlikely for county centroids, but safe)
set.seed(999)
coords_scaled <- coords_scaled + matrix(rnorm(nrow(coords_scaled) * 2,
                                               sd = 0.001),
                                         ncol = 2)
coords_scaled <- pmax(pmin(coords_scaled, 1), 0)

D <- as.matrix(dist(coords_scaled))
n <- nrow(lung_data)

cat(sprintf("  n = %d counties\n", n))
cat(sprintf("  Distance range: [%.4f, %.4f]\n", min(D[D > 0]), max(D)))


# ============================================================
# Build design matrix
# ============================================================
y      <- lung_data$Deaths
offset <- log(lung_data$Expected)

if (run_additive) {
  # Extract and rescale covariates to [0, 1]
  X_raw <- as.matrix(lung_data[, available_covs])
  p_covs <- ncol(X_raw)
  
  X <- X_raw
  for (j in 1:p_covs) {
    xj <- X_raw[, j]
    xj_range <- range(xj, na.rm = TRUE)
    if (xj_range[2] > xj_range[1]) {
      X[, j] <- (xj - xj_range[1]) / (xj_range[2] - xj_range[1])
    } else {
      X[, j] <- 0.5  # constant covariate — shouldn't happen
    }
  }
  colnames(X) <- available_covs
  
  # Adjust knots: with n~68 counties, can't use too many knots
  # Use min(M_C, floor(n/4)) to be safe
  M_use <- min(M_C, floor(n / 4))
  cat(sprintf("  Using M = %d knots per covariate (n=%d, p=%d)\n",
              M_use, n, p_covs))
  
  # Build LS basis
  des <- ls_additive_build(X, M_vec = M_use)
  H   <- cbind(1, des$W)
  col_map <- des$col_map
  var_names <- available_covs
  
} else {
  # Intercept-only model
  H <- matrix(1, n, 1)
  col_map <- list()  # empty — no spline terms
  p_covs <- 0
  var_names <- character(0)
  
  # For the sampler, we need at least one col_map entry.
  # Use a dummy approach: the sampler handles p_covs = 0
  # by skipping the tau2_s updates.
  # Actually, our sampler requires col_map to be non-empty.
  # Workaround: skip additive model and just run spatial-only.
  cat("  (Intercept-only model — spatial random effect only)\n")
}

cat(sprintf("  H dimensions: %d x %d\n", nrow(H), ncol(H)))


# ============================================================
# Run Poisson Gibbs sampler
# ============================================================
cat("\n--- Running Poisson spatial additive sampler ---\n")
t0 <- proc.time()

if (run_additive) {
  gs <- gibbs_poisson_sampler(
    y       = y,
    H       = H,
    D       = D,
    nu      = nu,
    col_map = col_map,
    offset  = offset,
    n_iter  = n_iter,
    n_burn  = n_burn,
    n_thin  = n_thin,
    a_sigma = 2, b_sigma = 0.5,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = -1.6,
    log_rho_sd = 1.0,
    iwls_scale_eta = 1.0,
    iwls_scale_b   = 1.0,
    mh_sd_log_rho  = 0.2,
    verbose = TRUE
  )
} else {
  # Intercept-only: need a modified call
  # Use a single dummy column in col_map
  # (the intercept is col 1, no spline blocks)
  # For now, just warn and skip
  cat("  Intercept-only model requires a modified sampler setup.\n")
  cat("  Please add Census covariates for the full geoadditive model.\n")
  stop("Add covariates to proceed.")
}

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Elapsed: %.1f seconds (%.1f minutes)\n",
            elapsed, elapsed / 60))


# ============================================================
# Posterior summaries
# ============================================================
cat("\n--- Posterior summaries ---\n")
cat(sprintf("  sigma2: mean=%.4f  median=%.4f  95%%CI=[%.4f, %.4f]\n",
            mean(gs$sigma2_samples),
            median(gs$sigma2_samples),
            quantile(gs$sigma2_samples, 0.025),
            quantile(gs$sigma2_samples, 0.975)))
cat(sprintf("  rho:    mean=%.4f  median=%.4f  95%%CI=[%.4f, %.4f]\n",
            mean(gs$rho_samples),
            median(gs$rho_samples),
            quantile(gs$rho_samples, 0.025),
            quantile(gs$rho_samples, 0.975)))

tau2_s_means <- colMeans(gs$tau2_s_samples)
cat("\n  Smoothing variances (tau2_s):\n")
for (j in seq_along(var_names)) {
  cat(sprintf("    %-25s: %.4f\n", var_names[j], tau2_s_means[j]))
}


# ============================================================
# MH acceptance rates
# ============================================================
cat(sprintf("\n  Acceptance rates: eta=%.3f  b=%.3f  rho=%.3f\n",
            gs$accept_rate["eta"],
            gs$accept_rate["b"],
            gs$accept_rate["rho"]))


# ============================================================
# Trace plots
# ============================================================
pdf("michigan_lung/traces.pdf", width = 10, height = 6)
plot_gibbs_trace_poisson(gs)
dev.off()
cat("  Saved: michigan_lung/traces.pdf\n")


# ============================================================
# Marginal covariate effect plots
# ============================================================
if (run_additive && length(var_names) > 0) {
  
  n_plots <- length(var_names)
  n_row <- ceiling(sqrt(n_plots))
  n_col <- ceiling(n_plots / n_row)
  
  pdf("michigan_lung/marginal_effects.pdf",
      width = 4 * n_col, height = 4 * n_row)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))
  
  eta_post_mean <- colMeans(gs$eta_samples)
  
  for (j in seq_along(var_names)) {
    idx_col <- 1 + col_map[[j]]
    W_grid_j <- des$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)
    
    # Posterior samples of f_j
    f_samples <- gs$eta_samples[, idx_col, drop = FALSE] %*% t(W_grid_j)
    f_samples <- f_samples - rowMeans(f_samples)
    
    f_mean <- colMeans(f_samples)
    f_lo   <- apply(f_samples, 2, quantile, 0.025)
    f_hi   <- apply(f_samples, 2, quantile, 0.975)
    
    ylim <- range(c(f_lo, f_hi))
    
    # Map x_grid back to original scale for x-axis labels
    xj_raw <- X_raw[, j]
    x_orig <- x_grid * (max(xj_raw) - min(xj_raw)) + min(xj_raw)
    
    plot(x_orig, f_mean, type = "l", col = "blue", lwd = 2,
         ylim = ylim,
         xlab = var_names[j],
         ylab = paste0("f(", var_names[j], ")"),
         main = var_names[j])
    polygon(c(x_orig, rev(x_orig)),
            c(f_lo, rev(f_hi)),
            col = rgb(0.3, 0.5, 0.9, 0.3), border = NA)
    abline(h = 0, lty = 2, col = "gray40")
    
    # Add rug of observed data points
    rug(xj_raw, col = rgb(0, 0, 0, 0.3))
  }
  
  dev.off()
  cat("  Saved: michigan_lung/marginal_effects.pdf\n")
}


# ============================================================
# Spatial random effect map
# ============================================================
b_post_mean <- colMeans(gs$b_samples)
b_post_lo   <- apply(gs$b_samples, 2, quantile, 0.025)
b_post_hi   <- apply(gs$b_samples, 2, quantile, 0.975)

# Posterior relative risk from spatial effect: exp(b)
RR_spatial <- exp(b_post_mean)

pdf("michigan_lung/spatial_effect.pdf", width = 12, height = 5)
par(mfrow = c(1, 3))

# Map 1: Observed SMR
color_breaks <- seq(0, max(c(lung_data$SMR, RR_spatial)) * 1.1,
                     length.out = 101)
obs_cols <- hcl.colors(100, "RdYlBu", rev = TRUE)[
  cut(lung_data$SMR, breaks = color_breaks, labels = FALSE)
]
plot(coords[, 1], coords[, 2], col = obs_cols,
     pch = 16, cex = 2.5,
     xlab = "Longitude", ylab = "Latitude",
     main = "Observed SMR")

# Map 2: Posterior spatial relative risk exp(b)
rr_cols <- hcl.colors(100, "RdYlBu", rev = TRUE)[
  cut(RR_spatial, breaks = color_breaks, labels = FALSE)
]
plot(coords[, 1], coords[, 2], col = rr_cols,
     pch = 16, cex = 2.5,
     xlab = "Longitude", ylab = "Latitude",
     main = "Spatial relative risk exp(b)")

# Map 3: Significance — counties where 95% CI for b excludes 0
sig <- ifelse(b_post_lo > 0, "high",
              ifelse(b_post_hi < 0, "low", "ns"))
sig_cols <- ifelse(sig == "high", "firebrick",
                   ifelse(sig == "low", "steelblue", "gray80"))
plot(coords[, 1], coords[, 2], col = sig_cols,
     pch = 16, cex = 2.5,
     xlab = "Longitude", ylab = "Latitude",
     main = "Spatial effect significance")
legend("bottomright",
       legend = c("Elevated (95% CI > 0)", "Reduced (95% CI < 0)", "Not sig."),
       col = c("firebrick", "steelblue", "gray80"),
       pch = 16, cex = 0.8)

dev.off()
cat("  Saved: michigan_lung/spatial_effect.pdf\n")


# ============================================================
# Model comparison: fitted vs observed
# ============================================================
eta_post_mean_full <- colMeans(gs$eta_samples)
linpred_post <- as.vector(H %*% eta_post_mean_full) + b_post_mean + offset
lambda_fitted <- exp(linpred_post)

pdf("michigan_lung/fitted_vs_observed.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Fitted vs Observed counts
plot(y, lambda_fitted,
     xlab = "Observed deaths", ylab = "Fitted E[Y]",
     main = "Fitted vs Observed",
     pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.5))
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Residuals
pearson_resid <- (y - lambda_fitted) / sqrt(lambda_fitted)
plot(lambda_fitted, pearson_resid,
     xlab = "Fitted E[Y]", ylab = "Pearson residual",
     main = "Pearson residuals",
     pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.5))
abline(h = 0, col = "red", lty = 2)
abline(h = c(-2, 2), col = "gray60", lty = 3)

dev.off()
cat("  Saved: michigan_lung/fitted_vs_observed.pdf\n")


# ============================================================
# Summary output
# ============================================================
cat("\n\n#########################################\n")
cat("# Michigan Lung Cancer — Model Summary\n")
cat("#########################################\n\n")

cat(sprintf("  Counties: %d\n", n))
cat(sprintf("  Covariates: %d (%s)\n", length(var_names),
            paste(var_names, collapse = ", ")))
cat(sprintf("  Knots per covariate: %d\n", M_use))
cat(sprintf("  MCMC: %d iterations, %d burn-in\n", n_iter, n_burn))
cat(sprintf("  Elapsed: %.1f minutes\n\n", elapsed / 60))

cat(sprintf("  sigma2 (spatial): %.4f [%.4f, %.4f]\n",
            mean(gs$sigma2_samples),
            quantile(gs$sigma2_samples, 0.025),
            quantile(gs$sigma2_samples, 0.975)))
cat(sprintf("  rho (range):      %.4f [%.4f, %.4f]\n",
            mean(gs$rho_samples),
            quantile(gs$rho_samples, 0.025),
            quantile(gs$rho_samples, 0.975)))
cat(sprintf("  Intercept:        %.4f [%.4f, %.4f]\n",
            mean(gs$eta_samples[, 1]),
            quantile(gs$eta_samples[, 1], 0.025),
            quantile(gs$eta_samples[, 1], 0.975)))

cat(sprintf("\n  Deviance: -2*loglik = %.1f\n",
            -2 * sum(dpois(y, lambda_fitted, log = TRUE))))

# Pearson GOF
pearson_X2 <- sum(pearson_resid^2)
cat(sprintf("  Pearson X2 = %.1f  (df ~ %d)\n", pearson_X2, n - ncol(H) - 1))

cat("\n=== Done ===\n")
