# ============================================================
# baby_bayes.R
# Stage A: Closed-form Bayesian posterior for eta = (mu, beta)
# with FIXED (rho, sigma2, tau2)
#
# Depends on:
#   - spatial_utils.R:  pairdist(), matern_cor(), solve_chol()
#   - ls_basis.R:       ls_additive_build()
#   - fit_spatial_reml.R: fit_ls_spatial()  [for REML comparison]
#   - marginal_utils.R: default_truth_f_list(), etc
# ============================================================

# ------------------------------------------------------------
# baby_bayes_fit()
# 
# Compute closed-form Gaussian posterior for eta = (mu, beta)
# given FIXED covariance Sigma.
#
# Prior:  eta ~ N(0, kappa2 * I_p)
# Likelihood: y | eta ~ N(H eta, Sigma)
# Posterior:  eta | y ~ N(eta_bar, V)
#
# INPUTS:
#   y      : n-vector
#   H      : n x p design matrix (= cbind(1, W))
#   Sigma  : n x n covariance matrix (= sigma2*R + tau2*I)
#   kappa2 : prior variance (scalar); large = vague
#   jitter : diagonal jitter for Cholesky stability
#
# RETURNS:
#   eta_bar  : p-vector, posterior mean
#   V        : p x p posterior covariance
#   eta_se   : p-vector, posterior std dev = sqrt(diag(V))
#   ci_lower : 95% credible interval lower bounds
#   ci_upper : 95% credible interval upper bounds
#   Q        : posterior precision matrix
# ------------------------------------------------------------
baby_bayes_fit <- function(y, H, Sigma, kappa2 = 1e6, jitter = 1e-8) {
  
  y <- as.numeric(y)
  H <- as.matrix(H)
  n <- length(y)
  p <- ncol(H)
  
  stopifnot(nrow(H) == n)
  stopifnot(nrow(Sigma) == n && ncol(Sigma) == n)
  
  # --- Cholesky of Sigma ---
  L <- chol(Sigma + diag(jitter, n))   # upper tri: Sigma ≈ t(L) %*% L
  
  # Whiten
  y_w <- forwardsolve(t(L), y)
  H_w <- forwardsolve(t(L), H)
  
  # --- Posterior precision and covariance ---
  HtSinvH <- crossprod(H_w)            # H' Sigma^{-1} H
  Lambda0_inv <- diag(1 / kappa2, p)   # prior precision
  
  Q <- HtSinvH + Lambda0_inv           # posterior precision
  U_Q <- chol(Q)
  V <- chol2inv(U_Q)                   # posterior covariance
  
  # --- Posterior mean ---
  HtSinv_y <- crossprod(H_w, y_w)      # H' Sigma^{-1} y
  eta_bar <- as.vector(V %*% HtSinv_y)
  
  # --- Credible intervals ---
  eta_se   <- sqrt(diag(V))
  ci_lower <- eta_bar - 1.96 * eta_se
  ci_upper <- eta_bar + 1.96 * eta_se
  
  list(
    eta_bar  = eta_bar,
    V        = V,
    eta_se   = eta_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    Q        = Q
  )
}


# ------------------------------------------------------------
# bayes_marginal_band()
#
# Compute pointwise posterior mean + 95% credible band
# for marginal f_j(x_grid).
#
# INPUTS:
#   bb     : output from baby_bayes_fit()
#   obj    : output from fit_ls_spatial()
#   j      : covariate index (1..p)
#   x_grid : grid of x values
#   clip   : clip to knot range?
#   level  : credible level (default 0.95)
#
# RETURNS:
#   data.frame with columns: x, f_hat, f_se, lower, upper
# ------------------------------------------------------------
bayes_marginal_band <- function(bb, obj, j, x_grid, clip = TRUE, level = 0.95) {
  
  z <- qnorm(1 - (1 - level) / 2)
  
  # Build W_grid for covariate j
  W_grid_j <- obj$des$objs[[j]]$design_new(x_grid, type = "W", clip = clip)
  
  # Which columns in the full eta correspond to covariate j?
  # +1 because eta[1] is intercept (mu)
  col_idx <- 1 + obj$des$col_map[[j]]
  
  # Posterior mean of f_j(x_grid)
  beta_j_bar <- bb$eta_bar[col_idx]
  f_hat <- as.vector(W_grid_j %*% beta_j_bar)
  
  # Posterior variance: diag(W_j V_jj W_j')
  V_j <- bb$V[col_idx, col_idx, drop = FALSE]
  f_var <- rowSums((W_grid_j %*% V_j) * W_grid_j)
  f_se  <- sqrt(pmax(f_var, 0))
  
  data.frame(
    x     = x_grid,
    f_hat = f_hat,
    f_se  = f_se,
    lower = f_hat - z * f_se,
    upper = f_hat + z * f_se
  )
}


# ------------------------------------------------------------
# plot_bayes_marginals()
#
# Plot marginal curves with credible bands.
# Overlay truth if provided.
#
# INPUTS:
#   bb            : output from baby_bayes_fit()
#   obj           : output from fit_ls_spatial()
#   x_grid        : evaluation grid
#   truth_f_list  : list of true functions (optional)
#   var_names     : variable names
#   X_raw         : raw covariate matrix (for centering truth)
#   level         : credible level
# ------------------------------------------------------------
plot_bayes_marginals <- function(bb, obj, x_grid = seq(0, 1, length.out = 101),
                                  truth_f_list = NULL, var_names = NULL,
                                  X_raw = NULL, level = 0.95,
                                  ylim_shared = TRUE, symmetric = TRUE, pad = 0.05) {
  
  p <- length(obj$des$objs)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  # Compute all bands first (for shared ylim)
  bands <- vector("list", p)
  for (j in 1:p) {
    bands[[j]] <- bayes_marginal_band(bb, obj, j, x_grid, level = level)
  }
  
  # Compute truth (centered)
  truth_curves <- NULL
  if (!is.null(truth_f_list)) {
    truth_curves <- vector("list", p)
    for (j in 1:p) {
      ytrue <- truth_f_list[[j]](x_grid)
      if (!is.null(X_raw)) {
        shift <- mean(truth_f_list[[j]](X_raw[, j]))
        ytrue <- ytrue - shift
      }
      truth_curves[[j]] <- ytrue
    }
  }
  
  # Shared ylim
  ylim <- NULL
  if (ylim_shared) {
    all_y <- unlist(lapply(bands, function(b) c(b$lower, b$upper)))
    if (!is.null(truth_curves)) all_y <- c(all_y, unlist(truth_curves))
    rng <- range(all_y, finite = TRUE)
    if (symmetric) {
      a <- max(abs(rng))
      rng <- c(-a, a)
    }
    span <- diff(rng); if (span == 0) span <- 1
    rng <- rng + c(-1, 1) * pad * span
    ylim <- rng
  }
  
  for (j in 1:p) {
    b <- bands[[j]]
    
    plot(b$x, b$f_hat, type = "n",
         main = paste0("Marginal: ", var_names[j]),
         xlab = "x", ylab = paste0("f_", j, "(x)"),
         ylim = ylim)
    
    # Credible band
    polygon(c(b$x, rev(b$x)),
            c(b$lower, rev(b$upper)),
            col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
    
    # Posterior mean
    lines(b$x, b$f_hat, lwd = 2, col = "blue")
    
    # Truth
    if (!is.null(truth_curves)) {
      lines(b$x, truth_curves[[j]], lty = 2, col = "red", lwd = 1.5)
      
      # Coverage: fraction of grid points where truth is inside band
      inside <- truth_curves[[j]] >= b$lower & truth_curves[[j]] <= b$upper
      cov_pct <- mean(inside) * 100
      
      # RMSE
      r <- sqrt(mean((b$f_hat - truth_curves[[j]])^2))
      mtext(sprintf("RMSE=%.4f  band cov=%.0f%%", r, cov_pct),
            side = 3, line = 0.2, cex = 0.8)
    }
    
    abline(h = 0, lty = 3, col = "gray50")
    
    leg_labels <- c("posterior mean", paste0(round(level*100), "% band"))
    leg_col    <- c("blue", rgb(0.2, 0.4, 0.8, 0.4))
    leg_lty    <- c(1, NA)
    leg_lwd    <- c(2, 10)
    if (!is.null(truth_curves)) {
      leg_labels <- c(leg_labels, "truth")
      leg_col    <- c(leg_col, "red")
      leg_lty    <- c(leg_lty, 2)
      leg_lwd    <- c(leg_lwd, 1.5)
    }
    legend("topleft", legend = leg_labels, col = leg_col,
           lty = leg_lty, lwd = leg_lwd, bty = "n", cex = 0.85)
  }
}


# ------------------------------------------------------------
# bayes_curve_coverage()
#
# For simulation: compute pointwise coverage of credible band
# over the grid for each covariate.
#
# Returns data.frame: var, rmse_curve, mean_coverage, min_coverage
# ------------------------------------------------------------
bayes_curve_coverage <- function(bb, obj, x_grid, truth_f_list,
                                  X_raw = NULL, var_names = NULL,
                                  level = 0.95) {
  p <- length(obj$des$objs)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  rmse_vec <- numeric(p)
  mean_cov <- numeric(p)
  
  for (j in 1:p) {
    b <- bayes_marginal_band(bb, obj, j, x_grid, level = level)
    
    ytrue <- truth_f_list[[j]](x_grid)
    if (!is.null(X_raw)) {
      shift <- mean(truth_f_list[[j]](X_raw[, j]))
      ytrue <- ytrue - shift
    }
    
    rmse_vec[j] <- sqrt(mean((b$f_hat - ytrue)^2))
    inside <- ytrue >= b$lower & ytrue <= b$upper
    mean_cov[j] <- mean(inside)
  }
  
  data.frame(
    var       = var_names,
    rmse_curve = rmse_vec,
    band_coverage = mean_cov
  )
}
