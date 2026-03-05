# ============================================================
# gibbs_bayes.R
# Stage B: Gibbs sampler for eta, b, sigma2, tau2
# with FIXED rho (and fixed nu)
#
# Model:
#   y = H eta + b + eps
#   b ~ N(0, sigma2 * R)       R = matern_cor(D, rho, nu)
#   eps ~ N(0, tau2 * I_n)
#
# Priors (Bundle B):
#   eta ~ N(0, kappa2 * I_p)          [vague]
#   sigma2 ~ IG(a_sigma, b_sigma)     [weakly informative]
#   tau2   ~ IG(a_tau, b_tau)         [weakly informative]
#
# Gibbs steps (all conjugate):
#   1) eta    | y, b, sigma2, tau2    ~ MVN
#   2) b      | y, eta, sigma2, tau2  ~ MVN
#   3) sigma2 | b                     ~ IG
#   4) tau2   | y, eta, b             ~ IG
#
# Depends on: spatial_utils.R, ls_basis.R, baby_bayes.R
# ============================================================

# ------------------------------------------------------------
# gibbs_sampler()
#
# INPUTS:
#   y       : n-vector
#   H       : n x p design (cbind(1, W))
#   R       : n x n Matern correlation matrix (fixed rho, nu)
#   n_iter  : total MCMC iterations
#   n_burn  : burn-in (discarded)
#   n_thin  : thinning interval
#   kappa2  : prior variance for eta
#   a_sigma, b_sigma : IG hyperparams for sigma2
#   a_tau, b_tau     : IG hyperparams for tau2
#   init    : optional list with starting values
#   verbose : print progress?
#
# RETURNS:
#   list with:
#     eta_samples   : (n_keep x p) matrix
#     b_samples     : (n_keep x n) matrix
#     sigma2_samples: n_keep vector
#     tau2_samples  : n_keep vector
#     accept_info   : diagnostics
# ------------------------------------------------------------
gibbs_sampler <- function(y, H, R,
                          n_iter  = 5000,
                          n_burn  = 1000,
                          n_thin  = 1,
                          kappa2  = 1e6,
                          a_sigma = 1, b_sigma = 0.005,
                          a_tau   = 1, b_tau   = 0.005,
                          init    = NULL,
                          jitter  = 1e-8,
                          verbose = TRUE) {
  
  y <- as.numeric(y)
  H <- as.matrix(H)
  R <- as.matrix(R)
  n <- length(y)
  p <- ncol(H)
  
  # --- Pre-compute Cholesky of R (fixed throughout) ---
  L_R <- chol(R + diag(jitter, n))   # R ≈ t(L_R) %*% L_R
  R_inv_via_chol <- chol2inv(L_R)    # R^{-1}
  logdet_R <- 2 * sum(log(diag(L_R)))
  
  # Prior precision for eta
  Lambda0_inv <- diag(1 / kappa2, p)
  
  # --- Initialize ---
  if (is.null(init)) {
    sigma2 <- 1.0
    tau2   <- 0.5
    eta    <- rep(0, p)
    b      <- rep(0, n)
  } else {
    sigma2 <- init$sigma2
    tau2   <- init$tau2
    eta    <- init$eta
    b      <- init$b
  }
  
  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples    <- matrix(NA, n_keep, p)
  b_samples      <- matrix(NA, n_keep, n)
  sigma2_samples <- numeric(n_keep)
  tau2_samples   <- numeric(n_keep)
  
  keep_idx <- 0
  
  # ============================================================
  # MAIN GIBBS LOOP
  # ============================================================
  for (iter in 1:n_iter) {
    
    # ----------------------------------------------------------
    # Step 1: Draw eta | y, b, sigma2, tau2
    #
    # Likelihood for eta: (y - b) = H eta + eps, eps ~ N(0, tau2 I)
    # Prior: eta ~ N(0, kappa2 I)
    #
    # Posterior precision: Q_eta = (1/tau2) H'H + Lambda0_inv
    # Posterior mean:      m_eta = Q_eta^{-1} (1/tau2) H'(y - b)
    # ----------------------------------------------------------
    r_eta <- y - b   # "working response" for eta
    
    Q_eta <- (1 / tau2) * crossprod(H) + Lambda0_inv
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% ((1 / tau2) * crossprod(H, r_eta)))
    
    # Draw from MVN(m_eta, V_eta) using: eta = m_eta + L^{-T} z
    z_eta <- rnorm(p)
    eta   <- m_eta + as.vector(backsolve(U_eta, z_eta))
    
    # ----------------------------------------------------------
    # Step 2: Draw b | y, eta, sigma2, tau2
    #
    # Likelihood for b: (y - H eta) = b + eps, eps ~ N(0, tau2 I)
    # Prior: b ~ N(0, sigma2 R)
    #
    # Posterior precision: Q_b = (1/tau2) I + (1/sigma2) R^{-1}
    # Posterior mean:      m_b = Q_b^{-1} (1/tau2)(y - H eta)
    # ----------------------------------------------------------
    r_b <- y - H %*% eta   # "working response" for b
    
    Q_b <- (1 / tau2) * diag(n) + (1 / sigma2) * R_inv_via_chol
    U_b <- chol(Q_b + diag(jitter, n))
    
    m_b <- as.vector(chol2inv(U_b) %*% ((1 / tau2) * r_b))
    
    # Draw from MVN(m_b, Q_b^{-1})
    z_b <- rnorm(n)
    b   <- m_b + as.vector(backsolve(U_b, z_b))
    
    # ----------------------------------------------------------
    # Step 3: Draw sigma2 | b
    #
    # Prior: sigma2 ~ IG(a_sigma, b_sigma)
    # Likelihood: b ~ N(0, sigma2 R) => b' R^{-1} b / sigma2
    #
    # Posterior: sigma2 ~ IG(a_sigma + n/2,
    #                        b_sigma + b' R^{-1} b / 2)
    # ----------------------------------------------------------
    bRinvb <- as.numeric(t(b) %*% R_inv_via_chol %*% b)
    
    a_post_sigma <- a_sigma + n / 2
    b_post_sigma <- b_sigma + bRinvb / 2
    
    sigma2 <- 1 / rgamma(1, shape = a_post_sigma, rate = b_post_sigma)
    
    # ----------------------------------------------------------
    # Step 4: Draw tau2 | y, eta, b
    #
    # Prior: tau2 ~ IG(a_tau, b_tau)
    # Likelihood: eps = y - H eta - b ~ N(0, tau2 I)
    #
    # Posterior: tau2 ~ IG(a_tau + n/2,
    #                      b_tau + ||y - H eta - b||^2 / 2)
    # ----------------------------------------------------------
    resid <- y - H %*% eta - b
    rss   <- sum(resid^2)
    
    a_post_tau <- a_tau + n / 2
    b_post_tau <- b_tau + rss / 2
    
    tau2 <- 1 / rgamma(1, shape = a_post_tau, rate = b_post_tau)
    
    # ----------------------------------------------------------
    # Store (after burn-in, with thinning)
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]  <- eta
      b_samples[keep_idx, ]    <- b
      sigma2_samples[keep_idx] <- sigma2
      tau2_samples[keep_idx]   <- tau2
    }
    
    # Progress
    if (verbose && (iter %% 1000 == 0)) {
      cat(sprintf("  iter %d/%d  sigma2=%.4f  tau2=%.4f\n",
                  iter, n_iter, sigma2, tau2))
    }
  }
  
  list(
    eta_samples    = eta_samples,
    b_samples      = b_samples,
    sigma2_samples = sigma2_samples,
    tau2_samples   = tau2_samples,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter
  )
}


# ------------------------------------------------------------
# gibbs_summary()
#
# Compute posterior summaries from Gibbs output.
# Returns: posterior means, SDs, credible intervals
# for eta, sigma2, tau2.
# ------------------------------------------------------------
gibbs_summary <- function(gs, level = 0.95) {
  
  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  p <- ncol(gs$eta_samples)
  
  # --- eta summaries ---
  eta_mean <- colMeans(gs$eta_samples)
  eta_sd   <- apply(gs$eta_samples, 2, sd)
  eta_q    <- apply(gs$eta_samples, 2, quantile, probs = probs)
  
  # --- variance component summaries ---
  s2_summary <- c(mean   = mean(gs$sigma2_samples),
                  sd     = sd(gs$sigma2_samples),
                  quantile(gs$sigma2_samples, probs))
  
  t2_summary <- c(mean   = mean(gs$tau2_samples),
                  sd     = sd(gs$tau2_samples),
                  quantile(gs$tau2_samples, probs))
  
  list(
    eta_mean    = eta_mean,
    eta_sd      = eta_sd,
    eta_ci_lower = eta_q[1, ],
    eta_ci_upper = eta_q[3, ],
    sigma2 = s2_summary,
    tau2   = t2_summary
  )
}


# ------------------------------------------------------------
# gibbs_marginal_band()
#
# Compute posterior mean + credible band for f_j(x_grid)
# using ALL posterior samples (propagates uncertainty in
# eta AND variance components).
#
# This is the key improvement over Baby Bayes Stage A.
# ------------------------------------------------------------
gibbs_marginal_band <- function(gs, obj, j, x_grid, 
                                 X_raw = NULL, truth_f = NULL,
                                 clip = TRUE, level = 0.95) {
  
  alpha <- (1 - level) / 2
  
  # Build W_grid for covariate j
  W_grid_j <- obj$des$objs[[j]]$design_new(x_grid, type = "W", clip = clip)
  
  # Which columns in eta correspond to covariate j?
  col_idx <- 1 + obj$des$col_map[[j]]
  
  # For each posterior sample, compute f_j(x_grid)
  n_samples <- nrow(gs$eta_samples)
  n_grid    <- length(x_grid)
  f_matrix  <- matrix(NA, n_samples, n_grid)
  
  beta_j_samples <- gs$eta_samples[, col_idx, drop = FALSE]
  for (s in 1:n_samples) {
    f_matrix[s, ] <- as.vector(W_grid_j %*% beta_j_samples[s, ])
  }
  
  # Pointwise summaries
  f_hat  <- colMeans(f_matrix)
  f_se   <- apply(f_matrix, 2, sd)
  f_lower <- apply(f_matrix, 2, quantile, probs = alpha)
  f_upper <- apply(f_matrix, 2, quantile, probs = 1 - alpha)
  
  # Truth (centered)
  f_true <- NULL
  if (!is.null(truth_f)) {
    f_true <- truth_f(x_grid)
    if (!is.null(X_raw)) {
      f_true <- f_true - mean(truth_f(X_raw[, j]))
    }
  }
  
  data.frame(
    x      = x_grid,
    f_hat  = f_hat,
    f_se   = f_se,
    lower  = f_lower,
    upper  = f_upper,
    f_true = if (!is.null(f_true)) f_true else NA
  )
}


# ------------------------------------------------------------
# gibbs_curve_coverage()
#
# Coverage + RMSE table across all covariates.
# ------------------------------------------------------------
gibbs_curve_coverage <- function(gs, obj, x_grid, truth_f_list,
                                  X_raw = NULL, var_names = NULL,
                                  level = 0.95) {
  p <- length(obj$des$objs)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  rmse_vec <- numeric(p)
  cov_vec  <- numeric(p)
  
  for (j in 1:p) {
    band <- gibbs_marginal_band(gs, obj, j, x_grid,
                                 X_raw = X_raw, truth_f = truth_f_list[[j]],
                                 level = level)
    
    rmse_vec[j] <- sqrt(mean((band$f_hat - band$f_true)^2))
    inside <- band$f_true >= band$lower & band$f_true <= band$upper
    cov_vec[j] <- mean(inside)
  }
  
  data.frame(
    var       = var_names,
    rmse_curve = rmse_vec,
    band_coverage = cov_vec
  )
}


# ------------------------------------------------------------
# plot_gibbs_marginals()
#
# Plot marginal curves with Gibbs credible bands.
# Same interface as plot_bayes_marginals() but uses Gibbs samples.
# ------------------------------------------------------------
plot_gibbs_marginals <- function(gs, obj, x_grid = seq(0, 1, length.out = 101),
                                  truth_f_list = NULL, var_names = NULL,
                                  X_raw = NULL, level = 0.95,
                                  ylim_shared = TRUE, symmetric = TRUE, pad = 0.05) {
  
  p <- length(obj$des$objs)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  # Compute all bands
  bands <- vector("list", p)
  for (j in 1:p) {
    tf <- if (!is.null(truth_f_list)) truth_f_list[[j]] else NULL
    bands[[j]] <- gibbs_marginal_band(gs, obj, j, x_grid,
                                       X_raw = X_raw, truth_f = tf,
                                       level = level)
  }
  
  # Shared ylim
  ylim <- NULL
  if (ylim_shared) {
    all_y <- unlist(lapply(bands, function(b) c(b$lower, b$upper)))
    if (!is.null(truth_f_list)) {
      all_y <- c(all_y, unlist(lapply(bands, function(b) b$f_true)))
    }
    all_y <- all_y[is.finite(all_y)]
    rng <- range(all_y)
    if (symmetric) { a <- max(abs(rng)); rng <- c(-a, a) }
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
    polygon(c(b$x, rev(b$x)), c(b$lower, rev(b$upper)),
            col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
    lines(b$x, b$f_hat, lwd = 2, col = "blue")
    
    # Truth
    if (!all(is.na(b$f_true))) {
      lines(b$x, b$f_true, lty = 2, col = "red", lwd = 1.5)
      inside <- b$f_true >= b$lower & b$f_true <= b$upper
      cov_pct <- mean(inside) * 100
      r <- sqrt(mean((b$f_hat - b$f_true)^2))
      mtext(sprintf("RMSE=%.4f  band cov=%.0f%%", r, cov_pct),
            side = 3, line = 0.2, cex = 0.8)
    }
    
    abline(h = 0, lty = 3, col = "gray50")
    legend("topleft",
           legend = c("posterior mean", paste0(round(level*100), "% band"),
                       if (!all(is.na(b$f_true))) "truth" else NULL),
           col    = c("blue", rgb(0.2, 0.4, 0.8, 0.4),
                       if (!all(is.na(b$f_true))) "red" else NULL),
           lty    = c(1, NA, if (!all(is.na(b$f_true))) 2 else NULL),
           lwd    = c(2, 10, if (!all(is.na(b$f_true))) 1.5 else NULL),
           bty = "n", cex = 0.85)
  }
}


# ------------------------------------------------------------
# plot_gibbs_trace()
#
# Trace plots + density for sigma2 and tau2.
# Basic MCMC diagnostic.
# ------------------------------------------------------------
plot_gibbs_trace <- function(gs, true_sigma2 = NULL, true_tau2 = NULL) {
  
  par(mfrow = c(2, 2))
  
  # sigma2 trace
  plot(gs$sigma2_samples, type = "l", col = "steelblue",
       main = expression(paste("Trace: ", sigma^2)),
       ylab = expression(sigma^2), xlab = "iteration (post burn-in)")
  if (!is.null(true_sigma2)) abline(h = true_sigma2, col = "red", lty = 2, lwd = 2)
  
  # sigma2 density
  plot(density(gs$sigma2_samples), col = "steelblue", lwd = 2,
       main = expression(paste("Density: ", sigma^2)))
  if (!is.null(true_sigma2)) abline(v = true_sigma2, col = "red", lty = 2, lwd = 2)
  
  # tau2 trace
  plot(gs$tau2_samples, type = "l", col = "darkorange",
       main = expression(paste("Trace: ", tau^2)),
       ylab = expression(tau^2), xlab = "iteration (post burn-in)")
  if (!is.null(true_tau2)) abline(h = true_tau2, col = "red", lty = 2, lwd = 2)
  
  # tau2 density
  plot(density(gs$tau2_samples), col = "darkorange", lwd = 2,
       main = expression(paste("Density: ", tau^2)))
  if (!is.null(true_tau2)) abline(v = true_tau2, col = "red", lty = 2, lwd = 2)
}
