# ============================================================
# gibbs_stage_c.R
# Stage C: Collapsed Gibbs + RW2 smoothing prior on beta
#
# Key changes from Stage B (v2):
#   1) Use MORE KNOTS (M=20 instead of 6) for flexibility
#   2) Replace vague N(0, kappa2*I) prior on beta_j with a
#      second-order random walk (RW2) prior:
#        beta_j ~ N(0, tau2_s_j * K_j^{-})
#      where K_j is the second-difference penalty matrix
#      and tau2_s_j is a per-covariate smoothing variance.
#   3) Sample tau2_s_j from its IG full conditional.
#
# The RW2 prior encodes: "adjacent spline coefficients should
# change smoothly" — the Bayesian analog of the P-spline penalty.
#
# Model:
#   y = H eta + b + eps
#   eta = (mu, beta_1, ..., beta_p)
#   b ~ N(0, sigma2 * R)       [marginalized out]
#   eps ~ N(0, tau2 * I)
#
# Priors:
#   mu ~ N(0, kappa2)                      [vague, intercept]
#   beta_j ~ N(0, tau2_s_j * K_j^{-})      [RW2 smoothing]
#   tau2_s_j ~ IG(a_s, b_s)                [smoothing variance]
#   sigma2 ~ IG(a_sigma, b_sigma)
#   tau2 ~ IG(a_tau, b_tau)
#
# Depends on: spatial_utils.R, ls_basis.R, gibbs_bayes_v2.R
# ============================================================


# ------------------------------------------------------------
# build_rw2_penalty(d)
#
# Build second-order difference penalty matrix K for
# a vector of d coefficients.
#
# RW2: beta_k = 2*beta_{k-1} - beta_{k-2} + u_k
# This gives K = D2' D2 where D2 is the (d-2) x d
# second-difference matrix.
#
# K is rank (d-2), so it's singular. That's fine —
# the RW2 prior is improper, but becomes proper when
# combined with the likelihood.
# ------------------------------------------------------------
build_rw2_penalty <- function(d) {
  if (d < 3) return(diag(d) * 0)  # no penalty for <3 coefs
  
  D2 <- matrix(0, d - 2, d)
  for (i in 1:(d - 2)) {
    D2[i, i]     <-  1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <-  1
  }
  crossprod(D2)  # K = D2' D2, symmetric positive semi-definite
}


# ------------------------------------------------------------
# build_block_prior_precision()
#
# Build the full prior precision matrix for eta = (mu, beta_1, ..., beta_p).
#
# - mu gets precision 1/kappa2  (vague)
# - beta_j gets precision (1/tau2_s_j) * K_j  (RW2)
#
# Since K_j is rank-deficient (rank d_j - 2), we add a tiny
# ridge (eps_ridge * I) to make it invertible.
# This is equivalent to a very weak proper prior on the
# null space of K_j (constant + linear).
# ------------------------------------------------------------
build_block_prior_precision <- function(col_map, tau2_s, kappa2 = 1e6,
                                         eps_ridge = 1e-6) {
  p_covs <- length(col_map)
  p_total <- 1 + max(unlist(col_map))  # intercept + all beta cols
  
  Q0 <- matrix(0, p_total, p_total)
  
  # Intercept precision
  Q0[1, 1] <- 1 / kappa2
  
  # Per-covariate RW2 penalty
  for (j in 1:p_covs) {
    idx <- 1 + col_map[[j]]  # +1 for intercept
    d_j <- length(idx)
    K_j <- build_rw2_penalty(d_j)
    
    # Prior precision: (1/tau2_s_j) * K_j + eps_ridge * I
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s[j]) * K_j + eps_ridge * diag(d_j)
  }
  
  Q0
}


# ------------------------------------------------------------
# gibbs_stage_c_sampler()
#
# Collapsed Gibbs with RW2 smoothing prior.
#
# Steps per iteration:
#   1) eta | y, sigma2, tau2, {tau2_s_j}  ~ MVN (collapsed over b)
#   2) sigma2 | y, eta, tau2              via MH
#   3) tau2 | y, eta, sigma2              via MH
#   4) tau2_s_j | beta_j                  ~ IG (conjugate)
#   5) b (derived posterior mean)
# ------------------------------------------------------------
gibbs_stage_c_sampler <- function(y, H, R, col_map,
                                   n_iter  = 5000,
                                   n_burn  = 1000,
                                   n_thin  = 1,
                                   kappa2  = 1e6,
                                   # IG priors on sigma2, tau2
                                   a_sigma = 2, b_sigma = 1,
                                   a_tau   = 2, b_tau   = 0.3,
                                   # IG priors on smoothing variances
                                   a_smooth = 1, b_smooth = 0.005,
                                   # MH tuning
                                   mh_sd_log_sigma2 = 0.3,
                                   mh_sd_log_tau2   = 0.3,
                                   eps_ridge = 1e-6,
                                   init    = NULL,
                                   jitter  = 1e-8,
                                   verbose = TRUE) {
  
  y <- as.numeric(y)
  H <- as.matrix(H)
  R <- as.matrix(R)
  n <- length(y)
  p <- ncol(H)
  p_covs <- length(col_map)
  
  # Pre-compute K_j matrices
  K_list <- vector("list", p_covs)
  d_list <- integer(p_covs)
  for (j in 1:p_covs) {
    d_j <- length(col_map[[j]])
    K_list[[j]] <- build_rw2_penalty(d_j)
    d_list[j] <- d_j
  }
  
  # --- Helpers ---
  log_marg_lik <- function(resid, sigma2, tau2) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- tryCatch(chol(Sigma + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L)) return(-Inf)
    logdet <- 2 * sum(log(diag(L)))
    alpha <- forwardsolve(t(L), resid)
    -0.5 * (logdet + sum(alpha^2))
  }
  
  log_ig_prior <- function(x, a, b) {
    -(a + 1) * log(x) - b / x
  }
  
  compute_b_postmean <- function(resid, sigma2, tau2) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    Sigma_inv <- chol2inv(chol(Sigma + diag(jitter, n)))
    as.vector(sigma2 * R %*% Sigma_inv %*% resid)
  }
  
  # --- Initialize ---
  if (is.null(init)) {
    sigma2 <- 1.0
    tau2   <- 0.5
    eta    <- rep(0, p)
    tau2_s <- rep(1.0, p_covs)
  } else {
    sigma2 <- init$sigma2
    tau2   <- init$tau2
    eta    <- init$eta
    tau2_s <- if (!is.null(init$tau2_s)) init$tau2_s else rep(1.0, p_covs)
  }
  
  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples    <- matrix(NA, n_keep, p)
  b_samples      <- matrix(NA, n_keep, n)
  sigma2_samples <- numeric(n_keep)
  tau2_samples   <- numeric(n_keep)
  tau2_s_samples <- matrix(NA, n_keep, p_covs)
  
  accept_sigma2 <- 0
  accept_tau2   <- 0
  keep_idx      <- 0
  
  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in 1:n_iter) {
    
    # ----------------------------------------------------------
    # Step 1: Draw eta | y, sigma2, tau2, {tau2_s_j}
    #
    # Marginal likelihood: y ~ N(H eta, sigma2*R + tau2*I)
    # Prior on eta: precision = Q0 (block diagonal with RW2)
    #
    # Posterior precision: Q_eta = H' Sigma^{-1} H + Q0
    # Posterior mean:      m_eta = Q_eta^{-1} H' Sigma^{-1} y
    # ----------------------------------------------------------
    Q0 <- build_block_prior_precision(col_map, tau2_s, kappa2, eps_ridge)
    
    Sigma <- sigma2 * R + tau2 * diag(n)
    L_Sigma <- chol(Sigma + diag(jitter, n))
    y_w <- forwardsolve(t(L_Sigma), y)
    H_w <- forwardsolve(t(L_Sigma), H)
    
    Q_eta <- crossprod(H_w) + Q0
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% crossprod(H_w, y_w))
    
    z <- rnorm(p)
    eta <- m_eta + as.vector(backsolve(U_eta, z))
    
    resid <- as.vector(y - H %*% eta)
    
    # ----------------------------------------------------------
    # Step 2: MH for sigma2
    # ----------------------------------------------------------
    log_sigma2_prop <- log(sigma2) + rnorm(1, 0, mh_sd_log_sigma2)
    sigma2_prop <- exp(log_sigma2_prop)
    
    lp_curr <- log_marg_lik(resid, sigma2, tau2) + log_ig_prior(sigma2, a_sigma, b_sigma)
    lp_prop <- log_marg_lik(resid, sigma2_prop, tau2) + log_ig_prior(sigma2_prop, a_sigma, b_sigma)
    
    if (log(runif(1)) < (lp_prop - lp_curr)) {
      sigma2 <- sigma2_prop
      accept_sigma2 <- accept_sigma2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 3: MH for tau2
    # ----------------------------------------------------------
    log_tau2_prop <- log(tau2) + rnorm(1, 0, mh_sd_log_tau2)
    tau2_prop <- exp(log_tau2_prop)
    
    lp_curr <- log_marg_lik(resid, sigma2, tau2) + log_ig_prior(tau2, a_tau, b_tau)
    lp_prop <- log_marg_lik(resid, sigma2, tau2_prop) + log_ig_prior(tau2_prop, a_tau, b_tau)
    
    if (log(runif(1)) < (lp_prop - lp_curr)) {
      tau2 <- tau2_prop
      accept_tau2 <- accept_tau2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 4: Draw tau2_s_j | beta_j  (conjugate IG)
    #
    # Prior: tau2_s_j ~ IG(a_smooth, b_smooth)
    # Likelihood: beta_j ~ N(0, tau2_s_j * K_j^{-})
    #   => beta_j' K_j beta_j / tau2_s_j
    #
    # Posterior: tau2_s_j ~ IG(a_smooth + rank(K_j)/2,
    #                          b_smooth + beta_j' K_j beta_j / 2)
    # ----------------------------------------------------------
    for (j in 1:p_covs) {
      idx <- 1 + col_map[[j]]
      beta_j <- eta[idx]
      K_j <- K_list[[j]]
      
      qform <- as.numeric(t(beta_j) %*% K_j %*% beta_j)
      rank_Kj <- d_list[j] - 2  # rank of second-difference penalty
      if (rank_Kj < 1) rank_Kj <- 1
      
      a_post <- a_smooth + rank_Kj / 2
      b_post <- b_smooth + qform / 2
      
      tau2_s[j] <- 1 / rgamma(1, shape = a_post, rate = b_post)
    }
    
    # ----------------------------------------------------------
    # Step 5: Compute b posterior mean (derived)
    # ----------------------------------------------------------
    b_mean <- compute_b_postmean(resid, sigma2, tau2)
    
    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]    <- eta
      b_samples[keep_idx, ]      <- b_mean
      sigma2_samples[keep_idx]   <- sigma2
      tau2_samples[keep_idx]     <- tau2
      tau2_s_samples[keep_idx, ] <- tau2_s
    }
    
    if (verbose && (iter %% 1000 == 0)) {
      cat(sprintf("  iter %d/%d  sigma2=%.4f  tau2=%.4f  tau2_s=%s\n",
                  iter, n_iter, sigma2, tau2,
                  paste(sprintf("%.3f", tau2_s), collapse=",")))
    }
  }
  
  cat(sprintf("  MH acceptance: sigma2=%.3f  tau2=%.3f\n",
              accept_sigma2 / n_iter, accept_tau2 / n_iter))
  
  list(
    eta_samples    = eta_samples,
    b_samples      = b_samples,
    sigma2_samples = sigma2_samples,
    tau2_samples   = tau2_samples,
    tau2_s_samples = tau2_s_samples,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
    accept_rate = c(sigma2 = accept_sigma2 / n_iter,
                    tau2   = accept_tau2 / n_iter)
  )
}
