# ============================================================
# gibbs_bayes_v2.R
# Stage B (v2): Gibbs sampler with reparameterized variance
#
# FIX: The v1 sampler had sigma2 explode because sigma2 and
# tau2 were sampled independently, allowing a runaway where
# b absorbs everything. 
#
# Solution: parameterize as
#   sigma2_total = sigma2 + tau2  (total marginal variance)
#   lambda = tau2 / sigma2        (noise-to-signal ratio)
#
# But this breaks conjugacy for lambda. Simpler practical fix:
# sample (sigma2 | b, tau2) with a proper IG prior whose
# scale is informed by data, and clamp sigma2 to a sensible
# range.
#
# ACTUALLY the cleanest fix: keep the original parameterization
# but use a prior on sigma (the SD, not the variance) and
# integrate b out when sampling variance components.
#
# SIMPLEST CORRECT FIX (what we implement):
#   - Do NOT sample b explicitly
#   - Marginalize b out: y ~ N(H eta, sigma2*R + tau2*I)
#   - Sample eta from this marginal (= your Baby Bayes step)
#   - Sample sigma2, tau2 via Metropolis-Hastings on the
#     marginal likelihood (profiling b out)
#   - THEN compute b posterior mean as a derived quantity
#
# This is the "marginalized" or "collapsed" Gibbs approach.
# It avoids the b-sigma2 coupling entirely.
# ============================================================

# ------------------------------------------------------------
# gibbs_collapsed_sampler()
#
# Collapsed Gibbs: marginalize over b.
#
# y | eta, sigma2, tau2 ~ N(H eta, sigma2*R + tau2*I)
#
# Steps:
#   1) eta | y, sigma2, tau2    ~ MVN  (same as Baby Bayes)
#   2) sigma2 | y, eta, tau2    via MH on marginal likelihood
#   3) tau2 | y, eta, sigma2    via MH on marginal likelihood
#   4) b (derived): E[b|y,eta,sigma2,tau2] = sigma2*R*Sigma^{-1}(y-H*eta)
#
# For steps 2-3, we use random walk MH on log(sigma2) and
# log(tau2) with the full marginal log-likelihood.
# ------------------------------------------------------------
gibbs_collapsed_sampler <- function(y, H, R,
                                     n_iter  = 5000,
                                     n_burn  = 1000,
                                     n_thin  = 1,
                                     kappa2  = 1e6,
                                     # IG priors on sigma2 and tau2
                                     a_sigma = 2, b_sigma = 1,
                                     a_tau   = 2, b_tau   = 0.3,
                                     # MH tuning
                                     mh_sd_log_sigma2 = 0.3,
                                     mh_sd_log_tau2   = 0.3,
                                     init    = NULL,
                                     jitter  = 1e-8,
                                     verbose = TRUE) {
  
  y <- as.numeric(y)
  H <- as.matrix(H)
  R <- as.matrix(R)
  n <- length(y)
  p <- ncol(H)
  
  Lambda0_inv <- diag(1 / kappa2, p)
  
  # --- Helper: log marginal likelihood for (sigma2, tau2) given eta ---
  # y - H*eta ~ N(0, sigma2*R + tau2*I)
  log_marg_lik <- function(resid, sigma2, tau2) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- tryCatch(chol(Sigma + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L)) return(-Inf)
    
    logdet <- 2 * sum(log(diag(L)))
    alpha  <- forwardsolve(t(L), resid)
    qform  <- sum(alpha^2)
    
    -0.5 * (logdet + qform)
  }
  
  # --- Helper: log IG prior ---
  log_ig_prior <- function(x, a, b) {
    # x ~ IG(a, b) means p(x) propto x^{-(a+1)} exp(-b/x)
    -(a + 1) * log(x) - b / x
  }
  
  # --- Helper: build Sigma and sample eta ---
  sample_eta <- function(y, sigma2, tau2) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- chol(Sigma + diag(jitter, n))
    y_w <- forwardsolve(t(L), y)
    H_w <- forwardsolve(t(L), H)
    
    Q_eta <- crossprod(H_w) + Lambda0_inv
    U_eta <- chol(Q_eta)
    V_eta <- chol2inv(U_eta)
    m_eta <- as.vector(V_eta %*% crossprod(H_w, y_w))
    
    z <- rnorm(p)
    eta <- m_eta + as.vector(backsolve(U_eta, z))
    list(eta = eta, m_eta = m_eta)
  }
  
  # --- Helper: compute b posterior mean (derived quantity) ---
  compute_b_postmean <- function(resid, sigma2, tau2) {
    Sigma <- sigma2 * R + tau2 * diag(n)
    L <- chol(Sigma + diag(jitter, n))
    alpha <- forwardsolve(t(L), forwardsolve(t(L), resid, upper.tri = FALSE))
    # Wait — need Sigma^{-1} resid properly
    alpha <- chol2inv(L) %*% resid
    as.vector(sigma2 * R %*% alpha)
  }
  
  # --- Initialize ---
  if (is.null(init)) {
    sigma2 <- 1.0
    tau2   <- 0.5
    eta    <- rep(0, p)
  } else {
    sigma2 <- init$sigma2
    tau2   <- init$tau2
    eta    <- init$eta
  }
  
  # --- Storage ---
  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples    <- matrix(NA, n_keep, p)
  b_samples      <- matrix(NA, n_keep, n)
  sigma2_samples <- numeric(n_keep)
  tau2_samples   <- numeric(n_keep)
  
  accept_sigma2 <- 0
  accept_tau2   <- 0
  keep_idx      <- 0
  
  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in 1:n_iter) {
    
    # ----------------------------------------------------------
    # Step 1: Draw eta | y, sigma2, tau2  (collapsed over b)
    # ----------------------------------------------------------
    res <- sample_eta(y, sigma2, tau2)
    eta <- res$eta
    
    # Residual for variance component updates
    resid <- as.vector(y - H %*% eta)
    
    # ----------------------------------------------------------
    # Step 2: MH for sigma2 | y, eta, tau2
    #
    # Propose log(sigma2*) = log(sigma2) + N(0, mh_sd^2)
    # Accept/reject based on marginal likelihood + prior
    # ----------------------------------------------------------
    log_sigma2_prop <- log(sigma2) + rnorm(1, 0, mh_sd_log_sigma2)
    sigma2_prop <- exp(log_sigma2_prop)
    
    log_target_curr <- log_marg_lik(resid, sigma2, tau2) + 
                       log_ig_prior(sigma2, a_sigma, b_sigma)
    log_target_prop <- log_marg_lik(resid, sigma2_prop, tau2) + 
                       log_ig_prior(sigma2_prop, a_sigma, b_sigma)
    
    # Jacobian for log-scale proposal: cancels (symmetric on log scale)
    log_alpha <- log_target_prop - log_target_curr
    
    if (log(runif(1)) < log_alpha) {
      sigma2 <- sigma2_prop
      accept_sigma2 <- accept_sigma2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 3: MH for tau2 | y, eta, sigma2
    # ----------------------------------------------------------
    log_tau2_prop <- log(tau2) + rnorm(1, 0, mh_sd_log_tau2)
    tau2_prop <- exp(log_tau2_prop)
    
    log_target_curr <- log_marg_lik(resid, sigma2, tau2) + 
                       log_ig_prior(tau2, a_tau, b_tau)
    log_target_prop <- log_marg_lik(resid, sigma2, tau2_prop) + 
                       log_ig_prior(tau2_prop, a_tau, b_tau)
    
    log_alpha <- log_target_prop - log_target_curr
    
    if (log(runif(1)) < log_alpha) {
      tau2 <- tau2_prop
      accept_tau2 <- accept_tau2 + 1
    }
    
    # ----------------------------------------------------------
    # Step 4: Compute b posterior mean (derived)
    # ----------------------------------------------------------
    b_mean <- compute_b_postmean(resid, sigma2, tau2)
    
    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]  <- eta
      b_samples[keep_idx, ]    <- b_mean
      sigma2_samples[keep_idx] <- sigma2
      tau2_samples[keep_idx]   <- tau2
    }
    
    if (verbose && (iter %% 1000 == 0)) {
      cat(sprintf("  iter %d/%d  sigma2=%.4f  tau2=%.4f  acc_s=%.2f  acc_t=%.2f\n",
                  iter, n_iter, sigma2, tau2,
                  accept_sigma2 / iter, accept_tau2 / iter))
    }
  }
  
  cat(sprintf("  Final MH acceptance: sigma2=%.3f  tau2=%.3f\n",
              accept_sigma2 / n_iter, accept_tau2 / n_iter))
  
  list(
    eta_samples    = eta_samples,
    b_samples      = b_samples,
    sigma2_samples = sigma2_samples,
    tau2_samples   = tau2_samples,
    n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
    accept_rate = c(sigma2 = accept_sigma2 / n_iter,
                    tau2   = accept_tau2 / n_iter)
  )
}
