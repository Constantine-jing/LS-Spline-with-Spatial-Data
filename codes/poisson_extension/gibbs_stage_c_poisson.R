# ============================================================
# gibbs_stage_c_poisson.R
# Poisson spatial additive model — Gibbs/MH sampler (v2 FIXED)
#
# Extension of gibbs_stage_c_full.R to Poisson responses.
# Key difference: b can no longer be marginalized out, so we
# sample it explicitly.  The collapsed Gibbs is replaced by a
# full Gibbs with IWLS-based MH proposals (Gamerman 1997;
# Lang & Brezger 2004; Fahrmeir & Lang 2001).
#
# Model:
#   y_i ~ Poisson(lambda_i)
#   log(lambda_i) = [H eta]_i + b_i + offset_i
#
#   b ~ N(0, sigma2 * R(rho, nu))
#
# Priors:
#   mu ~ N(0, kappa2)                  [vague]
#   beta_j ~ N(0, tau2_s_j * K_j^-)   [RW2 smoothing]
#   tau2_s_j ~ IG(a_s, b_s)
#   sigma2 ~ IG(a_sigma, b_sigma)
#   rho    ~ logN(log_rho_mu, log_rho_sd)
#   nu     : FIXED
#
# Gibbs steps:
#   1) eta | y, b, {tau2_s_j}       via MH (IWLS proposal)
#   2) b   | y, eta, sigma2, rho    via MH (IWLS proposal)
#   3) sigma2 | b, rho              ~ IG (conjugate!)
#   4) rho    | b, sigma2           via MH
#   5) tau2_s_j | beta_j            ~ IG (conjugate)
#
# v2 FIXES:
#   - Smart initialization via Poisson GLM (was: all zeros)
#   - Proper reverse IWLS proposal for correct MH ratio
#   - Clamped log_lambda to prevent exp() overflow
#   - Working response uses partial residual correctly
#
# Depends on: spatial_utils.R, ls_basis.R
# ============================================================


# --- Penalty builders (same as Gaussian version) ---

build_rw2_penalty <- function(d) {
  if (d < 3) return(diag(d) * 0)
  D2 <- matrix(0, d - 2, d)
  for (i in 1:(d - 2)) {
    D2[i, i]     <-  1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <-  1
  }
  crossprod(D2)
}

build_block_prior_precision <- function(col_map, tau2_s, kappa2 = 1e6,
                                         eps_ridge = 1e-6) {
  p_covs <- length(col_map)
  p_total <- 1 + max(unlist(col_map))
  Q0 <- matrix(0, p_total, p_total)
  Q0[1, 1] <- 1 / kappa2
  for (j in 1:p_covs) {
    idx <- 1 + col_map[[j]]
    d_j <- length(idx)
    K_j <- build_rw2_penalty(d_j)
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s[j]) * K_j + eps_ridge * diag(d_j)
  }
  Q0
}


# ============================================================
# Poisson helpers
# ============================================================

poisson_loglik <- function(y, log_lambda) {
  log_lambda <- pmin(log_lambda, 20)
  sum(y * log_lambda - exp(log_lambda))
}

poisson_iwls <- function(y, log_lambda) {
  log_lambda <- pmin(log_lambda, 20)
  lambda <- exp(log_lambda)
  lambda <- pmax(lambda, 1e-10)
  w <- lambda
  z <- log_lambda + (y - lambda) / lambda
  list(z = z, w = w)
}


# ============================================================
# Smart initialization: Poisson GLM via IWLS
# ============================================================
poisson_glm_init <- function(y, H, offset, max_iter = 25, tol = 1e-6) {
  n <- length(y)
  p <- ncol(H)

  eta <- rep(0, p)
  eta[1] <- log(max(mean(y / exp(offset)), 0.1))
  log_lam <- as.vector(H %*% eta) + offset

  for (it in 1:max_iter) {
    iwls <- poisson_iwls(y, log_lam)
    # Working response for eta: subtract offset
    z <- iwls$z - offset
    w <- iwls$w

    HtWH <- crossprod(H, w * H) + 1e-6 * diag(p)
    rhs  <- crossprod(H, w * z)

    eta_new <- as.vector(solve(HtWH, rhs))
    log_lam_new <- as.vector(H %*% eta_new) + offset

    change <- max(abs(eta_new - eta))
    eta <- eta_new
    log_lam <- log_lam_new
    if (change < tol) break
  }

  lambda_hat <- exp(pmin(log_lam, 20))
  resid_log <- (y - lambda_hat) / pmax(lambda_hat, 0.1)

  list(eta = eta, log_lam = log_lam, residuals = resid_log)
}


# ============================================================
# gibbs_poisson_sampler()
# ============================================================
gibbs_poisson_sampler <- function(y, H, D, nu = 1.5, col_map,
                                   offset = NULL,
                                   n_iter  = 10000,
                                   n_burn  = 2500,
                                   n_thin  = 1,
                                   kappa2  = 1e6,
                                   a_sigma = 2, b_sigma = 1,
                                   a_smooth = 1, b_smooth = 0.005,
                                   log_rho_mu = -1.6,
                                   log_rho_sd = 1.0,
                                   iwls_scale_eta = 1.0,
                                   iwls_scale_b   = 1.0,
                                   mh_sd_log_rho  = 0.2,
                                   eps_ridge = 1e-6,
                                   jitter    = 1e-8,
                                   init      = NULL,
                                   verbose   = TRUE) {

  y <- as.integer(y)
  H <- as.matrix(H)
  D <- as.matrix(D)
  n <- length(y)
  p <- ncol(H)
  p_covs <- length(col_map)

  if (is.null(offset)) offset <- rep(0, n)
  offset <- as.numeric(offset)

  K_list <- vector("list", p_covs)
  d_list <- integer(p_covs)
  for (j in 1:p_covs) {
    d_j <- length(col_map[[j]])
    K_list[[j]] <- build_rw2_penalty(d_j)
    d_list[j] <- d_j
  }

  compute_R <- function(rho) { matern_cor(D, rho = rho, nu = nu) }

  log_lognormal_prior <- function(rho, mu, sd) {
    dnorm(log(rho), mean = mu, sd = sd, log = TRUE) - log(rho)
  }

  log_prior_b <- function(b, sigma2, L_R) {
    alpha <- forwardsolve(t(L_R), b)
    Rinv_b <- backsolve(L_R, alpha)
    logdet_R <- 2 * sum(log(diag(L_R)))
    -0.5 * n * log(sigma2) - 0.5 * logdet_R -
      0.5 / sigma2 * sum(b * Rinv_b)
  }

  # ============================================================
  # SMART INITIALIZATION
  # ============================================================
  if (is.null(init)) {
    cat("  Initializing via Poisson GLM...\n")
    glm_init <- poisson_glm_init(y, H, offset)
    eta    <- glm_init$eta
    b      <- glm_init$residuals * 0.1
    sigma2 <- max(var(b), 0.1)
    rho    <- 0.2
    tau2_s <- rep(1.0, p_covs)

    cat(sprintf("  GLM init: intercept=%.3f, sigma2=%.4f\n", eta[1], sigma2))
    cat(sprintf("  GLM init: log_lam range=[%.2f, %.2f]\n",
                min(glm_init$log_lam), max(glm_init$log_lam)))
  } else {
    sigma2 <- init$sigma2
    rho    <- if (!is.null(init$rho)) init$rho else 0.2
    eta    <- init$eta
    b      <- if (!is.null(init$b)) init$b else rep(0, n)
    tau2_s <- if (!is.null(init$tau2_s)) init$tau2_s else rep(1.0, p_covs)
  }

  R   <- compute_R(rho)
  L_R <- chol(R + diag(jitter, n))
  log_lam <- as.vector(H %*% eta) + b + offset

  cat(sprintf("  Initial loglik: %.1f\n", poisson_loglik(y, log_lam)))

  n_keep <- floor((n_iter - n_burn) / n_thin)
  eta_samples    <- matrix(NA, n_keep, p)
  b_samples      <- matrix(NA, n_keep, n)
  sigma2_samples <- numeric(n_keep)
  rho_samples    <- numeric(n_keep)
  tau2_s_samples <- matrix(NA, n_keep, p_covs)

  accept_eta <- 0; accept_b <- 0; accept_rho <- 0; keep_idx <- 0

  # ============================================================
  # MAIN LOOP
  # ============================================================
  for (iter in 1:n_iter) {

    # ----------------------------------------------------------
    # Step 1: MH for eta | y, b, {tau2_s_j}
    #
    # IWLS proposal at current eta:
    #   Working response: z = H*eta + (y - lam)/lam
    #   Working weight:   w = lam
    #   Proposal precision: P = H'WH + Q0
    #   Proposal mean: m = P^{-1} H'Wz
    #
    # For correct MH ratio with state-dependent proposal,
    # we also compute the REVERSE proposal (at eta_prop).
    # ----------------------------------------------------------
    Q0 <- build_block_prior_precision(col_map, tau2_s, kappa2, eps_ridge)
    Heta <- as.vector(H %*% eta)

    iwls <- poisson_iwls(y, log_lam)
    z_eta <- Heta + (y - iwls$w) / iwls$w   # H*eta + (y-lam)/lam
    W_diag <- iwls$w

    P_eta <- crossprod(H, W_diag * H) + Q0
    L_P_eta <- tryCatch(chol(P_eta), error = function(e) NULL)

    if (!is.null(L_P_eta)) {
      V_eta <- chol2inv(L_P_eta)
      m_fwd <- as.vector(V_eta %*% crossprod(H, W_diag * z_eta))

      eta_prop <- m_fwd + iwls_scale_eta * as.vector(backsolve(L_P_eta, rnorm(p)))
      log_lam_prop <- as.vector(H %*% eta_prop) + b + offset

      # Reverse proposal at eta_prop
      iwls_r <- poisson_iwls(y, log_lam_prop)
      Heta_prop <- as.vector(H %*% eta_prop)
      z_eta_r <- Heta_prop + (y - iwls_r$w) / iwls_r$w
      P_eta_r <- crossprod(H, iwls_r$w * H) + Q0
      L_P_eta_r <- tryCatch(chol(P_eta_r), error = function(e) NULL)

      if (!is.null(L_P_eta_r)) {
        m_rev <- as.vector(chol2inv(L_P_eta_r) %*% crossprod(H, iwls_r$w * z_eta_r))

        # Log-posterior
        lp_curr <- poisson_loglik(y, log_lam) -
          0.5 * as.numeric(t(eta) %*% Q0 %*% eta)
        lp_prop <- poisson_loglik(y, log_lam_prop) -
          0.5 * as.numeric(t(eta_prop) %*% Q0 %*% eta_prop)

        # Log-proposal densities (including normalizing constants)
        s2 <- iwls_scale_eta^2
        d_fwd <- eta_prop - m_fwd
        lq_fwd <- -0.5/s2 * as.numeric(t(d_fwd) %*% P_eta %*% d_fwd) +
          sum(log(diag(L_P_eta)))

        d_rev <- eta - m_rev
        lq_rev <- -0.5/s2 * as.numeric(t(d_rev) %*% P_eta_r %*% d_rev) +
          sum(log(diag(L_P_eta_r)))

        log_alpha <- (lp_prop - lp_curr) + (lq_rev - lq_fwd)

        if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
          eta <- eta_prop
          log_lam <- log_lam_prop
          accept_eta <- accept_eta + 1
        }
      }
    }

    # ----------------------------------------------------------
    # Step 2: MH for b | y, eta, sigma2, rho
    # ----------------------------------------------------------
    iwls <- poisson_iwls(y, log_lam)
    z_b <- b + (y - iwls$w) / iwls$w   # b + (y-lam)/lam

    Rinv <- chol2inv(L_R)
    P_b <- diag(iwls$w) + (1/sigma2) * Rinv

    L_P_b <- tryCatch(chol(P_b), error = function(e) NULL)
    if (!is.null(L_P_b)) {
      m_fwd_b <- as.vector(chol2inv(L_P_b) %*% (iwls$w * z_b))
      b_prop <- m_fwd_b + iwls_scale_b * as.vector(backsolve(L_P_b, rnorm(n)))

      log_lam_prop <- as.vector(H %*% eta) + b_prop + offset

      # Reverse proposal
      iwls_r_b <- poisson_iwls(y, log_lam_prop)
      z_b_r <- b_prop + (y - iwls_r_b$w) / iwls_r_b$w
      P_b_r <- diag(iwls_r_b$w) + (1/sigma2) * Rinv
      L_P_b_r <- tryCatch(chol(P_b_r), error = function(e) NULL)

      if (!is.null(L_P_b_r)) {
        m_rev_b <- as.vector(chol2inv(L_P_b_r) %*% (iwls_r_b$w * z_b_r))

        lp_curr_b <- poisson_loglik(y, log_lam) + log_prior_b(b, sigma2, L_R)
        lp_prop_b <- poisson_loglik(y, log_lam_prop) + log_prior_b(b_prop, sigma2, L_R)

        s2b <- iwls_scale_b^2
        d_fwd_b <- b_prop - m_fwd_b
        lq_fwd_b <- -0.5/s2b * as.numeric(t(d_fwd_b) %*% P_b %*% d_fwd_b) +
          sum(log(diag(L_P_b)))

        d_rev_b <- b - m_rev_b
        lq_rev_b <- -0.5/s2b * as.numeric(t(d_rev_b) %*% P_b_r %*% d_rev_b) +
          sum(log(diag(L_P_b_r)))

        log_alpha_b <- (lp_prop_b - lp_curr_b) + (lq_rev_b - lq_fwd_b)

        if (is.finite(log_alpha_b) && log(runif(1)) < log_alpha_b) {
          b <- b_prop
          log_lam <- log_lam_prop
          accept_b <- accept_b + 1
        }
      }
    }

    # ----------------------------------------------------------
    # Step 3: sigma2 | b, rho  (conjugate IG)
    # ----------------------------------------------------------
    Rinv_b <- as.vector(backsolve(L_R, forwardsolve(t(L_R), b)))
    qform_b <- sum(b * Rinv_b)
    sigma2 <- 1/rgamma(1, shape = a_sigma + n/2, rate = b_sigma + qform_b/2)

    # ----------------------------------------------------------
    # Step 4: MH for rho | b, sigma2
    # ----------------------------------------------------------
    log_rho_prop <- log(rho) + rnorm(1, 0, mh_sd_log_rho)
    rho_prop <- exp(log_rho_prop)

    R_prop   <- compute_R(rho_prop)
    L_R_prop <- tryCatch(chol(R_prop + diag(jitter, n)), error = function(e) NULL)

    if (!is.null(L_R_prop)) {
      lp_rho_curr <- log_prior_b(b, sigma2, L_R) +
        log_lognormal_prior(rho, log_rho_mu, log_rho_sd)
      lp_rho_prop <- log_prior_b(b, sigma2, L_R_prop) +
        log_lognormal_prior(rho_prop, log_rho_mu, log_rho_sd)

      if (is.finite(lp_rho_prop - lp_rho_curr) &&
          log(runif(1)) < (lp_rho_prop - lp_rho_curr)) {
        rho <- rho_prop; R <- R_prop; L_R <- L_R_prop
        accept_rho <- accept_rho + 1
      }
    }

    # ----------------------------------------------------------
    # Step 5: tau2_s_j | beta_j  (conjugate IG)
    # ----------------------------------------------------------
    for (j in 1:p_covs) {
      idx <- 1 + col_map[[j]]
      beta_j <- eta[idx]
      qform <- as.numeric(t(beta_j) %*% K_list[[j]] %*% beta_j)
      rank_Kj <- max(d_list[j] - 2, 1)
      tau2_s[j] <- 1/rgamma(1, shape = a_smooth + rank_Kj/2,
                              rate = b_smooth + qform/2)
    }

    # ----------------------------------------------------------
    # Store
    # ----------------------------------------------------------
    if (iter > n_burn && ((iter - n_burn) %% n_thin == 0)) {
      keep_idx <- keep_idx + 1
      eta_samples[keep_idx, ]    <- eta
      b_samples[keep_idx, ]      <- b
      sigma2_samples[keep_idx]   <- sigma2
      rho_samples[keep_idx]      <- rho
      tau2_s_samples[keep_idx, ] <- tau2_s
    }

    if (verbose && (iter %% 500 == 0)) {
      cat(sprintf("  iter %d/%d  rho=%.4f  sig2=%.4f  ll=%.1f  acc(eta=%.2f,b=%.2f,rho=%.2f)\n",
                  iter, n_iter, rho, sigma2,
                  poisson_loglik(y, log_lam),
                  accept_eta/iter, accept_b/iter, accept_rho/iter))
    }
  }

  cat(sprintf("  MH acceptance: eta=%.3f  b=%.3f  rho=%.3f\n",
              accept_eta/n_iter, accept_b/n_iter, accept_rho/n_iter))

  list(eta_samples = eta_samples, b_samples = b_samples,
       sigma2_samples = sigma2_samples, rho_samples = rho_samples,
       tau2_s_samples = tau2_s_samples,
       n_burn = n_burn, n_thin = n_thin, n_iter = n_iter,
       accept_rate = c(eta = accept_eta/n_iter, b = accept_b/n_iter,
                       rho = accept_rho/n_iter))
}


# ============================================================
# plot_gibbs_trace_poisson()
# ============================================================
plot_gibbs_trace_poisson <- function(gs, true_sigma2 = NULL, true_rho = NULL) {
  par(mfrow = c(2, 2))
  plot(gs$sigma2_samples, type="l", col="steelblue",
       main=expression(paste("Trace: ",sigma^2)), ylab=expression(sigma^2), xlab="iteration")
  if (!is.null(true_sigma2)) abline(h=true_sigma2, col="red", lty=2, lwd=2)
  plot(density(gs$sigma2_samples), col="steelblue", lwd=2,
       main=expression(paste("Density: ",sigma^2)))
  if (!is.null(true_sigma2)) abline(v=true_sigma2, col="red", lty=2, lwd=2)
  plot(gs$rho_samples, type="l", col="forestgreen",
       main=expression(paste("Trace: ",rho)), ylab=expression(rho), xlab="iteration")
  if (!is.null(true_rho)) abline(h=true_rho, col="red", lty=2, lwd=2)
  plot(density(gs$rho_samples), col="forestgreen", lwd=2,
       main=expression(paste("Density: ",rho)))
  if (!is.null(true_rho)) abline(v=true_rho, col="red", lty=2, lwd=2)
}
