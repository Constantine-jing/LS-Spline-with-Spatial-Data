# ============================================================
# fit_mgcv_wrapper.R
# 
# Wraps mgcv::gam to produce the canonical fit_result schema
# for the Scenario A cross-method comparison.
# 
# Model:
#   y = mu + sum_j s_j(X_j; bs="ps", k=20) + s_spatial(lon, lat; bs="gp", m=c(2,...)) + eps
# 
# Smooth bases (locked for fairness):
#   - Main effects:  bs = "ps"  (P-splines with 2nd-order RW penalty)
#                    k  = 20    (matches our M=20)
#                    m  = c(2, 2)  (cubic B-spline, 2nd-order penalty)
#   - Spatial:       bs = "gp"  (Gaussian process / Matern)
#                    m  = c(2, ...)   (m[1] = 2 -> nu = 1.0 in 2D)
#                    rho estimated by REML (matches our sampler's
#                    sampling rho from the posterior)
#                    nu = 1.0 fixed (matches Scenario A truth)
# 
# Joint posterior draws of each smooth on the eval grid are
# constructed via MASS::mvrnorm(n_draws, mu_g, V_g) where
#   mu_g = X_pred %*% coef(gam)
#   V_g  = X_pred %*% gam$Vp %*% t(X_pred)
# X_pred is the lp matrix from predict(..., type="lpmatrix")
# restricted to the columns of the smooth in question.
# 
# This gives statistically honest joint posterior draws under
# mgcv's empirical Bayes view (Wood 2006). For our pointwise
# RMSE / pointwise 95%-coverage metrics the result is
# numerically identical to mgcv's analytic intervals; joint
# structure costs nothing extra.
# ============================================================

source("canonical_schema.R")

suppressMessages({
  library(mgcv)
  library(MASS)
})


fit_mgcv <- function(sim, settings = NULL) {
  
  # ---- locked-in Stage 1 settings ----
  defaults <- list(
    k_main      = 20L,        # basis count per main-effect smooth
    bs_main     = "ps",       # P-spline (RW2-penalized B-spline)
    m_main      = c(2L, 2L),  # cubic spline, 2nd-order difference penalty
    
    bs_spatial  = "gp",       # Gaussian process (Matern)
    m_spatial   = c(2),       # m[1]=2 -> Matern nu=1.0 in 2D
                              # NOTE: leaving m as length 1 lets mgcv
                              # estimate rho by REML, matching how our
                              # Gibbs samples rho.
    k_spatial   = 100L,       # bigger basis for spatial term (mgcv default ok)
    
    method      = "REML",     # standard for mgcv
    n_draws     = 2000L       # MUST match the locked-in canonical n_draws
  )
  if (!is.null(settings)) {
    for (nm in names(settings)) defaults[[nm]] <- settings[[nm]]
  }
  S <- defaults
  
  if (abs(1.0 - sim$truth_params$nu) > 1e-12) {
    warning(sprintf(
      "fit_mgcv hardcodes nu=1 (m_spatial[1]=2) but sim$truth_params$nu=%g; result is biased",
      sim$truth_params$nu))
  }
  
  t0 <- proc.time()
  
  # ---- assemble the data frame ----
  df <- sim$data
  smooth_names <- paste0("X", seq_len(sim$p))
  stopifnot(all(smooth_names %in% colnames(df)))
  stopifnot(all(c("lon", "lat", "y") %in% colnames(df)))
  
  # ---- build formula programmatically ----
  smooth_terms <- sapply(smooth_names, function(nm) {
    sprintf("s(%s, bs='%s', k=%d, m=c(%d,%d))",
            nm, S$bs_main, S$k_main, S$m_main[1], S$m_main[2])
  })
  spatial_term <- sprintf("s(lon, lat, bs='%s', k=%d, m=c(%d))",
                          S$bs_spatial, S$k_spatial, S$m_spatial[1])
  rhs <- paste(c(smooth_terms, spatial_term), collapse = " + ")
  fmla <- as.formula(paste("y ~", rhs))
  
  # ---- fit ----
  t_fit_start <- proc.time()
  gam_fit <- gam(fmla, data = df, method = S$method)
  fit_sec <- (proc.time() - t_fit_start)[3]
  
  # ---- joint posterior draws from N(beta_hat, Vp) ----
  # We sample once for the FULL coefficient vector, then push the
  # same draws through each smooth's lpmatrix. This preserves
  # cross-smooth correlation if anyone ever cares; and it's the
  # right thing to do.
  t_post_start <- proc.time()
  
  beta_hat <- coef(gam_fit)            # full coef vector
  Vp       <- gam_fit$Vp               # full posterior covariance
  
  set.seed(sim$seed + 1000L)           # reproducible mgcv draws
  beta_draws <- MASS::mvrnorm(n = S$n_draws, mu = beta_hat, Sigma = Vp)
  # beta_draws: (n_draws x p_total)
  
  # ---- main-effect smooths on x_grid_1d ----
  # Trick: build a synthetic newdata where one column varies
  # over the eval grid and all other columns are at their training
  # means. Then predict(type="lpmatrix") gives the lp matrix for
  # the WHOLE model. We zero out columns that don't belong to the
  # j-th main-effect smooth so f_j(grid) = X_pred[, cols_j] %*% beta[cols_j].
  # 
  # Use the smooth-block label parsing from gam_fit$smooth so we
  # don't have to count columns by hand.
  smooth_blocks <- vector("list", length(gam_fit$smooth))
  names(smooth_blocks) <- sapply(gam_fit$smooth, function(s) s$label)
  for (k in seq_along(gam_fit$smooth)) {
    sm <- gam_fit$smooth[[k]]
    smooth_blocks[[k]] <- sm$first.para:sm$last.para
  }
  
  n_grid_1d <- length(sim$x_grid_1d)
  
  # build a mean-of-training newdata template
  template <- df[1, , drop = FALSE]
  for (nm in colnames(df)) {
    if (is.numeric(df[[nm]])) template[[nm]] <- mean(df[[nm]])
  }
  template <- template[, , drop = FALSE]
  
  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    nm_j <- smooth_names[j]
    nd <- template[rep(1, n_grid_1d), , drop = FALSE]
    nd[[nm_j]] <- sim$x_grid_1d
    Xp <- predict(gam_fit, newdata = nd, type = "lpmatrix")
    
    # find the smooth-block columns for s(X_j)
    label_j <- sprintf("s(%s)", nm_j)
    if (!label_j %in% names(smooth_blocks)) {
      stop(sprintf("Could not find smooth block '%s' in gam_fit", label_j))
    }
    cols_j <- smooth_blocks[[label_j]]
    
    # f_j on grid for each draw: (n_grid_1d x n_draws)
    Xp_j <- Xp[, cols_j, drop = FALSE]                       # (n_grid x |cols|)
    bj   <- beta_draws[, cols_j, drop = FALSE]               # (n_draws x |cols|)
    f_main[[j]] <- Xp_j %*% t(bj)                            # (n_grid x n_draws)
    
    # IDENTIFIABILITY: mgcv's smooths are constructed sum-to-zero
    # by absorbing identifiability constraints, so f_j is already
    # roughly mean-zero on the data points.  We re-center on the
    # eval grid for a like-for-like comparison with our wrapper.
    col_means <- colMeans(f_main[[j]])
    f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
  }
  
  # ---- spatial random effect at observed locations ----
  # Same pattern: build lpmatrix at the observed (lon, lat),
  # restrict to the spatial smooth's columns.
  Xp_obs <- predict(gam_fit, newdata = df, type = "lpmatrix")
  label_s <- "s(lon,lat)"
  if (!label_s %in% names(smooth_blocks)) {
    stop("Could not find spatial smooth block 's(lon,lat)'")
  }
  cols_s <- smooth_blocks[[label_s]]
  Xp_s   <- Xp_obs[, cols_s, drop = FALSE]
  bs_d   <- beta_draws[, cols_s, drop = FALSE]
  s_obs  <- Xp_s %*% t(bs_d)                                  # (n x n_draws)
  
  # ---- variance components ----
  # mgcv parameterizes smooths via smoothing parameters lambda_j;
  # the implied Bayesian variance is sigma^2 / lambda_j (where
  # sigma^2 is the residual variance scale parameter).
  # We report:
  #   sigma2: gam_fit$sig2 plus a draws-shaped vector of
  #           Vp-implied draws if available (we use a constant vec
  #           because mgcv treats sig2 as a point estimate).
  #   tau2_j: sigma2 / lambda_j  for each main-effect smooth.
  #   tau2_s: sigma2 / lambda_spatial
  # 
  # mgcv does NOT give a draws-based posterior of these; we fill
  # length-n_draws vectors with the REML point estimate. The
  # canonical schema then computes "RMSE of point estimates" with
  # zero variance, which matches mgcv's frequentist nature.
  sig2 <- as.numeric(gam_fit$sig2)
  
  sp_named <- gam_fit$sp                  # named numeric vector of lambdas
  # match smooth labels to sp names. mgcv names sp like "s(X1)", "s(lon,lat)"
  get_lambda <- function(label) {
    if (label %in% names(sp_named)) sp_named[[label]]
    else NA_real_
  }
  
  tau2_j_pt <- numeric(sim$p)
  for (j in seq_len(sim$p)) {
    lam <- get_lambda(sprintf("s(%s)", smooth_names[j]))
    tau2_j_pt[j] <- sig2 / lam
  }
  lam_s <- get_lambda("s(lon,lat)")
  tau2_s_pt <- sig2 / lam_s
  
  rep_n <- function(x) rep(x, S$n_draws)
  var_comp <- list(
    sigma2   = rep_n(sig2),
    tau2     = lapply(tau2_j_pt, rep_n),
    tau2_s   = rep_n(tau2_s_pt),
    tau2_int = list(),
    rho      = NULL,    # mgcv's gp basis estimates rho but it is
                        # NOT a sampled posterior; reported via
                        # settings$rho_hat below as a point estimate
    nu       = NULL     # fixed via m_spatial[1]=2
  )
  
  # extract mgcv's point-estimated rho for the spatial smooth, if
  # we can. mgcv stores it in gam_fit$smooth[[k]]$gp.smooth.spec or
  # gam_fit$smooth[[k]]$xt$rho depending on version. Best-effort.
  rho_hat <- NA_real_
  for (sm in gam_fit$smooth) {
    if (identical(sm$label, "s(lon,lat)")) {
      cand <- c(sm$gp.theta, sm$xt$rho, sm$rho)
      cand <- cand[!is.null(cand)]
      if (length(cand) >= 1) rho_hat <- as.numeric(cand[1])
      break
    }
  }
  
  post_sec <- (proc.time() - t_post_start)[3]
  total_sec <- (proc.time() - t0)[3]
  
  # ---- assemble ----
  fit <- list(
    method       = "mgcv",
    scenario     = sim$scenario,
    seed         = sim$seed,
    f_main       = f_main,
    f_int        = list(),
    s_obs        = s_obs,
    var_comp     = var_comp,
    x_grid_1d    = sim$x_grid_1d,
    x_grid_2d    = sim$x_grid_2d,
    timing       = list(total_sec = as.numeric(total_sec),
                        fit_sec   = as.numeric(fit_sec),
                        post_sec  = as.numeric(post_sec)),
    convergence  = list(
      rhat = NULL,
      ess  = NULL,
      reml_score = as.numeric(gam_fit$gcv.ubre),  # REML score
      converged  = isTRUE(gam_fit$converged)
    ),
    settings     = c(S, list(
      scenario = sim$scenario,
      formula  = deparse1(fmla),
      rho_hat  = rho_hat,
      sig2_hat = sig2,
      sp_hat   = as.list(sp_named)
    ))
  )
  
  validate_fit_result(
    fit,
    n_obs       = sim$n,
    p           = sim$p,
    n_grid_1d   = length(sim$x_grid_1d),
    n_grid_2d   = length(sim$x_grid_2d$u) * length(sim$x_grid_2d$v),
    n_draws     = S$n_draws,
    has_spatial = TRUE,
    int_keys    = character(0)
  )
  
  fit
}
