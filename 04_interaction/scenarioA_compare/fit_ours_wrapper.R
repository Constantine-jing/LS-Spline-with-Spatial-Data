# ============================================================
# fit_ours_wrapper.R
# 
# Wraps gibbs_stage_c_full.R to produce the canonical
# fit_result schema for the Scenario A cross-method comparison.
# 
# This is a THIN shell: it does NOT modify the sampler. It only
# (a) builds the LS design matrix and distance matrix,
# (b) calls gibbs_full_sampler() with the locked Stage 1 settings,
# (c) maps each posterior sample of the coefficient vector into
#     joint draws of f_j on the canonical 1D evaluation grid,
# (d) renames variables to match the canonical schema,
# (e) validates the output.
# 
# IMPORTANT NAMING NOTE
# ---------------------
# gibbs_stage_c_full uses the convention:
#     y = H eta + b + eps
#     b   ~ N(0, sigma2 * R(rho, nu))   [SPATIAL variance]
#     eps ~ N(0, tau2 * I)              [RESIDUAL variance]
# The canonical schema uses the convention:
#     b   ~ N(0, tau2_s * R(rho, nu))   [SPATIAL variance]
#     eps ~ N(0, sigma2 * I)            [RESIDUAL variance]
# i.e. sigma2 and tau2 are SWAPPED.  This wrapper performs the
# rename so downstream metric code sees a consistent meaning.
# ============================================================

# These three sources are expected in the working directory.
# Adjust paths if your folder layout differs.
source("ls_basis.R")
source("spatial_utils.R")
source("gibbs_stage_c_full.R")
source("canonical_schema.R")


# ------------------------------------------------------------
# fit_ours(sim, settings = NULL)
# 
# Args:
#   sim       : object returned by simulate_scenario_A().
#   settings  : optional named list overriding any of the
#               locked-in defaults below. Provided so we can
#               sweep settings later without editing this file.
# 
# Returns:
#   fit_result conforming to canonical_schema.R.
# ------------------------------------------------------------
fit_ours <- function(sim, settings = NULL) {
  
  # ---- locked-in Stage 1 settings ----
  defaults <- list(
    M           = 20L,        # knot count per smooth (matches Sim1-Sim4+ default)
    nu          = 1.0,        # MUST match sim$truth_params$nu (Scenario A: 1.0)
    
    # MCMC budget chosen so post-burn, post-thin we get exactly 2000 draws.
    # 2500 burn-in + 4 * 2000 = 10500 total iter, thin by 4 -> 2000 draws.
    n_iter      = 10500L,
    n_burn      = 2500L,
    n_thin      = 4L,
    n_draws     = 2000L,      # MUST equal floor((n_iter - n_burn) / n_thin)
    
    # Priors (matches Sim1-Sim4+ baseline)
    kappa2      = 1e6,
    a_sigma     = 2.0, b_sigma = 1.0,
    a_tau       = 2.0, b_tau   = 0.3,
    a_smooth    = 1.0, b_smooth = 0.005,
    log_rho_mu  = log(0.2), log_rho_sd = 1.0,
    
    # MH proposal sds (matches Sim1-Sim4+ tuning)
    mh_sd_log_sigma2 = 0.5,
    mh_sd_log_tau2   = 0.3,
    mh_sd_log_rho    = 0.3,
    
    eps_ridge = 1e-6,
    jitter    = 1e-8,
    verbose   = TRUE
  )
  if (!is.null(settings)) {
    for (nm in names(settings)) defaults[[nm]] <- settings[[nm]]
  }
  S <- defaults
  
  # Sanity: the post-burn, post-thin count must equal n_draws
  expected_draws <- floor((S$n_iter - S$n_burn) / S$n_thin)
  if (expected_draws != S$n_draws) {
    stop(sprintf(
      "n_iter=%d, n_burn=%d, n_thin=%d gives %d draws but n_draws=%d",
      S$n_iter, S$n_burn, S$n_thin, expected_draws, S$n_draws))
  }
  if (abs(S$nu - sim$truth_params$nu) > 1e-12) {
    stop(sprintf("nu mismatch: settings$nu=%g vs sim$truth_params$nu=%g",
                 S$nu, sim$truth_params$nu))
  }
  
  t0 <- proc.time()
  
  # ---- design matrices ----
  X_raw  <- as.matrix(sim$data[, paste0("X", seq_len(sim$p)), drop = FALSE])
  coords <- as.matrix(sim$data[, c("lon", "lat")])
  y      <- sim$data$y
  D      <- pairdist(coords)
  
  des <- ls_additive_build(X_raw, M_vec = rep(S$M, sim$p))
  H   <- cbind(1, des$W)
  
  # ---- design matrix for the 1D evaluation grid (per smooth) ----
  # For smooth j we want a length-(n_grid_1d) vector of f_j values
  # for every posterior draw. We build a per-smooth design block
  # W_grid_j by setting all OTHER covariates to a fixed dummy value
  # (irrelevant -- we only multiply by the columns of beta_j) and
  # taking the j-th block of columns from ls_additive_design_new.
  # 
  # Equivalent and simpler: use the per-smooth recipe stored in
  # des$objs[[j]]$design_new directly.
  n_grid_1d <- length(sim$x_grid_1d)
  W_grid_list <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    W_grid_list[[j]] <- des$objs[[j]]$design_new(sim$x_grid_1d,
                                                 type = "W", clip = TRUE)
  }
  
  # ---- REML init (same pattern as Sim1) ----
  source("fit_spatial_reml.R")
  obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                         M_vec = rep(6, sim$p), nu = S$nu,
                         rho_init    = sim$truth_params$rho,
                         lambda_init = sim$truth_params$sigma2 /
                                       max(sim$truth_params$tau2_s, 0.01),
                         verbose = FALSE)
  R_init     <- matern_cor(D, rho = obj6$fit$rho, nu = S$nu)
  Sigma_init <- max(obj6$fit$sigma2, 0.01) * R_init +
                max(obj6$fit$tau2,   0.01) * diag(sim$n)
  L_init     <- chol(Sigma_init + diag(1e-8, sim$n))
  y_w        <- forwardsolve(t(L_init), y)
  H_w        <- forwardsolve(t(L_init), H)
  eta_init   <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H)),
                                 crossprod(H_w, y_w)))
  init <- list(eta    = eta_init,
               sigma2 = max(obj6$fit$sigma2, 0.01),
               tau2   = max(obj6$fit$tau2,   0.01),
               rho    = obj6$fit$rho,
               tau2_s = rep(1.0, sim$p))
  
  # ---- Gibbs ----
  t_fit_start <- proc.time()
  gs <- gibbs_full_sampler(
    y       = y,
    H       = H,
    D       = D,
    nu      = S$nu,
    col_map = des$col_map,
    n_iter  = S$n_iter,
    n_burn  = S$n_burn,
    n_thin  = S$n_thin,
    kappa2  = S$kappa2,
    a_sigma = S$a_sigma, b_sigma = S$b_sigma,
    a_tau   = S$a_tau,   b_tau   = S$b_tau,
    a_smooth = S$a_smooth, b_smooth = S$b_smooth,
    log_rho_mu = S$log_rho_mu, log_rho_sd = S$log_rho_sd,
    mh_sd_log_sigma2 = S$mh_sd_log_sigma2,
    mh_sd_log_tau2   = S$mh_sd_log_tau2,
    mh_sd_log_rho    = S$mh_sd_log_rho,
    eps_ridge = S$eps_ridge,
    init      = init,
    jitter    = S$jitter,
    verbose   = S$verbose
  )
  fit_sec <- (proc.time() - t_fit_start)[3]
  
  # ---- Build joint posterior draws of f_j on the eval grid ----
  # gs$eta_samples is (n_draws x p_total) where p_total = 1 + sum(M_j - 1).
  # The intercept is column 1; smooth-j coefficients are columns
  # 1 + des$col_map[[j]] (note the +1 offset because mu is column 1).
  # 
  # f_j_grid[k, ] = W_grid_j %*% beta_j_k  (row vector across grid pts)
  # We want (n_grid_1d x n_draws) so transpose.
  t_post_start <- proc.time()
  
  n_draws <- S$n_draws
  stopifnot(nrow(gs$eta_samples) == n_draws)
  
  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    cols_j <- 1L + des$col_map[[j]]
    beta_j_draws <- gs$eta_samples[, cols_j, drop = FALSE]   # (n_draws x M_j-1)
    Wj <- W_grid_list[[j]]                                   # (n_grid_1d x M_j-1)
    # f_j on grid: (n_grid_1d x n_draws)
    f_main[[j]] <- Wj %*% t(beta_j_draws)
    
    # IDENTIFIABILITY: the LS basis for each smooth has its own
    # implicit intercept absorbed into mu, but per-draw f_j is not
    # automatically sum-to-zero on the eval grid. Center per-draw
    # so RMSE/coverage are computed against the CENTERED truth.
    col_means <- colMeans(f_main[[j]])
    f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
  }
  
  # ---- Spatial RE draws at observed locations ----
  # gs$b_samples is (n_draws x n).  Transpose to (n x n_draws).
  s_obs <- t(gs$b_samples)
  
  # ---- Variance components, with the sigma2/tau2 RENAME ----
  # gs$sigma2_samples = SPATIAL variance  -> canonical tau2_s
  # gs$tau2_samples   = RESIDUAL variance -> canonical sigma2
  var_comp <- list(
    sigma2   = as.numeric(gs$tau2_samples),         # residual
    tau2     = lapply(seq_len(sim$p),
                      function(j) as.numeric(gs$tau2_s_samples[, j])),
    tau2_s   = as.numeric(gs$sigma2_samples),       # spatial
    tau2_int = list(),                              # empty for Scenario A
    rho      = as.numeric(gs$rho_samples),
    nu       = NULL                                 # fixed, not sampled
  )
  
  # ---- Convergence diagnostics (R-hat needs >=2 chains; skip
  # for Stage 1 single-chain setup. Compute ESS, basic checks) ----
  ess_safe <- function(x) {
    if (length(x) < 2) return(NA_real_)
    # simple Geyer initial-positive-sequence ESS via coda if available
    if (requireNamespace("coda", quietly = TRUE)) {
      as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
    } else {
      NA_real_
    }
  }
  ess_named <- c(
    sigma2 = ess_safe(var_comp$sigma2),
    tau2_s = ess_safe(var_comp$tau2_s),
    rho    = ess_safe(var_comp$rho)
  )
  for (j in seq_len(sim$p)) {
    ess_named[paste0("tau2_", j)] <- ess_safe(var_comp$tau2[[j]])
  }
  
  post_sec <- (proc.time() - t_post_start)[3]
  total_sec <- (proc.time() - t0)[3]
  
  # ---- assemble canonical fit_result ----
  fit <- list(
    method       = "ours",
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
      rhat = NULL,        # single chain, not computed here
      ess  = ess_named,
      mh_accept = gs$accept_rate
    ),
    settings     = c(S, list(scenario = sim$scenario))
  )
  
  # ---- validate before returning ----
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
