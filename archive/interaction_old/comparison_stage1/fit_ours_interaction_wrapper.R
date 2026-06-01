# ============================================================
# fit_ours_interaction_wrapper.R
# 
# Wraps gibbs_interaction_sampler (from gibbs_interaction.R) to
# produce the canonical fit_result schema for Scenario B.
# 
# Model:
#   y = mu + sum_j f_j(X_j) + sum_{u<v} f_{uv}(X_u, X_v) + b + eps
#   b ~ N(0, tau2_s * Matern(rho, nu=1))     [marginalized out]
#   eps ~ N(0, sigma2 * I)
#   f_j     : LS basis (M=20) with RW2 prior
#   f_{uv}  : LS tensor-product basis (Khatri-Rao) with 2D Kronecker-sum RW2
# 
# Same naming flip as Scenario A:
#   gibbs_interaction.R uses sigma2=spatial, tau2=residual
#   canonical schema uses sigma2=residual, tau2_s=spatial
# Wrapper performs the rename when assembling var_comp.
# 
# All p(p-1)/2 = 3 interaction blocks are included in the model
# (full Phase 2 of run_interaction_2x2_v2.R). The truth has only
# f_{1,2} non-zero -- the model needs to find the others null.
# ============================================================

source("ls_basis.R")
source("ls_interaction.R")
source("spatial_utils.R")
source("gibbs_stage_c_full.R")        # build_rw2_penalty()
source("gibbs_interaction.R")
source("canonical_schema.R")


fit_ours_interaction <- function(sim, settings = NULL) {
  
  # ---- locked-in Stage 1 settings ----
  defaults <- list(
    M           = 20L,
    nu          = 1.0,
    
    # MCMC budget. Same as Scenario A (10500 / 2500 burn / thin 4 -> 2000 draws).
    # Interaction model has more parameters but per-iteration cost is similar
    # (the bottleneck is the n x n Cholesky inside the marginal likelihood).
    n_iter      = 10500L,
    n_burn      = 2500L,
    n_thin      = 4L,
    n_draws     = 2000L,
    
    # Priors -- match Scenario A baseline.
    kappa2      = 1e6,
    a_sigma     = 2.0, b_sigma = 1.0,
    a_tau       = 2.0, b_tau   = 0.3,
    a_smooth    = 1.0, b_smooth = 0.005,
    log_rho_mu  = log(0.2), log_rho_sd = 1.0,
    
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
  
  # ---- inputs ----
  X_raw  <- as.matrix(sim$data[, paste0("X", seq_len(sim$p)), drop = FALSE])
  coords <- as.matrix(sim$data[, c("lon", "lat")])
  y      <- sim$data$y
  D      <- pairdist(coords)
  
  # ---- main-effect bases (need the FULL objs for ls_build_interaction) ----
  # ls_additive_build returns compact objs WITHOUT $W stored, but
  # ls_build_interaction needs $W. We rebuild via ls_build_one_full
  # for each smooth -- same M as the additive build, so design matrices
  # are bit-identical.
  full_objs <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    full_objs[[j]] <- ls_build_one_full(X_raw[, j], M = S$M)
  }
  
  # Main-effect design block (n x sum(M-1)) and 1D RW2 penalties (M x M each).
  W_main_blocks <- lapply(full_objs, function(o) o$W)
  W_main <- do.call(cbind, W_main_blocks)
  K_main_list <- lapply(full_objs, function(o) {
    # K is the identified RW2 penalty: T^T K_raw T  where K_raw is M x M.
    # In gibbs_stage_c_full this is computed inside the sampler from
    # build_rw2_penalty(M-1) on the (M-1)-dim coefficient. To match,
    # we use the (M-1)-dim RW2 penalty directly.
    # (ls_basis.R doesn't store K; we build it here.)
    build_rw2_penalty(S$M - 1L)
  })
  
  # ---- interaction designs and 2D RW2 penalties ----
  int_list <- ls_build_all_interactions(full_objs)
  # int_list keys are "u_v" with u<v -- matches sim$int_keys.
  stopifnot(setequal(names(int_list), sim$int_keys))
  
  W_int_blocks <- lapply(int_list, function(x) x$W_uv)
  K_int_list   <- lapply(int_list, function(x) x$K_uv)
  
  # ---- assemble H and col_maps (0-indexed offsets after intercept) ----
  W_int_concat <- do.call(cbind, W_int_blocks)
  H <- cbind(1, W_main, W_int_concat)
  
  # main effect col offsets
  col_map_main <- vector("list", sim$p)
  off <- 0L
  for (j in seq_len(sim$p)) {
    nj <- ncol(W_main_blocks[[j]])
    col_map_main[[j]] <- off + seq_len(nj) - 1L
    off <- off + nj
  }
  
  # interaction col offsets (continue from where main left off)
  col_map_int <- vector("list", length(int_list))
  names(col_map_int) <- names(int_list)
  for (k in seq_along(int_list)) {
    nk <- ncol(W_int_blocks[[k]])
    col_map_int[[k]] <- off + seq_len(nk) - 1L
    off <- off + nk
  }
  
  # ---- prediction-grid design matrices ----
  # 1D: per-smooth design at sim$x_grid_1d
  W_grid_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    W_grid_main[[j]] <- full_objs[[j]]$design_new(sim$x_grid_1d,
                                                   type = "W", clip = TRUE)
  }
  
  # 2D: per-pair design at sim$flat_grid_2d (row-major u then v)
  flat2d <- sim$flat_grid_2d                # n_grid_2d x 2 with cols u,v
  n_grid_2d <- nrow(flat2d)
  W_grid_int <- vector("list", length(int_list))
  names(W_grid_int) <- names(int_list)
  for (k in seq_along(int_list)) {
    pair <- strsplit(names(int_list)[k], "_", fixed = TRUE)[[1]]
    u_idx <- as.integer(pair[1])
    v_idx <- as.integer(pair[2])
    rec   <- int_list[[k]]$recipe
    W_grid_int[[k]] <- ls_interaction_design_new(
      X_u_new = flat2d[, "u"],
      X_v_new = flat2d[, "v"],
      recipe  = rec,
      clip    = TRUE
    )
    # NB: the recipe maps (X_u, X_v) -> Khatri-Rao basis using the SAME
    # tau knots and contrasts as the training W_uv block. So column
    # ordering of W_grid_int[[k]] matches col_map_int[[k]] exactly.
  }
  
  # ---- REML init from non-interaction model ----
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
  init <- list(eta         = eta_init,
               sigma2      = max(obj6$fit$sigma2, 0.01),
               tau2        = max(obj6$fit$tau2,   0.01),
               rho         = obj6$fit$rho,
               tau2_s_main = rep(1.0, sim$p),
               tau2_s_int  = rep(1.0, length(int_list)))
  
  # ---- Gibbs ----
  t_fit_start <- proc.time()
  gs <- gibbs_interaction_sampler(
    y = y, H = H, D = D, nu = S$nu,
    col_map_main = col_map_main,
    K_main_list  = K_main_list,
    col_map_int  = col_map_int,
    K_int_list   = K_int_list,
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
  
  # ---- joint posterior draws of f_j on x_grid_1d ----
  t_post_start <- proc.time()
  
  n_draws <- S$n_draws
  stopifnot(nrow(gs$eta_samples) == n_draws)
  
  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    cols_j <- 1L + col_map_main[[j]]
    beta_j_draws <- gs$eta_samples[, cols_j, drop = FALSE]   # (n_draws x M-1)
    f_main[[j]] <- W_grid_main[[j]] %*% t(beta_j_draws)
    
    # ANOVA-type centering on the eval grid (same convention as Scenario A)
    col_means <- colMeans(f_main[[j]])
    f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
  }
  
  # ---- joint posterior draws of f_{uv} on flat_grid_2d ----
  f_int <- vector("list", length(int_list))
  names(f_int) <- names(int_list)
  for (k in seq_along(int_list)) {
    cols_k <- 1L + col_map_int[[k]]
    beta_k_draws <- gs$eta_samples[, cols_k, drop = FALSE]   # (n_draws x (M-1)^2)
    f_int[[k]] <- W_grid_int[[k]] %*% t(beta_k_draws)        # (n_grid_2d x n_draws)
    
    # ANOVA-type centering on the 2D eval grid: subtract per-draw global mean.
    # This is the 2D analog of the 1D centering above and is required for
    # like-for-like comparison with the centered truth.
    col_means_2d <- colMeans(f_int[[k]])
    f_int[[k]]   <- sweep(f_int[[k]], 2, col_means_2d, FUN = "-")
  }
  
  # ---- spatial RE draws ----
  s_obs <- t(gs$b_samples)
  
  # ---- variance components, naming flip applied ----
  # gs$sigma2_samples = SPATIAL variance  -> canonical tau2_s
  # gs$tau2_samples   = RESIDUAL variance -> canonical sigma2
  var_comp <- list(
    sigma2   = as.numeric(gs$tau2_samples),         # residual
    tau2     = lapply(seq_len(sim$p),
                      function(j) as.numeric(gs$tau2_s_main_samp[, j])),
    tau2_s   = as.numeric(gs$sigma2_samples),       # spatial
    tau2_int = lapply(seq_along(int_list),
                      function(k) as.numeric(gs$tau2_s_int_samp[, k])),
    rho      = as.numeric(gs$rho_samples),
    nu       = NULL
  )
  names(var_comp$tau2_int) <- names(int_list)
  
  # ---- ESS diagnostics ----
  ess_safe <- function(x) {
    if (length(x) < 2) return(NA_real_)
    if (requireNamespace("coda", quietly = TRUE)) {
      as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
    } else NA_real_
  }
  ess_named <- c(
    sigma2 = ess_safe(var_comp$sigma2),
    tau2_s = ess_safe(var_comp$tau2_s),
    rho    = ess_safe(var_comp$rho)
  )
  for (j in seq_len(sim$p)) {
    ess_named[paste0("tau2_", j)] <- ess_safe(var_comp$tau2[[j]])
  }
  for (k in names(int_list)) {
    ess_named[paste0("tau2_int_", k)] <- ess_safe(var_comp$tau2_int[[k]])
  }
  
  post_sec <- (proc.time() - t_post_start)[3]
  total_sec <- (proc.time() - t0)[3]
  
  # ---- assemble canonical fit_result ----
  fit <- list(
    method       = "ours",
    scenario     = sim$scenario,
    seed         = sim$seed,
    f_main       = f_main,
    f_int        = f_int,
    s_obs        = s_obs,
    var_comp     = var_comp,
    x_grid_1d    = sim$x_grid_1d,
    x_grid_2d    = sim$x_grid_2d,
    timing       = list(total_sec = as.numeric(total_sec),
                        fit_sec   = as.numeric(fit_sec),
                        post_sec  = as.numeric(post_sec)),
    convergence  = list(
      rhat      = NULL,
      ess       = ess_named,
      mh_accept = gs$accept_rate
    ),
    settings     = c(S, list(
      scenario = sim$scenario,
      int_keys = names(int_list)
    ))
  )
  
  validate_fit_result(
    fit,
    n_obs       = sim$n,
    p           = sim$p,
    n_grid_1d   = length(sim$x_grid_1d),
    n_grid_2d   = nrow(flat2d),
    n_draws     = S$n_draws,
    has_spatial = TRUE,
    int_keys    = sim$int_keys
  )
  
  fit
}
