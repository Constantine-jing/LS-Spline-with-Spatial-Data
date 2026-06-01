# ============================================================
# fit_inla_wrapper.R
# 
# Wraps R-INLA to produce the canonical fit_result schema for
# the Scenario A cross-method comparison.
# 
# Model:
#   y_i = mu + sum_j f_j(X_ij) + s(loc_i) + eps_i
#   eps_i ~ N(0, sigma2)
#   f_j      : RW2 on a 50-bin grid of [0,1]   (standard INLA P-spline analog)
#   s(loc_i) : SPDE(alpha=2)  -> Matern with nu = alpha - d/2 = 1 in 2D
# 
# Settings (locked for fairness):
#   bins_per_smooth = 50                      # bin count for RW2
#   alpha           = 2                       # nu = 1 in 2D
#   mesh            = auto via inla.mesh.2d   # see settings doc
#   pc.prec for f_j : P(sd > 1) = 0.01        # weakly informative
#   pc.matern        : P(range < 0.10) = 0.05
#                      P(sigma > 1.0)  = 0.05
#   n_draws         = 2000                    # via inla.posterior.sample
# 
# Joint posterior draws are obtained from inla.posterior.sample(),
# which gives true joint draws of latent field + hyperparameters.
# We then project draws onto the canonical eval grids.
# 
# WARNING: First-time calls into INLA can be slow (few minutes)
# because INLA caches mesh / FEM matrices. Subsequent fits at
# the same mesh are much faster.
# ============================================================

source("canonical_schema.R")

suppressMessages({
  library(INLA)
  library(Matrix)
})


fit_inla <- function(sim, settings = NULL) {
  
  # ---- locked-in Stage 1 settings ----
  defaults <- list(
    bins_per_smooth = 50L,
    alpha           = 2,        # nu = 1 in 2D
    
    # mesh tuning. Inner mesh covers [0,1]^2; outer offset prevents
    # boundary effects from biasing the posterior at the unit-square edges.
    mesh_max_edge_inner = 0.10,
    mesh_max_edge_outer = 0.30,
    mesh_offset         = c(0.10, 0.30),
    mesh_cutoff         = 0.05,
    
    # PC priors
    pc_prec_f_u   = 1.0,    # P(sd > pc_prec_f_u) = pc_prec_f_alpha
    pc_prec_f_alpha = 0.01,
    pc_range_u    = 0.10,   # P(range < pc_range_u) = pc_range_alpha
    pc_range_alpha = 0.05,
    pc_sigma_u    = 1.0,    # P(spatial sd > pc_sigma_u) = pc_sigma_alpha
    pc_sigma_alpha = 0.05,
    
    # Likelihood precision: PC prior on residual sd
    pc_prec_y_u     = 1.0,
    pc_prec_y_alpha = 0.01,
    
    n_draws = 2000L,
    verbose = FALSE
  )
  if (!is.null(settings)) {
    for (nm in names(settings)) defaults[[nm]] <- settings[[nm]]
  }
  S <- defaults
  
  if (abs(2 - S$alpha) > 1e-12) {
    warning(sprintf(
      "fit_inla hardcodes alpha=2 (nu=1 in 2D); got alpha=%g; result is biased",
      S$alpha))
  }
  
  t0 <- proc.time()
  
  # ---- assemble inputs ----
  X_raw <- as.matrix(sim$data[, paste0("X", seq_len(sim$p)), drop = FALSE])
  loc   <- as.matrix(sim$data[, c("lon", "lat")])
  y     <- sim$data$y
  n     <- sim$n
  
  # ---- bin each X_j onto a 50-grid of [0,1] for RW2 ----
  # The bin centers ARE the discrete locations of the RW2.
  # We map data points -> nearest bin index, fit RW2 there,
  # then map posterior back to sim$x_grid_1d via bin centers.
  bin_centers <- seq(1 / (2 * S$bins_per_smooth), 1, length.out = S$bins_per_smooth)
  bin_edges   <- seq(0, 1, length.out = S$bins_per_smooth + 1)
  
  bin_idx_of <- function(x) {
    # bin index 1..bins_per_smooth, both endpoints inclusive at 1 and last bin
    idx <- findInterval(x, bin_edges, rightmost.closed = TRUE)
    pmax(1L, pmin(S$bins_per_smooth, idx))
  }
  
  X_bin <- apply(X_raw, 2, bin_idx_of)        # (n x p) integer bin indices
  colnames(X_bin) <- paste0("Xbin", seq_len(sim$p))
  
  # eval-grid bin indices (for predicting f_j on sim$x_grid_1d)
  grid_bin_idx <- bin_idx_of(sim$x_grid_1d)   # length n_grid_1d, in 1..bins
  
  # ---- mesh, SPDE, projector matrix ----
  # Build a 2D mesh covering [0,1]^2 with an outer offset.
  loc_for_mesh <- loc
  mesh <- inla.mesh.2d(
    loc      = loc_for_mesh,
    max.edge = c(S$mesh_max_edge_inner, S$mesh_max_edge_outer),
    offset   = S$mesh_offset,
    cutoff   = S$mesh_cutoff
  )
  
  spde <- inla.spde2.pcmatern(
    mesh      = mesh,
    alpha     = S$alpha,
    prior.range = c(S$pc_range_u, S$pc_range_alpha),
    prior.sigma = c(S$pc_sigma_u, S$pc_sigma_alpha)
  )
  
  # Projector matrices: A_obs maps mesh nodes -> observation locations
  A_obs  <- inla.spde.make.A(mesh = mesh, loc = loc)
  spde_idx <- inla.spde.make.index("s_field", n.spde = spde$n.spde)
  
  # ---- inla.stack: covariates as RW2 indices, spatial via projector ----
  # Effects:
  #   - intercept
  #   - per-smooth bin index (one effect per smooth, named Xbin1..Xbinp)
  #   - spatial latent field via SPDE indexing
  effects_obs <- list(
    list(intercept = rep(1, n)),
    spde_idx
  )
  for (j in seq_len(sim$p)) {
    nm_j <- paste0("Xbin", j)
    eff_j <- list()
    eff_j[[nm_j]] <- as.integer(X_bin[, j])
    effects_obs <- c(effects_obs, list(eff_j))
  }
  
  # A list: identity for the scalar effects and bin-index effects,
  # A_obs for the spatial field.
  A_list <- list(1, A_obs)
  for (j in seq_len(sim$p)) A_list <- c(A_list, list(1))
  
  stack_obs <- inla.stack(
    tag      = "obs",
    data     = list(y = y),
    A        = A_list,
    effects  = effects_obs,
    compress = TRUE
  )
  
  # ---- formula ----
  pc_prec_f <- list(prec = list(prior = "pc.prec",
                                param = c(S$pc_prec_f_u, S$pc_prec_f_alpha)))
  smooth_terms <- sapply(seq_len(sim$p), function(j) {
    sprintf("f(Xbin%d, model='rw2', constr=TRUE, scale.model=TRUE, hyper=hyper_pc_prec_f)",
            j)
  })
  rhs <- paste(c("intercept",
                 smooth_terms,
                 "f(s_field, model=spde)"),
               collapse = " + ")
  fmla <- as.formula(paste("y ~ -1 +", rhs))
  
  # Pass the hyperprior list into the fit env via assign in this scope
  hyper_pc_prec_f <- pc_prec_f
  
  # ---- fit ----
  t_fit_start <- proc.time()
  fit_inla_obj <- inla(
    formula = fmla,
    data    = inla.stack.data(stack_obs),
    family  = "gaussian",
    control.predictor = list(A = inla.stack.A(stack_obs), compute = TRUE),
    control.compute   = list(config = TRUE),     # required for posterior.sample
    control.family    = list(hyper = list(prec = list(
      prior = "pc.prec",
      param = c(S$pc_prec_y_u, S$pc_prec_y_alpha)
    ))),
    verbose = S$verbose
  )
  fit_sec <- (proc.time() - t_fit_start)[3]
  
  # ---- joint posterior draws ----
  t_post_start <- proc.time()
  
  # set.seed inside INLA for reproducibility of the draws
  inla.seed <- as.integer(sim$seed) + 7000L
  set.seed(inla.seed)
  pst <- inla.posterior.sample(n = S$n_draws, result = fit_inla_obj,
                               seed = inla.seed)
  
  # ---- extract main-effect smooths on sim$x_grid_1d ----
  # Each draw's `latent` is a single big named vector. Bin-j RW2
  # values are stored under names like "Xbin1:1", "Xbin1:2", ...
  # We pull them, then index by grid_bin_idx to get f_j on the grid.
  pull_block <- function(draw, prefix) {
    nm  <- rownames(draw$latent)
    msk <- grepl(paste0("^", prefix, ":"), nm)
    as.numeric(draw$latent[msk, 1])
  }
  
  n_grid_1d <- length(sim$x_grid_1d)
  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    blk <- matrix(NA_real_, nrow = S$bins_per_smooth, ncol = S$n_draws)
    for (k in seq_len(S$n_draws)) {
      blk[, k] <- pull_block(pst[[k]], sprintf("Xbin%d", j))
    }
    # Map bin-level draws -> eval-grid draws via nearest bin
    f_main[[j]] <- blk[grid_bin_idx, , drop = FALSE]   # (n_grid_1d x n_draws)
    
    # Center per-draw on the eval grid, matching the convention used
    # by fit_ours and fit_mgcv. RW2 in INLA already has constr=TRUE
    # which forces sum-to-zero across BIN positions (uniform weight),
    # not across the eval grid (where one bin gets ~5 grid points).
    # Re-center to the eval grid for like-for-like comparison.
    col_means <- colMeans(f_main[[j]])
    f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
  }
  
  # ---- spatial RE at observed locations ----
  # SPDE field draws are stored under "s_field:k" for mesh nodes 1..n.spde.
  # Project to observed locations via A_obs (n x n.spde).
  s_obs <- matrix(NA_real_, nrow = n, ncol = S$n_draws)
  for (k in seq_len(S$n_draws)) {
    s_field_draw <- pull_block(pst[[k]], "s_field")
    s_obs[, k] <- as.numeric(A_obs %*% s_field_draw)
  }
  
  # ---- variance components ----
  # INLA returns precisions; convert each draw to variance.
  # We pull them from the hyperpar block of each draw.
  hpar_names <- function(draw) names(draw$hyperpar)
  hp_first <- hpar_names(pst[[1]])
  
  find_hp <- function(pattern) {
    idx <- grep(pattern, hp_first)
    if (length(idx) == 0) return(NA_integer_)
    idx[1]
  }
  
  i_prec_y    <- find_hp("Precision for the Gaussian observations")
  i_prec_f_j  <- sapply(seq_len(sim$p),
                        function(j) find_hp(sprintf("Precision for Xbin%d", j)))
  i_range_s   <- find_hp("Range for s_field")
  i_stdev_s   <- find_hp("Stdev for s_field")
  
  if (is.na(i_prec_y))   stop("Could not locate 'Precision for the Gaussian observations' hyperpar")
  if (any(is.na(i_prec_f_j))) {
    msg <- paste(hp_first, collapse = "  |  ")
    stop("Could not locate per-smooth precisions; available hyperpars:\n  ", msg)
  }
  if (is.na(i_range_s) || is.na(i_stdev_s)) {
    msg <- paste(hp_first, collapse = "  |  ")
    stop("Could not locate SPDE range/stdev; available hyperpars:\n  ", msg)
  }
  
  hp_mat <- sapply(pst, function(d) as.numeric(d$hyperpar))
  # hp_mat is (n_hp x n_draws)
  
  sigma2_draws  <- 1 / hp_mat[i_prec_y, ]
  tau2_j_draws  <- lapply(i_prec_f_j, function(i) 1 / hp_mat[i, ])
  tau2_s_draws  <- (hp_mat[i_stdev_s, ])^2
  rho_draws     <- hp_mat[i_range_s, ]
  
  # NOTE on rho parameterization mismatch with our wrapper:
  # INLA-SPDE's "range" is the "practical range" ~= sqrt(8*nu)/kappa,
  # i.e. distance at which the Matern correlation is roughly 0.13.
  # Our matern_cor() uses rho such that cor = (z^nu * K_nu(z)) form
  # with z = d/rho. These are NOT the same parameter; conversion
  # documented in the settings doc. Reported as-is here.
  
  var_comp <- list(
    sigma2   = as.numeric(sigma2_draws),
    tau2     = lapply(tau2_j_draws, as.numeric),
    tau2_s   = as.numeric(tau2_s_draws),
    tau2_int = list(),
    rho      = as.numeric(rho_draws),
    nu       = NULL          # fixed via alpha=2
  )
  
  post_sec <- (proc.time() - t_post_start)[3]
  total_sec <- (proc.time() - t0)[3]
  
  # ---- assemble ----
  fit <- list(
    method       = "inla",
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
      rhat        = NULL,
      ess         = NULL,
      mlik        = as.numeric(fit_inla_obj$mlik[1, 1]),
      n_mesh      = mesh$n
    ),
    settings     = c(S, list(
      scenario   = sim$scenario,
      formula    = deparse1(fmla),
      n_mesh     = mesh$n,
      bin_centers = bin_centers
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
