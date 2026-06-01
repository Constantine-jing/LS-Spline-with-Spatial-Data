# ============================================================
# fit_ours_interaction_wrapper.R   (patched: orthogonalize = TRUE)
#
# Same as before EXCEPT:
#   - Adds an `orthogonalize` setting (default TRUE) that is forwarded to
#     ls_build_all_interactions().
#   - Saves the setting into the returned fit$settings for traceability.
#
# To run the OLD behavior (basis-correlation-vulnerable), pass
#     settings = list(orthogonalize = FALSE)
# ============================================================

source("ls_basis.R")
source("ls_interaction.R")
source("spatial_utils.R")
source("gibbs_stage_c_full.R")
source("gibbs_interaction.R")
source("canonical_schema.R")


fit_ours_interaction <- function(sim, settings = NULL) {

  defaults <- list(
    M           = 20L,
    nu          = 1.0,

    n_iter      = 10500L,
    n_burn      = 2500L,
    n_thin      = 4L,
    n_draws     = 2000L,

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
    verbose   = TRUE,

    # NEW: design-level ANOVA orthogonalization of W_uv against [1, W_u, W_v].
    # TRUE  : Project 1 fix (recommended)
    # FALSE : original (basis-correlation-vulnerable) behavior
    orthogonalize = TRUE,

    # NEW (May 2026): clamp the per-smooth main-effect smoothing variance.
    # NULL or all-NA  : conjugate IG update (default behaviour, no clamp).
    # numeric vec     : element j finite -> tau2_s_main[j] held fixed at that
    #                   value throughout MCMC (skips Step 5 IG draw for j).
    # Used for Hypothesis 3 confirmation: tau2_s_main_fixed = c(0.005, NA, NA)
    tau2_s_main_fixed = NULL,

    # NEW (May 2026): same idea for interaction blocks. Length = number of
    # interaction pairs (3 for p=3 if all pairs included).
    tau2_s_int_fixed = NULL,

    # NEW (May 2026): scalar clamps for variance components. NA -> sample.
    sigma2_fixed = NA_real_,
    tau2_fixed   = NA_real_,
    rho_fixed    = NA_real_,

    # NEW (May 2026, P1 cleanup):
    # TRUE  -> all-pairs interaction blocks are built and sampled (default).
    # FALSE -> the interaction-basis construction is skipped, H = [1 | W_main]
    #          only, and the inner sampler's Step 5b is bypassed.
    # `sim$int_keys` is still consulted (downstream code may reference it)
    # but no posterior is sampled for those keys; fit$f_int and
    # fit$var_comp$tau2_int will be empty lists.
    fit_interactions = TRUE
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

  X_raw  <- as.matrix(sim$data[, paste0("X", seq_len(sim$p)), drop = FALSE])
  coords <- as.matrix(sim$data[, c("lon", "lat")])
  y      <- sim$data$y
  D      <- pairdist(coords)

  full_objs <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    full_objs[[j]] <- ls_build_one_full(X_raw[, j], M = S$M)
  }

  W_main_blocks <- lapply(full_objs, function(o) o$W)
  W_main <- do.call(cbind, W_main_blocks)
  K_main_list <- lapply(full_objs, function(o) {
    build_rw2_penalty(S$M - 1L)
  })

  # ----- Interaction-basis construction (skipped if fit_interactions=FALSE) -----
  # Pass orthogonalize through to ls_build_all_interactions.
  if (S$fit_interactions) {
    int_list <- ls_build_all_interactions(full_objs,
                                          orthogonalize = S$orthogonalize)
    stopifnot(setequal(names(int_list), sim$int_keys))

    if (S$verbose) {
      cat(sprintf("[fit_ours_interaction] orthogonalize = %s\n", S$orthogonalize))
      # Quick post-build sanity: print the cross-Gram norms for the first pair
      pair1 <- int_list[[1]]
      W_uv_1 <- pair1$W_uv
      W_u_1  <- full_objs[[ as.integer(strsplit(names(int_list)[1],"_")[[1]][1]) ]]$W
      W_v_1  <- full_objs[[ as.integer(strsplit(names(int_list)[1],"_")[[1]][2]) ]]$W
      cat(sprintf("  pair %s: max |1'W_uv|=%.2e  max |W_u'W_uv|=%.2e  max |W_v'W_uv|=%.2e\n",
                  names(int_list)[1],
                  max(abs(colSums(W_uv_1))),
                  max(abs(crossprod(W_u_1, W_uv_1))),
                  max(abs(crossprod(W_v_1, W_uv_1)))))
    }

    W_int_blocks <- lapply(int_list, function(x) x$W_uv)
    K_int_list   <- lapply(int_list, function(x) x$K_uv)

    W_int_concat <- do.call(cbind, W_int_blocks)
    H <- cbind(1, W_main, W_int_concat)
  } else {
    if (S$verbose) {
      cat("[fit_ours_interaction] fit_interactions = FALSE, skipping interaction Gibbs step\n")
    }
    int_list     <- list()
    W_int_blocks <- list()
    K_int_list   <- list()
    H            <- cbind(1, W_main)
  }

  col_map_main <- vector("list", sim$p)
  off <- 1L  # start after intercept (η[1] is intercept; β starts at η[2])
  for (j in seq_len(sim$p)) {
    nj <- ncol(W_main_blocks[[j]])
    col_map_main[[j]] <- off + seq_len(nj) - 1L
    off <- off + nj
  }

  col_map_int <- vector("list", length(int_list))
  names(col_map_int) <- names(int_list)
  for (k in seq_along(int_list)) {
    nk <- ncol(W_int_blocks[[k]])
    col_map_int[[k]] <- off + seq_len(nk) - 1L
    off <- off + nk
  }

  W_grid_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    W_grid_main[[j]] <- full_objs[[j]]$design_new(sim$x_grid_1d,
                                                   type = "W", clip = TRUE)
  }

  flat2d <- sim$flat_grid_2d
  n_grid_2d <- nrow(flat2d)
  W_grid_int <- vector("list", length(int_list))
  names(W_grid_int) <- names(int_list)
  if (S$fit_interactions) {
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
    }
  }

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
               tau2_s_int  = if (length(int_list) > 0L)
                               rep(1.0, length(int_list))
                             else
                               numeric(0))

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
    verbose   = S$verbose,
    tau2_s_main_fixed = S$tau2_s_main_fixed,
    tau2_s_int_fixed  = S$tau2_s_int_fixed,
    sigma2_fixed      = S$sigma2_fixed,
    tau2_fixed        = S$tau2_fixed,
    rho_fixed         = S$rho_fixed,
    fit_interactions  = S$fit_interactions
  )
  fit_sec <- (proc.time() - t_fit_start)[3]

  t_post_start <- proc.time()
  n_draws <- S$n_draws
  stopifnot(nrow(gs$eta_samples) == n_draws)

  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    cols_j <- 1L + col_map_main[[j]]
    beta_j_draws <- gs$eta_samples[, cols_j, drop = FALSE]
    f_main[[j]] <- W_grid_main[[j]] %*% t(beta_j_draws)

    col_means <- colMeans(f_main[[j]])
    f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
  }

  f_int <- vector("list", length(int_list))
  names(f_int) <- names(int_list)
  for (k in seq_along(int_list)) {
    cols_k <- 1L + col_map_int[[k]]
    beta_k_draws <- gs$eta_samples[, cols_k, drop = FALSE]
    f_int[[k]] <- W_grid_int[[k]] %*% t(beta_k_draws)

    col_means_2d <- colMeans(f_int[[k]])
    f_int[[k]]   <- sweep(f_int[[k]], 2, col_means_2d, FUN = "-")
  }

  s_obs <- t(gs$b_samples)

  var_comp <- list(
    sigma2   = as.numeric(gs$tau2_samples),
    tau2     = lapply(seq_len(sim$p),
                      function(j) as.numeric(gs$tau2_s_main_samp[, j])),
    tau2_s   = as.numeric(gs$sigma2_samples),
    tau2_int = lapply(seq_along(int_list),
                      function(k) as.numeric(gs$tau2_s_int_samp[, k])),
    rho      = as.numeric(gs$rho_samples),
    nu       = NULL
  )
  names(var_comp$tau2_int) <- names(int_list)

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
      int_keys = if (length(int_list) > 0L) names(int_list) else character(0)
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
    int_keys    = if (S$fit_interactions) sim$int_keys else character(0)
  )

  fit
}
