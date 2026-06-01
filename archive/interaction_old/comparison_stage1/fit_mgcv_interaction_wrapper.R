# ============================================================
# fit_mgcv_interaction_wrapper.R
# 
# Wraps mgcv::gam with tensor-product interactions for the
# Scenario B comparison.
# 
# Model:
#   y = mu + sum_j s(X_j; bs="ps", k=k_main)
#         + sum_{u<v} ti(X_u, X_v; bs="ps", k=c(k_int_per, k_int_per))
#         + s(lon, lat; bs="gp", m=c(2))
#         + eps
# 
# Smooth bases (locked for fairness):
#   - Main effects:  bs = "ps", k = k_main, m = c(2,2)
#   - Interactions:  ti() with bs = "ps", k = c(k_int_per, k_int_per),
#                    m = list(c(2,2), c(2,2))
#                    -> pure interaction (main effects subtracted by mgcv);
#                       gives ANOVA-orthogonal f_uv (matches our ANOVA-style
#                       interpretation of LS tensor product).
#   - Spatial:       bs = "gp", m = c(2) -> Matern nu = 1.0
# 
# Joint posterior draws via MASS::mvrnorm(coef, Vp).
# 
# Note on knot counts:
#   For comparison with our M=20 result, set k_main=20, k_int_per=20.
#   For comparison with our M=8 result,  set k_main=8,  k_int_per=8.
#   ti() with k=c(K,K) gives (K-1)^2 basis columns -- matches our LS
#   tensor product which has (M-1)^2 columns.
# ============================================================

source("canonical_schema.R")

suppressMessages({
  library(mgcv)
  library(MASS)
})


fit_mgcv_interaction <- function(sim, settings = NULL) {
  
  # ---- locked-in Stage 1 settings ----
  defaults <- list(
    k_main      = 20L,
    k_int_per   = 20L,        # ti(X1,X2,k=c(20,20)) gives 19^2 = 361 cols
    bs_main     = "ps",
    m_main      = c(2L, 2L),
    bs_int      = "ps",
    m_int_per   = c(2L, 2L),
    
    bs_spatial  = "gp",
    m_spatial   = c(2),
    k_spatial   = 100L,
    
    method      = "REML",
    n_draws     = 2000L
  )
  if (!is.null(settings)) {
    for (nm in names(settings)) defaults[[nm]] <- settings[[nm]]
  }
  S <- defaults
  
  if (abs(1.0 - sim$truth_params$nu) > 1e-12) {
    warning(sprintf(
      "fit_mgcv_interaction hardcodes nu=1; got truth nu=%g",
      sim$truth_params$nu))
  }
  
  t0 <- proc.time()
  
  # ---- assemble data ----
  df <- sim$data
  smooth_names <- paste0("X", seq_len(sim$p))
  
  # ---- build formula ----
  smooth_terms <- sapply(smooth_names, function(nm) {
    sprintf("s(%s, bs='%s', k=%d, m=c(%d,%d))",
            nm, S$bs_main, S$k_main, S$m_main[1], S$m_main[2])
  })
  
  # ti() interaction terms for every pair (u,v) with u<v
  pair_keys <- sim$int_keys     # e.g. "1_2", "1_3", "2_3"
  ti_terms <- sapply(pair_keys, function(k) {
    parts <- strsplit(k, "_", fixed = TRUE)[[1]]
    u <- parts[1]; v <- parts[2]
    sprintf(
      "ti(X%s, X%s, bs='%s', k=c(%d,%d), m=list(c(%d,%d), c(%d,%d)))",
      u, v, S$bs_int, S$k_int_per, S$k_int_per,
      S$m_int_per[1], S$m_int_per[2], S$m_int_per[1], S$m_int_per[2]
    )
  })
  names(ti_terms) <- pair_keys
  
  spatial_term <- sprintf("s(lon, lat, bs='%s', k=%d, m=c(%d))",
                          S$bs_spatial, S$k_spatial, S$m_spatial[1])
  rhs <- paste(c(smooth_terms, ti_terms, spatial_term), collapse = " + ")
  fmla <- as.formula(paste("y ~", rhs))
  
  # ---- fit ----
  t_fit_start <- proc.time()
  gam_fit <- gam(fmla, data = df, method = S$method)
  fit_sec <- (proc.time() - t_fit_start)[3]
  
  # ---- joint posterior draws ----
  t_post_start <- proc.time()
  
  beta_hat <- coef(gam_fit)
  Vp       <- gam_fit$Vp
  set.seed(sim$seed + 1000L)
  beta_draws <- MASS::mvrnorm(n = S$n_draws, mu = beta_hat, Sigma = Vp)
  
  # ---- locate smooth blocks by label ----
  smooth_blocks <- vector("list", length(gam_fit$smooth))
  names(smooth_blocks) <- sapply(gam_fit$smooth, function(s) s$label)
  for (k in seq_along(gam_fit$smooth)) {
    sm <- gam_fit$smooth[[k]]
    smooth_blocks[[k]] <- sm$first.para:sm$last.para
  }
  
  n_grid_1d <- length(sim$x_grid_1d)
  n_grid_2d <- nrow(sim$flat_grid_2d)
  
  # ---- main-effect smooths on x_grid_1d ----
  template <- df[1, , drop = FALSE]
  for (nm in colnames(df)) {
    if (is.numeric(df[[nm]])) template[[nm]] <- mean(df[[nm]])
  }
  
  f_main <- vector("list", sim$p)
  for (j in seq_len(sim$p)) {
    nm_j <- smooth_names[j]
    nd <- template[rep(1, n_grid_1d), , drop = FALSE]
    nd[[nm_j]] <- sim$x_grid_1d
    Xp <- predict(gam_fit, newdata = nd, type = "lpmatrix")
    
    label_j <- sprintf("s(%s)", nm_j)
    if (!label_j %in% names(smooth_blocks)) {
      stop(sprintf("Could not find smooth block '%s'", label_j))
    }
    cols_j <- smooth_blocks[[label_j]]
    
    Xp_j <- Xp[, cols_j, drop = FALSE]
    bj   <- beta_draws[, cols_j, drop = FALSE]
    f_main[[j]] <- Xp_j %*% t(bj)
    
    # NOTE on centering: matching our wrapper's convention (now also
    # NO per-draw centering, after the diagnostic showed it didn't
    # change anything materially). mgcv's s() is sum-to-zero by
    # construction at the *training* points (not at the eval grid),
    # so we do NOT re-center here. Pointwise RMSE is computed
    # against the truth, which is centered analytically -- both
    # mgcv and ours use this same convention.
  }
  
  # ---- interaction surfaces on flat_grid_2d ----
  flat2d <- sim$flat_grid_2d
  f_int <- vector("list", length(pair_keys))
  names(f_int) <- pair_keys
  
  for (k in pair_keys) {
    parts <- strsplit(k, "_", fixed = TRUE)[[1]]
    u <- as.integer(parts[1])
    v <- as.integer(parts[2])
    nm_u <- paste0("X", u)
    nm_v <- paste0("X", v)
    
    nd <- template[rep(1, n_grid_2d), , drop = FALSE]
    nd[[nm_u]] <- flat2d[, "u"]
    nd[[nm_v]] <- flat2d[, "v"]
    Xp <- predict(gam_fit, newdata = nd, type = "lpmatrix")
    
    # ti() label is "ti(Xu,Xv)" -- check version-dependent spacing
    label_uv <- sprintf("ti(%s,%s)", nm_u, nm_v)
    if (!label_uv %in% names(smooth_blocks)) {
      # try other possible spacings
      cand <- sprintf("ti(%s, %s)", nm_u, nm_v)
      if (cand %in% names(smooth_blocks)) {
        label_uv <- cand
      } else {
        stop(sprintf(
          "Could not find ti() block; available labels: %s",
          paste(names(smooth_blocks), collapse = " | ")))
      }
    }
    cols_uv <- smooth_blocks[[label_uv]]
    
    Xp_uv <- Xp[, cols_uv, drop = FALSE]
    b_uv  <- beta_draws[, cols_uv, drop = FALSE]
    f_int[[k]] <- Xp_uv %*% t(b_uv)
  }
  
  # ---- spatial RE at observed locations ----
  Xp_obs <- predict(gam_fit, newdata = df, type = "lpmatrix")
  label_s <- "s(lon,lat)"
  if (!label_s %in% names(smooth_blocks)) {
    cand <- "s(lon, lat)"
    if (cand %in% names(smooth_blocks)) {
      label_s <- cand
    } else {
      stop(sprintf("Could not find spatial smooth block; available: %s",
                    paste(names(smooth_blocks), collapse = " | ")))
    }
  }
  cols_s <- smooth_blocks[[label_s]]
  Xp_s   <- Xp_obs[, cols_s, drop = FALSE]
  bs_d   <- beta_draws[, cols_s, drop = FALSE]
  s_obs  <- Xp_s %*% t(bs_d)
  
  # ---- variance components (mgcv reports point estimates) ----
  sig2 <- as.numeric(gam_fit$sig2)
  sp_named <- gam_fit$sp
  
  get_lambda <- function(label) {
    if (label %in% names(sp_named)) sp_named[[label]] else NA_real_
  }
  
  tau2_j_pt <- numeric(sim$p)
  for (j in seq_len(sim$p)) {
    tau2_j_pt[j] <- sig2 / get_lambda(sprintf("s(%s)", smooth_names[j]))
  }
  tau2_int_pt <- numeric(length(pair_keys))
  names(tau2_int_pt) <- pair_keys
  for (k in pair_keys) {
    parts <- strsplit(k, "_", fixed = TRUE)[[1]]
    u <- parts[1]; v <- parts[2]
    # mgcv assigns a separate smoothing param per marginal in ti();
    # report the average (or NA if any missing). This is approximate
    # and not directly comparable to our scalar tau2_int.
    label_pat <- sprintf("ti\\(X%s.*X%s\\)", u, v)
    matched <- grep(label_pat, names(sp_named), value = TRUE)
    if (length(matched) > 0) {
      tau2_int_pt[k] <- sig2 / mean(sp_named[matched])
    } else {
      tau2_int_pt[k] <- NA_real_
    }
  }
  
  lam_s <- get_lambda(label_s)
  tau2_s_pt <- sig2 / lam_s
  
  rep_n <- function(x) rep(x, S$n_draws)
  var_comp <- list(
    sigma2   = rep_n(sig2),
    tau2     = lapply(tau2_j_pt, rep_n),
    tau2_s   = rep_n(tau2_s_pt),
    tau2_int = lapply(tau2_int_pt, rep_n),
    rho      = NULL,
    nu       = NULL
  )
  
  post_sec <- (proc.time() - t_post_start)[3]
  total_sec <- (proc.time() - t0)[3]
  
  fit <- list(
    method       = "mgcv",
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
      rhat = NULL, ess = NULL,
      reml_score = as.numeric(gam_fit$gcv.ubre),
      converged  = isTRUE(gam_fit$converged)
    ),
    settings     = c(S, list(
      scenario = sim$scenario,
      formula  = deparse1(fmla),
      sig2_hat = sig2,
      sp_hat   = as.list(sp_named)
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
