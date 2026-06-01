# ============================================================
# canonical_schema.R
# 
# Defines the canonical fit-result object that every method
# wrapper (ours, mgcv, INLA, BayesX) must return.
# 
# Downstream code (compute_recovery_metrics.R) is method-
# agnostic: it consumes only this schema.
# 
# All methods produce JOINT posterior draws of the smooth
# functions evaluated on a fixed grid. For Bayesian methods
# these are real MCMC draws; for mgcv they are mvrnorm samples
# from the empirical-Bayes posterior covariance Vp.
# 
# Why draws-only (no "summary" type):
#   - Pointwise quantiles, RMSE-of-mean, and 95% coverage are
#     computed identically for all methods.
#   - For mgcv this matches the analytic intervals exactly when
#     n_draws is large; the cost is ~1MB per fit and is trivial.
#   - Keeps compute_recovery_metrics.R fully method-agnostic.
# 
# All draws matrices have rows = grid points, cols = draws.
# ============================================================

# ------------------------------------------------------------
# Schema definition (documentation; not enforced as a class).
# 
# fit_result <- list(
#   # ---- identification ----
#   method        : character scalar; one of
#                   "ours" | "mgcv" | "inla" | "bayesx"
#   scenario      : character scalar; e.g. "A"
#   seed          : integer scalar; the simulation seed
# 
#   # ---- main effects ----
#   # List of length p. Element j is an (n_grid_1d x n_draws)
#   # matrix of posterior draws of f_j evaluated on x_grid_1d.
#   # Functions must be CENTERED to mean zero on [0,1] (matches
#   # the truth in simulate_scenario_A) so RMSE is on a like-for-
#   # like scale. The wrapper is responsible for centering.
#   f_main        : list of p matrices
# 
#   # ---- interaction surfaces ----
#   # Named list. Names are like "1_2" meaning f_{1,2}.
#   # Each element is an (n_grid_2d x n_draws) matrix of draws
#   # on the row-major flattening of x_grid_2d (see below).
#   # Empty list when the model has no interaction terms.
#   f_int         : named list of matrices  (possibly empty)
# 
#   # ---- spatial random effect at observed locations ----
#   # (n x n_draws) matrix of draws of s(loc_i).
#   # NULL if the method has no spatial term.
#   s_obs         : matrix or NULL
# 
#   # ---- variance components ----
#   # Each element is a numeric vector of length n_draws.
#   # tau2 is a list of length p (one per main-effect smooth).
#   # tau2_int is a named list keyed like f_int (possibly empty).
#   # rho and nu may be NULL if fixed (mgcv with bs="gp" fixes them).
#   var_comp      : list(
#                     sigma2   = numeric(n_draws),
#                     tau2     = list of p numeric(n_draws),
#                     tau2_s   = numeric(n_draws) or NULL,
#                     tau2_int = named list (possibly empty),
#                     rho      = numeric(n_draws) or NULL,
#                     nu       = numeric(n_draws) or NULL
#                   )
# 
#   # ---- evaluation grids (echoed for self-containment) ----
#   x_grid_1d     : numeric vector, length n_grid_1d, in [0,1]
#   x_grid_2d     : list(u = numeric, v = numeric); the 2D grid
#                   is the Cartesian product, flattened row-major
#                   (u varies fastest, then v) so it matches the
#                   flattening used by f_int.
# 
#   # ---- diagnostics ----
#   timing        : list(total_sec, fit_sec, post_sec)
#   convergence   : list(rhat = named numeric, ess = named numeric)
#                   NULL for non-MCMC methods (mgcv).
# 
#   # ---- traceability ----
#   settings      : list of all settings actually used
#                   (knot count, prior hyperparameters, n_draws,
#                    nu_fixed, software version, etc.)
# )
# ------------------------------------------------------------


# ------------------------------------------------------------
# Validator. Call this at the END of every wrapper's
# fit_<method>_wrapper() to fail fast on shape mistakes.
# 
# Args:
#   fit       : the fit_result list to validate
#   n_obs     : expected number of observations
#   p         : expected number of main-effect smooths
#   n_grid_1d : expected length of x_grid_1d
#   n_grid_2d : expected total length of flattened 2D grid
#               (i.e. length(x_grid_2d$u) * length(x_grid_2d$v))
#   n_draws   : expected number of posterior draws
#   has_spatial   : logical, whether s_obs should be present
#   int_keys      : character vector of expected keys for f_int
#                   (e.g. c() for Scenario A, c("1_2") for B)
# ------------------------------------------------------------
validate_fit_result <- function(fit,
                                n_obs, p,
                                n_grid_1d, n_grid_2d,
                                n_draws,
                                has_spatial = TRUE,
                                int_keys    = character(0)) {
  
  err <- function(msg) stop("[validate_fit_result] ", msg, call. = FALSE)
  
  # ---- top-level fields ----
  required <- c("method", "scenario", "seed",
                "f_main", "f_int", "s_obs",
                "var_comp",
                "x_grid_1d", "x_grid_2d",
                "timing", "convergence", "settings")
  missing_fields <- setdiff(required, names(fit))
  if (length(missing_fields) > 0L) {
    err(sprintf("missing top-level fields: %s",
                paste(missing_fields, collapse = ", ")))
  }
  
  if (!is.character(fit$method) || length(fit$method) != 1L) {
    err("method must be a character scalar")
  }
  if (!fit$method %in% c("ours", "mgcv", "inla", "bayesx")) {
    err(sprintf("method '%s' not in allowed set", fit$method))
  }
  
  # ---- f_main ----
  if (!is.list(fit$f_main) || length(fit$f_main) != p) {
    err(sprintf("f_main must be a list of length p=%d (got %d)",
                p, length(fit$f_main)))
  }
  for (j in seq_len(p)) {
    M <- fit$f_main[[j]]
    if (!is.matrix(M) || !is.numeric(M)) {
      err(sprintf("f_main[[%d]] must be a numeric matrix", j))
    }
    if (nrow(M) != n_grid_1d || ncol(M) != n_draws) {
      err(sprintf("f_main[[%d]] has dim (%d, %d); expected (%d, %d)",
                  j, nrow(M), ncol(M), n_grid_1d, n_draws))
    }
    if (any(!is.finite(M))) {
      err(sprintf("f_main[[%d]] contains non-finite values", j))
    }
  }
  
  # ---- f_int ----
  if (!is.list(fit$f_int)) err("f_int must be a list (possibly empty)")
  if (!setequal(names(fit$f_int), int_keys)) {
    err(sprintf("f_int keys mismatch: got {%s}, expected {%s}",
                paste(names(fit$f_int), collapse = ","),
                paste(int_keys,         collapse = ",")))
  }
  for (k in int_keys) {
    M <- fit$f_int[[k]]
    if (!is.matrix(M) || !is.numeric(M)) {
      err(sprintf("f_int[['%s']] must be a numeric matrix", k))
    }
    if (nrow(M) != n_grid_2d || ncol(M) != n_draws) {
      err(sprintf("f_int[['%s']] has dim (%d, %d); expected (%d, %d)",
                  k, nrow(M), ncol(M), n_grid_2d, n_draws))
    }
  }
  
  # ---- s_obs ----
  if (has_spatial) {
    if (is.null(fit$s_obs)) err("s_obs is NULL but has_spatial=TRUE")
    if (!is.matrix(fit$s_obs) || !is.numeric(fit$s_obs)) {
      err("s_obs must be a numeric matrix")
    }
    if (nrow(fit$s_obs) != n_obs || ncol(fit$s_obs) != n_draws) {
      err(sprintf("s_obs has dim (%d, %d); expected (%d, %d)",
                  nrow(fit$s_obs), ncol(fit$s_obs),
                  n_obs, n_draws))
    }
  } else {
    if (!is.null(fit$s_obs)) err("s_obs must be NULL when has_spatial=FALSE")
  }
  
  # ---- var_comp ----
  vc <- fit$var_comp
  if (!is.list(vc)) err("var_comp must be a list")
  
  must_have <- c("sigma2", "tau2", "tau2_s", "tau2_int", "rho", "nu")
  miss_vc <- setdiff(must_have, names(vc))
  if (length(miss_vc) > 0L) {
    err(sprintf("var_comp missing fields: %s",
                paste(miss_vc, collapse = ", ")))
  }
  
  check_draws_vec <- function(v, label) {
    if (is.null(v)) return(invisible())
    if (!is.numeric(v) || length(v) != n_draws) {
      err(sprintf("var_comp$%s must be numeric of length n_draws=%d (got length %d)",
                  label, n_draws, length(v)))
    }
  }
  check_draws_vec(vc$sigma2, "sigma2")
  check_draws_vec(vc$tau2_s, "tau2_s")
  check_draws_vec(vc$rho,    "rho")
  check_draws_vec(vc$nu,     "nu")
  
  if (!is.list(vc$tau2) || length(vc$tau2) != p) {
    err(sprintf("var_comp$tau2 must be a list of length p=%d", p))
  }
  for (j in seq_len(p)) check_draws_vec(vc$tau2[[j]], sprintf("tau2[[%d]]", j))
  
  if (!is.list(vc$tau2_int)) err("var_comp$tau2_int must be a list")
  if (!setequal(names(vc$tau2_int), int_keys)) {
    err("var_comp$tau2_int keys must match f_int keys")
  }
  for (k in int_keys) check_draws_vec(vc$tau2_int[[k]], sprintf("tau2_int[['%s']]", k))
  
  # ---- grids ----
  if (!is.numeric(fit$x_grid_1d) || length(fit$x_grid_1d) != n_grid_1d) {
    err(sprintf("x_grid_1d must be numeric of length %d", n_grid_1d))
  }
  if (!is.list(fit$x_grid_2d) ||
      !all(c("u", "v") %in% names(fit$x_grid_2d))) {
    err("x_grid_2d must be a list with elements u and v")
  }
  if (length(fit$x_grid_2d$u) * length(fit$x_grid_2d$v) != n_grid_2d) {
    err("x_grid_2d$u %x% x_grid_2d$v size does not match n_grid_2d")
  }
  
  # ---- timing & settings ----
  if (!is.list(fit$timing) ||
      !all(c("total_sec", "fit_sec", "post_sec") %in% names(fit$timing))) {
    err("timing must be a list with total_sec, fit_sec, post_sec")
  }
  if (!is.list(fit$settings) || length(fit$settings) == 0L) {
    err("settings must be a non-empty list")
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------
# Helper: row-major flatten of a 2D grid. Used by every wrapper
# that produces interaction surfaces, and by the truth function
# in simulate_scenario_*. Keeping this in ONE place ensures
# every method orders the grid identically.
# 
# Returns a matrix with columns u, v and
#   length(u_grid) * length(v_grid) rows,
# with u varying fastest (i.e. for each value of v, all values
# of u).  Equivalent to expand.grid(u=u_grid, v=v_grid) which
# also has u varying fastest -- we wrap it just to lock the
# convention down.
# ------------------------------------------------------------
flatten_grid_2d <- function(u_grid, v_grid) {
  G <- expand.grid(u = u_grid, v = v_grid, KEEP.OUT.ATTRS = FALSE)
  as.matrix(G)
}
