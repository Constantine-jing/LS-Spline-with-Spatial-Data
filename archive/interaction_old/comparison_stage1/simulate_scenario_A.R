# ============================================================
# simulate_scenario_A.R
# 
# Scenario A — baseline cross-method comparison.
#   * smooth main effects (no interaction)
#   * modest spatial dependence
#   * n = 500, p = 3
# 
# This is the SIMULATION-ONLY DGP for the recovery comparison.
# The truth functions, knot count expectation, and spatial
# parameters are documented in:
#   comparison/output/settings/scenarioA_settings.md
# 
# Anchored on the validated centered-truth approach used in
# Sim1-Sim4+: each f_j is shifted to integrate to zero on
# [0,1] so ANOVA-type identifiability holds and per-smooth
# RMSE is on a like-for-like scale across methods.
# 
# Truth functions (centered on [0,1]):
#   f1(x) = 2 (sin(pi x) - 2/pi)
#   f2(x) = 1.5 (exp(x - 0.5) - (exp(0.5) - exp(-0.5)))
#   f3(x) = 0.7 (x^2 - 1/3)
# 
# Spatial:
#   Locations s_i ~ Uniform([0,1]^2)
#   b ~ N(0, tau2_s * Matern(rho_true, nu_true))
# 
# Default parameter calibration (locked in for Scenario A):
#   sigma2_true = 1.0
#   tau2_s_true = 1.0       (signal-to-noise ratio ~ 1)
#   rho_true    = 0.06      (Matern nu=1 effective range ~ 0.24)
#   nu_true     = 1.0       (matches what every method can fix)
# 
# Usage:
#   source("simulate_scenario_A.R")
#   sim <- simulate_scenario_A(seed = 1L)
#   sim$data           # data.frame with X1..X3, lon, lat, y
#   sim$truth_f_grid   # list of 3 numeric vectors on x_grid_1d
#   sim$truth_s_obs    # numeric vector length n (TRUE b at locs)
#   sim$truth_params   # list(sigma2, tau2_s, rho, nu)
#   sim$x_grid_1d      # numeric, length n_grid_1d
#   sim$x_grid_2d      # list(u, v); kept for API symmetry
#   sim$flat_grid_2d   # row-major flatten (n_grid_2d_total x 2)
#   sim$truth_f_int    # named list (empty for Scenario A)
# ============================================================

source("spatial_utils.R")  # for matern_cor() and pairdist()

# ------------------------------------------------------------
# Centered truth functions on [0,1].
# Hard-coded analytic constants so they don't drift between
# runs; verified to integrate to ~1e-5 on [0,1] (numerical
# limit of trapezoidal integration over a fine grid).
# ------------------------------------------------------------
.scenarioA_f1 <- function(x) 2.0  * (sin(pi * x) - 2 / pi)
.scenarioA_f2 <- function(x) 1.5  * (exp(x - 0.5) - (exp(0.5) - exp(-0.5)))
.scenarioA_f3 <- function(x) 0.7  * (x^2 - 1/3)

scenarioA_truth_f_list <- function() {
  list(.scenarioA_f1, .scenarioA_f2, .scenarioA_f3)
}

# ------------------------------------------------------------
# Main simulator.
# ------------------------------------------------------------
simulate_scenario_A <- function(seed         = 1L,
                                n            = 500L,
                                p            = 3L,
                                sigma2_true  = 1.0,
                                tau2_s_true  = 1.0,
                                rho_true     = 0.06,
                                nu_true      = 1.0,
                                n_grid_1d    = 101L,
                                n_grid_2d_per = 30L) {
  
  if (p != 3L) stop("Scenario A is fixed at p=3.")
  stopifnot(seed == as.integer(seed))
  
  set.seed(seed)
  
  # ---- covariates: iid Uniform(0,1) ----
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  
  # ---- spatial locations: iid Uniform([0,1]^2) ----
  loc <- cbind(lon = runif(n), lat = runif(n))
  
  # ---- spatial random effect: b ~ N(0, tau2_s * R(rho, nu)) ----
  D <- pairdist(loc)
  R <- matern_cor(D, rho = rho_true, nu = nu_true)
  # numerical jitter for cholesky stability
  R_chol <- chol(R + 1e-8 * diag(n))
  z <- rnorm(n)
  b <- sqrt(tau2_s_true) * as.numeric(crossprod(R_chol, z))
  
  # ---- signal: sum of centered main effects ----
  truth_f <- scenarioA_truth_f_list()
  eta_main <- numeric(n)
  for (j in seq_len(p)) {
    eta_main <- eta_main + truth_f[[j]](X[, j])
  }
  
  # ---- response: y = eta + b + eps ----
  eps <- rnorm(n, mean = 0, sd = sqrt(sigma2_true))
  y   <- eta_main + b + eps
  
  # ---- evaluation grids ----
  x_grid_1d <- seq(0, 1, length.out = n_grid_1d)
  
  # 2D grid kept for API symmetry; not used in Scenario A
  u_grid <- seq(0, 1, length.out = n_grid_2d_per)
  v_grid <- seq(0, 1, length.out = n_grid_2d_per)
  flat2d <- expand.grid(u = u_grid, v = v_grid, KEEP.OUT.ATTRS = FALSE)
  flat2d <- as.matrix(flat2d)
  
  # ---- truth on the grids ----
  truth_f_grid <- lapply(truth_f, function(f) f(x_grid_1d))
  
  # ---- assemble data.frame for downstream wrappers ----
  df <- data.frame(X)
  df$lon <- loc[, "lon"]
  df$lat <- loc[, "lat"]
  df$y   <- y
  
  list(
    scenario      = "A",
    seed          = as.integer(seed),
    n             = n,
    p             = p,
    
    data          = df,
    
    # truth on grids
    truth_f_grid  = truth_f_grid,           # list of p
    truth_s_obs   = b,                      # length n
    truth_f_int   = list(),                 # empty for Scenario A
    
    # truth parameters
    truth_params  = list(
      sigma2 = sigma2_true,
      tau2   = NA_real_,    # main-effect smoothing variances are
                            # not part of the DGP -- the truth is
                            # the function, not a smoothing penalty
      tau2_s = tau2_s_true,
      rho    = rho_true,
      nu     = nu_true
    ),
    
    # grids
    x_grid_1d     = x_grid_1d,
    x_grid_2d     = list(u = u_grid, v = v_grid),
    flat_grid_2d  = flat2d,
    
    # for traceability
    settings      = list(
      n            = n,
      p            = p,
      sigma2_true  = sigma2_true,
      tau2_s_true  = tau2_s_true,
      rho_true     = rho_true,
      nu_true      = nu_true,
      n_grid_1d    = n_grid_1d,
      n_grid_2d_per = n_grid_2d_per,
      seed         = as.integer(seed)
    )
  )
}


# ------------------------------------------------------------
# Sanity check (run once when sourcing).
# Confirms centering of truth functions on the eval grid.
# ------------------------------------------------------------
.scenarioA_centering_check <- function(tol = 1e-3) {
  x <- seq(0, 1, length.out = 10001)
  flist <- scenarioA_truth_f_list()
  ok <- TRUE
  for (j in seq_along(flist)) {
    m <- mean(flist[[j]](x))
    if (abs(m) > tol) {
      warning(sprintf("scenarioA truth f%d not centered: mean = %g", j, m))
      ok <- FALSE
    }
  }
  invisible(ok)
}
.scenarioA_centering_check()
