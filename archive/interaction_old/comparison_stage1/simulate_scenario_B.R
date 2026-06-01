# ============================================================
# simulate_scenario_B.R
# 
# Scenario B -- interaction recovery comparison.
#   * smooth main effects (same as Scenario A)
#   * ONE non-zero pairwise interaction: f12(X1, X2)
#   * the other p(p-1)/2 - 1 pairwise interactions are zero
#   * modest spatial dependence (same as Scenario A)
#   * n = 500, p = 3, c_int = 1.5
# 
# DGP locked to match the validated 2x2 interaction experiment
# (`run_interaction_2x2_v2.R`): same centered f1, f2, f3, same
# centered f12, same nu=1, same spatial calibration as Scenario A.
# 
# Truth functions (centered on [0,1]):
#   f1(x) = 2 (sin(pi x) - 2/pi)
#   f2(x) = 1.5 (exp(x - 0.5) - C2)
#   f3(x) = 0.7 (x^2 - 1/3)
#   f12(x1,x2) = c_int * (sin(pi x1) - 2/pi) * (exp(x2 - 0.5) - C2)
# where C2 = exp(0.5) - exp(-0.5).
# 
# All four functions integrate to ~0 on [0,1] (or [0,1]^2 for f12),
# AND the row/column means of f12 are ~0 -- i.e. f12 is orthogonal
# to the marginals so ANOVA-type identifiability holds.
# 
# Spatial:
#   Locations s_i ~ Uniform([0,1]^2)
#   b ~ N(0, tau2_s * Matern(rho_true, nu_true))
#   sigma2_true = 1.0,  tau2_s_true = 1.0
#   rho_true    = 0.06, nu_true     = 1.0
# 
# Usage:
#   source("simulate_scenario_B.R")
#   sim <- simulate_scenario_B(seed = 1L)
# ============================================================

source("spatial_utils.R")

# centering constants
.C2_F2 <- exp(0.5) - exp(-0.5)

# ------------------------------------------------------------
# Centered truth functions (matched to run_interaction_2x2_v2.R)
# ------------------------------------------------------------
.scenarioB_f1  <- function(x) 2.0  * (sin(pi * x) - 2 / pi)
.scenarioB_f2  <- function(x) 1.5  * (exp(x - 0.5) - .C2_F2)
.scenarioB_f3  <- function(x) 0.7  * (x^2 - 1/3)

.scenarioB_f12 <- function(x1, x2, c_int = 1.5) {
  c_int * (sin(pi * x1) - 2 / pi) * (exp(x2 - 0.5) - .C2_F2)
}

scenarioB_truth_f_list <- function() {
  list(.scenarioB_f1, .scenarioB_f2, .scenarioB_f3)
}


# ------------------------------------------------------------
# Main simulator.
# ------------------------------------------------------------
simulate_scenario_B <- function(seed         = 1L,
                                n            = 500L,
                                p            = 3L,
                                c_int        = 1.5,
                                sigma2_true  = 1.0,
                                tau2_s_true  = 1.0,
                                rho_true     = 0.06,
                                nu_true      = 1.0,
                                n_grid_1d    = 101L,
                                n_grid_2d_per = 30L) {
  
  if (p != 3L) stop("Scenario B is fixed at p=3.")
  stopifnot(seed == as.integer(seed))
  
  set.seed(seed)
  
  # ---- covariates ----
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  
  # ---- locations ----
  loc <- cbind(lon = runif(n), lat = runif(n))
  
  # ---- spatial random effect ----
  D <- pairdist(loc)
  R <- matern_cor(D, rho = rho_true, nu = nu_true)
  R_chol <- chol(R + 1e-8 * diag(n))
  z <- rnorm(n)
  b <- sqrt(tau2_s_true) * as.numeric(crossprod(R_chol, z))
  
  # ---- main effects ----
  truth_f <- scenarioB_truth_f_list()
  eta_main <- numeric(n)
  for (j in seq_len(p)) {
    eta_main <- eta_main + truth_f[[j]](X[, j])
  }
  
  # ---- interaction signal: only X1 x X2 is non-zero ----
  f12_obs <- .scenarioB_f12(X[, 1], X[, 2], c_int = c_int)
  
  # ---- response ----
  eps <- rnorm(n, mean = 0, sd = sqrt(sigma2_true))
  y   <- eta_main + f12_obs + b + eps
  
  # ---- evaluation grids ----
  x_grid_1d <- seq(0, 1, length.out = n_grid_1d)
  
  u_grid <- seq(0, 1, length.out = n_grid_2d_per)
  v_grid <- seq(0, 1, length.out = n_grid_2d_per)
  flat2d <- as.matrix(expand.grid(u = u_grid, v = v_grid,
                                  KEEP.OUT.ATTRS = FALSE))
  
  # ---- truth on grids ----
  truth_f_grid <- lapply(truth_f, function(f) f(x_grid_1d))
  
  # 2D truth: f12(u, v) on the flattened (u,v) grid.
  # The other pairs (1,3) and (2,3) are identically zero on the grid.
  truth_f_int <- list()
  truth_f_int[["1_2"]] <- .scenarioB_f12(flat2d[, "u"], flat2d[, "v"],
                                          c_int = c_int)
  if (p >= 3L) {
    truth_f_int[["1_3"]] <- numeric(nrow(flat2d))    # all zeros
    truth_f_int[["2_3"]] <- numeric(nrow(flat2d))    # all zeros
  }
  
  # ---- assemble ----
  df <- data.frame(X)
  df$lon <- loc[, "lon"]
  df$lat <- loc[, "lat"]
  df$y   <- y
  
  list(
    scenario      = "B",
    seed          = as.integer(seed),
    n             = n,
    p             = p,
    
    data          = df,
    
    truth_f_grid  = truth_f_grid,
    truth_s_obs   = b,
    truth_f_int   = truth_f_int,
    truth_f12_obs = f12_obs,    # at observed (X1, X2) -- diagnostic only
    
    truth_params  = list(
      sigma2 = sigma2_true,
      tau2   = NA_real_,
      tau2_s = tau2_s_true,
      rho    = rho_true,
      nu     = nu_true,
      c_int  = c_int
    ),
    
    x_grid_1d     = x_grid_1d,
    x_grid_2d     = list(u = u_grid, v = v_grid),
    flat_grid_2d  = flat2d,
    
    int_keys      = names(truth_f_int),    # canonical f_int keys for this scenario
    
    settings      = list(
      n            = n,
      p            = p,
      c_int        = c_int,
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
# Self-check on source.
# ------------------------------------------------------------
.scenarioB_centering_check <- function(tol = 1e-3) {
  x <- seq(0, 1, length.out = 10001)
  flist <- scenarioB_truth_f_list()
  ok <- TRUE
  for (j in seq_along(flist)) {
    m <- mean(flist[[j]](x))
    if (abs(m) > tol) {
      warning(sprintf("scenarioB truth f%d not centered: mean = %g", j, m))
      ok <- FALSE
    }
  }
  
  # 2D centering: integrate f12 over [0,1]^2 and verify ANOVA orthogonality
  ng <- 401
  xg <- seq(0, 1, length.out = ng)
  G  <- expand.grid(u = xg, v = xg)
  Fv <- .scenarioB_f12(G$u, G$v, c_int = 1.5)
  Fmat <- matrix(Fv, ng, ng)
  if (abs(mean(Fmat)) > tol) {
    warning(sprintf("scenarioB f12 not globally centered: mean = %g",
                    mean(Fmat)))
    ok <- FALSE
  }
  if (max(abs(rowMeans(Fmat))) > 5e-3 ||
      max(abs(colMeans(Fmat))) > 5e-3) {
    warning("scenarioB f12 row/col means not ~0 (ANOVA identifiability)")
    ok <- FALSE
  }
  invisible(ok)
}
.scenarioB_centering_check()
