# ============================================================
# fit_spatial_reml.R
# Shared fitter: Profile-REML for (rho, lambda) + GLS beta
# plus BLUP b_hat = R * Sigma0^{-1} (y - X beta)
#
# Depends on:
#   - spatial_utils.R: pairdist(), matern_cor(), solve_chol()
#   - ls_basis.R:      ls_additive_build()
# ============================================================

spatial_reml_fit <- function(y, X, coords, nu = 1.5,
                             rho_init = 0.2, lambda_init = 0.1,
                             rho_lower = 1e-4, rho_upper = NULL,
                             lambda_lower = 1e-8, lambda_upper = Inf,
                             jitter = 1e-8, verbose = TRUE) {
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  coords <- as.matrix(coords)
  
  n <- length(y)
  p <- ncol(X)
  stopifnot(nrow(X) == n)
  stopifnot(n > p)
  
  D <- pairdist(coords)
  maxD <- max(D)
  if (is.null(rho_upper)) rho_upper <- 10 * maxD   # default (Sim2/3)
  
  # REML objective (negative log-REML up to constant)
  reml_obj <- function(par) {
    rho    <- exp(par[1])
    lambda <- exp(par[2])
    
    if (rho < rho_lower || rho > rho_upper) return(1e30)
    if (lambda < lambda_lower || lambda > lambda_upper) return(1e30)
    
    R <- matern_cor(D, rho = rho, nu = nu)
    Sigma0 <- R + diag(lambda, n)
    
    U0 <- tryCatch(chol(Sigma0 + diag(jitter, n)), error = function(e) NULL)
    if (is.null(U0)) return(1e30)
    
    # whiten: y_t = L^{-1} y where L = t(U0)
    y_t <- forwardsolve(t(U0), y)
    X_t <- forwardsolve(t(U0), X)
    
    XtX <- crossprod(X_t)
    UX  <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(UX)) return(1e30)
    
    beta_hat <- solve(XtX, crossprod(X_t, y_t))
    resid_t  <- y_t - X_t %*% beta_hat
    rss <- sum(resid_t^2)
    
    sigma2_hat <- rss / (n - p)
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) return(1e30)
    
    logdet_Sigma0 <- 2 * sum(log(diag(U0)))
    logdet_XtX    <- 2 * sum(log(diag(UX)))
    
    0.5 * ((n - p) * log(sigma2_hat) + logdet_Sigma0 + logdet_XtX)
  }
  
  par0 <- c(log(rho_init), log(lambda_init))
  opt <- optim(par0, reml_obj, method = "L-BFGS-B",
               lower = c(log(rho_lower), log(lambda_lower)),
               upper = c(log(rho_upper), log(lambda_upper)))
  
  rho_hat    <- exp(opt$par[1])
  lambda_hat <- exp(opt$par[2])
  
  # Recompute at optimum
  R_hat <- matern_cor(D, rho = rho_hat, nu = nu)
  Sigma0_hat <- R_hat + diag(lambda_hat, n)
  
  U0 <- chol(Sigma0_hat + diag(jitter, n))
  y_t <- forwardsolve(t(U0), y)
  X_t <- forwardsolve(t(U0), X)
  
  XtX <- crossprod(X_t)
  beta_hat <- solve(XtX, crossprod(X_t, y_t))
  resid_t  <- y_t - X_t %*% beta_hat
  rss <- sum(resid_t^2)
  
  sigma2_hat <- rss / (n - p)
  tau2_hat   <- lambda_hat * sigma2_hat
  
  if (verbose) {
    cat("REML status:", opt$convergence, "\n")
    cat("rho_hat    =", rho_hat, "\n")
    cat("lambda_hat =", lambda_hat, "\n")
    cat("sigma2_hat =", sigma2_hat, "\n")
    cat("tau2_hat   =", tau2_hat, "\n")
  }
  
  list(beta = as.vector(beta_hat),
       rho = rho_hat, nu = nu, lambda = lambda_hat,
       sigma2 = sigma2_hat, tau2 = tau2_hat,
       opt = opt, D = D, R = R_hat, Sigma0 = Sigma0_hat)
}


fit_ls_spatial <- function(y, X_raw, coords, M_vec, nu = 1.5,
                           rho_init = 0.2, lambda_init = 0.1,
                           rho_upper = NULL,
                           jitter = 1e-8, verbose = TRUE) {
  
  y <- as.numeric(y)
  X_raw <- as.matrix(X_raw)
  coords <- as.matrix(coords)
  
  des <- ls_additive_build(X_raw, M_vec = M_vec)
  W <- des$W
  X_fix <- cbind(1, W)  # intercept + identified spline blocks
  
  fit <- spatial_reml_fit(
    y = y, X = X_fix, coords = coords, nu = nu,
    rho_init = rho_init, lambda_init = lambda_init,
    rho_upper = rho_upper, # <<< pass through
    jitter = jitter, verbose = verbose)
  
  resid <- as.numeric(y - X_fix %*% fit$beta)
  
  # BLUP: b_hat = R * Sigma0^{-1} resid  (sigma2 cancels)
  U0 <- chol(fit$Sigma0 + diag(jitter, length(y)))
  alpha <- solve_chol(U0, resid)         # alpha = Sigma0^{-1} resid
  b_hat <- as.numeric(fit$R %*% alpha)
  y_hat <- as.numeric(X_fix %*% fit$beta + b_hat)
  
  list(des = des, fit = fit, coords = coords,
       X_fix = X_fix, resid = resid,
       alpha = alpha, b_hat = b_hat, y_hat = y_hat)
}
