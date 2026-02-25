# ------------------------------------------------------------
# Additive truth 
# ------------------------------------------------------------
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X)
  stopifnot(ncol(X) >= 4)
  
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  
  f1 <- 2 * sin(pi * x1)
  f2 <- 1.5 * exp(x2 - 0.5)
  f3 <- 0.7 * (x3^2)
  f4 <- 0.5 * sin(2 * pi * x4)
  
  mu + f1 + f2 + f3 + f4
}


spatial_reml_fit <- function(y, X, coords, nu = 1.5,
                             rho_init = 0.2, lambda_init = 0.1,
                             rho_lower = 1e-4, rho_upper = NULL,
                             lambda_lower = 1e-8, lambda_upper = Inf,
                             jitter = 1e-8, verbose = TRUE) {
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y); p <- ncol(X)
  stopifnot(nrow(X) == n)
  stopifnot(n > p)
  
  D <- pairdist(coords)
  if (is.null(rho_upper)) rho_upper <- 10 * max(D)
  
  reml_obj <- function(par) {
    rho    <- exp(par[1])
    lambda <- exp(par[2])
    if (rho < rho_lower || rho > rho_upper) return(1e30)
    if (lambda < lambda_lower || lambda > lambda_upper) return(1e30)
    
    R <- matern_cor(D, rho = rho, nu = nu)
    Sigma0 <- R + diag(lambda, n)
    
    L0 <- tryCatch(chol(Sigma0 + diag(jitter, n)), error = function(e) NULL)
    if (is.null(L0)) return(1e30)
    
    y_t <- forwardsolve(t(L0), y)
    X_t <- forwardsolve(t(L0), X)
    
    XtX <- crossprod(X_t)
    LX <- tryCatch(chol(XtX), error = function(e) NULL)
    if (is.null(LX)) return(1e30)
    
    beta_hat <- solve(XtX, crossprod(X_t, y_t))
    resid_t <- y_t - X_t %*% beta_hat
    rss <- sum(resid_t^2)
    
    sigma2_hat <- rss / (n - p)
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) return(1e30)
    
    logdet_Sigma0 <- 2 * sum(log(diag(L0)))
    logdet_XtX    <- 2 * sum(log(diag(LX)))
    
    0.5 * ((n - p) * log(sigma2_hat) + logdet_Sigma0 + logdet_XtX)
  }
  
  par0 <- c(log(rho_init), log(lambda_init))
  opt <- optim(par0, reml_obj, method = "L-BFGS-B",
               lower = c(log(rho_lower), log(lambda_lower)),
               upper = c(log(rho_upper), log(lambda_upper)))
  
  rho_hat    <- exp(opt$par[1])
  lambda_hat <- exp(opt$par[2])
  
  R_hat <- matern_cor(D, rho = rho_hat, nu = nu)
  Sigma0_hat <- R_hat + diag(lambda_hat, n)
  
  L0 <- chol(Sigma0_hat + diag(jitter, n))
  y_t <- forwardsolve(t(L0), y)
  X_t <- forwardsolve(t(L0), X)
  XtX <- crossprod(X_t)
  beta_hat <- solve(XtX, crossprod(X_t, y_t))
  resid_t <- y_t - X_t %*% beta_hat
  rss <- sum(resid_t^2)
  
  sigma2_hat <- rss / (n - p)
  tau2_hat <- lambda_hat * sigma2_hat
  
  if (verbose) {
    cat("REML status:", opt$convergence, "\n")
    cat("rho_hat   =", rho_hat, "\n")
    cat("lambda_hat=", lambda_hat, "\n")
    cat("sigma2_hat=", sigma2_hat, "\n")
    cat("tau2_hat  =", tau2_hat, "\n")
  }
  
  list(beta = as.vector(beta_hat),
       rho = rho_hat, nu = nu, lambda = lambda_hat,
       sigma2 = sigma2_hat, tau2 = tau2_hat,
       opt = opt, D = D, R = R_hat, Sigma0 = Sigma0_hat)
}




fit_ls_spatial <- function(y, X_raw, coords, M_vec, nu = 1.5,
                           rho_init = 0.2, lambda_init = 0.1,
                           jitter = 1e-8, verbose = TRUE) {
  
  des <- ls_additive_build(X_raw, M_vec = M_vec)
  W <- des$W
  X_fix <- cbind(1, W)  # intercept + identified spline blocks
  
  fit <- spatial_reml_fit(y = y, X = X_fix, coords = coords, nu = nu,
                          rho_init = rho_init, lambda_init = lambda_init,
                          jitter = jitter, verbose = verbose)
  
  resid <- as.numeric(y - X_fix %*% fit$beta)
  
  L0 <- chol(fit$Sigma0 + diag(jitter, length(y)))
  alpha <- solve_chol(L0, resid)         # alpha = Sigma0^{-1} resid
  b_hat <- as.numeric(fit$R %*% alpha)   # b_hat = R * alpha
  y_hat <- as.numeric(X_fix %*% fit$beta + b_hat)
  
  list(des = des, fit = fit, coords = coords,
       X_fix = X_fix, resid = resid, alpha = alpha,
       b_hat = b_hat, y_hat = y_hat)
}




# ============================================================
# Stage 1b: LS mean + spatial GP b(s) + nugget
# ============================================================

simulate_stage1b <- function(n = 400,
                             domain = c(0,1,0,1),
                             design = c("random","grid")[1],
                             p = 4,
                             mu = 0,
                             # spatial GP params
                             sigma2 = 0.8, rho = 0.2, nu = 1.5,
                             # nugget
                             tau2 = 0.15,
                             seed = 42,
                             jitter = 1e-8) {
  set.seed(seed)
  
  # 1) coords
  if (design == "grid") {
    m <- ceiling(sqrt(n)); n <- m*m
    gx <- seq(domain[1], domain[2], length.out = m)
    gy <- seq(domain[3], domain[4], length.out = m)
    coords <- as.matrix(expand.grid(gx, gy))
    colnames(coords) <- c("x","y")
  } else {
    coords <- cbind(
      x = runif(n, domain[1], domain[2]),
      y = runif(n, domain[3], domain[4])
    )
  }
  
  # 2) covariates: iid U(0,1)  (IMPORTANT: Stage 1b keeps X indep)
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  # 3) truth mean eta (reuse Stage 1a truth)
  # expects eta_truth_additive() is available (from stage1_ls_only.R or copy here)
  eta <- eta_truth_additive(X, mu = mu)
  
  # 4) spatial GP b(s)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R + diag(jitter, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  
  # 5) nugget eps
  eps <- rnorm(n, mean = 0, sd = sqrt(tau2))
  
  # 6) response
  y <- eta + b + eps
  
  data.frame(
    x = coords[,1],
    y_coord = coords[,2],
    Y = y,
    eta = eta,
    b = b,
    eps = eps,
    X
  )
}

# run simu 1b
# dat <- simulate_stage1b(n=400, sigma2=0.8, rho=0.2, nu=1.5, tau2=0.15, seed=42)
# c(var(dat$b), var(dat$eps), cor(dat$b, dat$eps), cor(dat$eta, dat$b))


# ============================================================
# Stage 1b: ONE full run (simulate -> fit -> print diagnostics)
# ============================================================

stage1b_run_once <- function(n = 400, M = 6,
                             sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                             seed = 42) {
  
  # --- simulate ---
  dat <- simulate_stage1b(
    n = n, design = "random", p = 4, mu = 0,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  # --- quick truth checks (simulation sanity) ---
  cat("\n====================\n")
  cat("Stage 1b: SIMULATION sanity checks\n")
  cat("====================\n")
  cat("target sqrt(tau2) =", sqrt(tau2), "\n")
  cat("sd(eps)           =", sd(dat$eps), "\n")
  cat("var(b)            =", var(dat$b), " (target marginal sigma2 =", sigma2, ")\n")
  cat("cor(b, eps)       =", cor(dat$b, dat$eps), "\n")
  cat("cor(eta, b)       =", cor(dat$eta, dat$b), " (should be near 0)\n")
  
  # --- fit LS + spatial (you already have these functions) ---
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = TRUE
  )
  
  # --- decompose fitted mean + fitted spatial effect ---

  
  eta_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  b_hat   <- as.numeric(obj$b_hat)
  
  rmse_eta <- sqrt(mean((eta_hat - dat$eta)^2))
  rmse_b   <- sqrt(mean((b_hat  - dat$b)^2))
  
  cat("RMSE(eta_hat, eta_true) =", rmse_eta, "\n")
  cat("RMSE(b_hat,  b_true)   =", rmse_b,   "\n")
  
  
  
  
  # --- fit checks (recovery + separation) ---
  cat("\n====================\n")
  cat("Stage 1b: FIT diagnostics\n")
  cat("====================\n")
  cat("rho_hat    =", obj$fit$rho,    " (true", rho,    ")\n")
  cat("sigma2_hat =", obj$fit$sigma2, " (true", sigma2, ")\n")
  cat("tau2_hat   =", obj$fit$tau2,   " (true", tau2,   ")\n")
  
  cat("\n--- separation (this is the KEY part) ---\n")
  cat("cor(eta_hat, eta_true) =", cor(eta_hat, dat$eta), "\n")
  cat("cor(b_hat,  b_true)   =", cor(b_hat,  dat$b),   "\n")
  cat("cor(eta_hat, b_true)   =", cor(eta_hat, dat$b),
      "  # should be near 0 if no confounding\n")
  
  invisible(list(dat = dat, obj = obj, eta_hat = eta_hat, b_hat = b_hat))
}

source("ls_basis.R")
source("spatial_utils.R")
out <- stage1b_run_once(n=400, M=6, sigma2=0.8, rho=0.2, nu=1.5, tau2=0.15, seed=42)





