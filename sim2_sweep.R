# ============================================================
# Sim2_sweep.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")

# --- truth mean (shared) ---
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

# --- Sim2 data generator ---
simulate_sim2 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
  set.seed(seed)
  
  # coords
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
  
  # iid covariates
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  eta <- eta_truth_additive(X, mu = mu)
  
  # spatial GP b(s)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R + diag(jitter, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps
  
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# --- one run => one summary row (good for tables) ---
sim2_run_once <- function(n, M = 6,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                          seed = 42, verbose_fit = FALSE) {
  
  stopifnot(ls_tests())
  
  dat <- simulate_sim2(n = n, p = 4, sigma2 = sigma2, rho = rho, nu = nu,
                       tau2 = tau2, seed = seed)
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = verbose_fit
  )
  
  mu_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  b_hat  <- as.numeric(obj$b_hat)
  
  data.frame(
    n = n, seed = seed, M = M,
    rho_true = rho, rho_hat = obj$fit$rho,
    sigma2_true = sigma2, sigma2_hat = obj$fit$sigma2,
    tau2_true = tau2, tau2_hat = obj$fit$tau2,
    rmse_eta = rmse(mu_hat, dat$eta),
    cor_eta  = cor(mu_hat, dat$eta),
    cor_b    = cor(b_hat, dat$b),
    cor_eta_btrue = cor(mu_hat, dat$b)
  )
}

# --- sweep driver ---
sim2_sweep <- function(n_vec = c(100, 400, 1000,10000),
                       M = 6,
                       sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                       seed = 42,
                       max_n_dense = 10000) {
  
  out <- list()
  for (n in n_vec) {
    if (n > max_n_dense) {
      message(sprintf("[Sim2] skip n=%d (dense GP O(n^3) too expensive).", n))
      next
    }
    message(sprintf("[Sim2] running n=%d ...", n))
    out[[as.character(n)]] <- sim2_run_once(
      n = n, M = M,
      sigma2 = sigma2, rho = rho, nu = nu, tau2 = tau2,
      seed = seed, verbose_fit = FALSE
    )
  }
  do.call(rbind, out)
}

# Example:
library(mvtnorm)
res2 <- sim2_sweep(n_vec = c(100, 400, 1000, 10000), seed = 42)
print(res2)

