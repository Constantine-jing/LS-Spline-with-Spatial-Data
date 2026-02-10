# ============================================================
# Sim1.R  (no spatial effect in truth, but fit the same spatial model)
# Sweep over n and collect a results table
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")

# ------------------------------------------------------------
# 1) Additive truth
# ------------------------------------------------------------
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

# ------------------------------------------------------------
# 2) Simulator for Sim1 (b(s) ≡ 0)
# ------------------------------------------------------------
simulate_sim1 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p = 4, mu = 0, tau2 = 0.15, seed = 42) {
  set.seed(seed)
  
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
  
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  eta <- eta_truth_additive(X, mu = mu)
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + eps
  
  data.frame(
    x = coords[,1], y_coord = coords[,2],
    Y = y, eta = eta, eps = eps, X
  )
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

# ------------------------------------------------------------
# 3) One run: fit spatial model, return one-row summary
# ------------------------------------------------------------
sim1_run_once <- function(n = 400, M = 6, tau2 = 0.15, seed = 42,
                          nu = 1.5, rho_init = 0.2, lambda_init = 0.1,
                          rho_upper_mult = 1.0) {
  
  stopifnot(ls_tests())
  
  dat <- simulate_sim1(n = n, p = 4, tau2 = tau2, seed = seed)
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  # Sim1-friendly: prevent rho -> huge (flat R ~ all-ones)
  D <- pairdist(coords)
  maxD <- max(D)
  rho_upper <- rho_upper_mult * maxD
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho_init,
    lambda_init = lambda_init,
    rho_upper = rho_upper,
    verbose = FALSE
  )
  
  eta_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  
  data.frame(
    n = n,
    M = M,
    tau2_true = tau2,
    RMSE_eta = rmse(eta_hat, dat$eta),
    cor_eta  = cor(eta_hat, dat$eta),
    rho_hat  = obj$fit$rho,
    lambda_hat = obj$fit$lambda,
    sigma2_hat = obj$fit$sigma2,
    tau2_hat   = obj$fit$tau2,
    total_var  = obj$fit$sigma2 + obj$fit$tau2,
    maxD = maxD,
    rho_upper = rho_upper,
    seed = seed
  )
}

# ------------------------------------------------------------
# 4) Sweep over n
# ------------------------------------------------------------
sim1_sweep_n <- function(n_list = c(100, 400, 1000, 10000),
                         M = 6, tau2 = 0.15, seed = 42,
                         nu = 1.5, rho_init = 0.2, lambda_init = 0.1,
                         rho_upper_mult = 1.0) {
  
  res <- do.call(rbind, lapply(n_list, function(n) {
    sim1_run_once(
      n = n, M = M, tau2 = tau2, seed = seed,
      nu = nu, rho_init = rho_init, lambda_init = lambda_init,
      rho_upper_mult = rho_upper_mult
    )
  }))
  rownames(res) <- NULL
  res
}

# ------------------------------------------------------------
# Example
# ------------------------------------------------------------
out_tab <- sim1_sweep_n(n_list = c(100, 400, 1000, 10000))
print(out_tab)
#write.csv(out_tab, "sim1_sweep.csv", row.names = FALSE)


