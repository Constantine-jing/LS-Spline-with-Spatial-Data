# ============================================================
# Stage 1a: LS-only (no spatial b)
# ============================================================

# depends on ls_basis.R already sourced:
# - ls_additive_build()
# - ls_additive_design_new()
# - ls_tests()

# ------------------------------------------------------------
# 1) Additive truth 
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

# ------------------------------------------------------------
# 2) Simulator: coords + iid covariates + eta + iid noise
# ------------------------------------------------------------
simulate_ls_only <- function(n = 400,
                             domain = c(0,1,0,1),  # [xmin,xmax,ymin,ymax]
                             design = c("random","grid")[1],
                             p = 4,
                             mu = 0,
                             tau2 = 0.15,
                             seed = 42) {
  set.seed(seed)
  
  # coords kept for plug-and-play compatibility
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
  
  # iid covariates in [0,1]
  X <- matrix(runif(n*p, 0, 1), n, p)
  colnames(X) <- paste0("X", 1:p)
  
  # truth mean
  eta <- eta_truth_additive(X, mu = mu)
  
  # iid noise
  eps <- rnorm(n, mean = 0, sd = sqrt(tau2))
  
  # response
  y <- eta + eps
  
  data.frame(
    x = coords[,1],
    y_coord = coords[,2],
    Y = y,
    eta = eta,
    eps = eps,
    X
  )
}

# ------------------------------------------------------------
# 3) Fit: LS additive design + OLS
# ------------------------------------------------------------
fit_ls_only <- function(y, X_raw, M_vec) {
  y <- as.numeric(y)
  X_raw <- as.matrix(X_raw)
  
  des <- ls_additive_build(X_raw, M_vec = M_vec)
  W <- des$W
  X_fix <- cbind(1, W)   # intercept + identified LS blocks
  
  # OLS: beta_hat = (X'X)^{-1} X'y
  XtX <- crossprod(X_fix)
  Xty <- crossprod(X_fix, y)
  beta_hat <- solve(XtX, Xty)
  
  y_hat <- as.numeric(X_fix %*% beta_hat)
  resid <- y - y_hat
  
  list(beta = beta_hat, des = des, X_fix = X_fix,
       y_hat = y_hat, resid = resid)
}

# ------------------------------------------------------------
# 4) Predict: build W_new using stored objs and predict mean
# ------------------------------------------------------------
predict_ls_only <- function(obj, X_new_raw, clip = TRUE) {
  X_new_raw <- as.matrix(X_new_raw)
  W_new <- ls_additive_design_new(X_new_raw, obj$des$objs, clip = clip)
  X_new_fix <- cbind(1, W_new)
  as.numeric(X_new_fix %*% obj$beta)
}

# ------------------------------------------------------------
# 5) stage1a simulation
# ------------------------------------------------------------
stage1_demo <- function(n = 400, M = 6, tau2 = 0.15, seed = 42) {
  stopifnot(ls_tests())
  
  dat <- simulate_ls_only(n = n, p = 4, tau2 = tau2, seed = seed)
  X_raw <- as.matrix(dat[, paste0("X", 1:4)])
  y <- dat$Y
  
  obj <- fit_ls_only(y, X_raw, M_vec = rep(M, 4))
  
  rmse_eta <- sqrt(mean((obj$y_hat - dat$eta)^2))
  cat("RMSE(eta_hat, eta_true) =", rmse_eta, "\n")
  cat("corr(y_hat, eta_true) =", cor(obj$y_hat, dat$eta), "\n")
  cat("sd(resid) should be ~ sqrt(tau2) =", sd(obj$resid), "\n")
  
  # optional quick plot if you want:
  plot(dat$eta, obj$y_hat); abline(0,1)
  plot(dat$eta, obj$y_hat, xlab="eta (true mean)", ylab="eta_hat (fitted mean)")
  
  invisible(list(dat = dat, obj = obj))
}

# IN CONSOLE START TEST
out <- stage1_demo(n=400, M=6, tau2=0.15, seed=42)

















