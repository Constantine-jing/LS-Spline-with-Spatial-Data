# ============================================================
# Stage 2 (simple): spatial covariates + spatial residual GP
# ============================================================

simulate_stage2 <- function(n = 400,
                            domain = c(0,1,0,1),
                            design = c("random","grid")[1],
                            p = 4,
                            mu = 0,
                            # X(s) spatial structure (mild)
                            rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                            # residual spatial GP b(s)
                            sigma2 = 0.8, rho = 0.2, nu = 1.5,
                            # nugget
                            tau2 = 0.15,
                            seed = 42,
                            jitter_b = 1e-8) {
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
  
  # 2) spatial covariates X_j(s)
  D <- pairdist(coords)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  
  # draw p GP fields then map to (0,1) by ranks (keeps marginal ~Uniform)
  X <- sapply(1:p, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X <- as.matrix(X)
  colnames(X) <- paste0("X", 1:p)
  
  # 3) truth mean eta(s) = mu + sum f_j(X_j(s))
  eta <- eta_truth_additive(X, mu = mu)
  
  # 4) residual spatial GP b(s)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R_b + diag(jitter_b, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  
  # 5) nugget
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

# ============================================================
# Stage 2: ONE full run (simulate -> fit -> diagnostics)
# ============================================================

stage2_run_once <- function(n = 400, M = 6,
                            rho_X = 0.10, nu_X = 1.0,
                            sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                            seed = 42) {
  
  dat <- simulate_stage2(
    n = n, design = "random", p = 4, mu = 0,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  cat("\n====================\n")
  cat("Stage 2: SIMULATION sanity checks\n")
  cat("====================\n")
  cat("sd(eps) ~ sqrt(tau2):", sd(dat$eps), "vs", sqrt(tau2), "\n")
  cat("cor(eta, b):", cor(dat$eta, dat$b), " (will be > Stage 1b sometimes)\n")
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = TRUE
  )
  
  mu_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  b_hat  <- as.numeric(obj$b_hat)
  
  cat("\n====================\n")
  cat("Stage 2: FIT diagnostics\n")
  cat("====================\n")
  cat("rho_hat    =", obj$fit$rho,    " (true", rho,    ")\n")
  cat("sigma2_hat =", obj$fit$sigma2, " (true", sigma2, ")\n")
  cat("tau2_hat   =", obj$fit$tau2,   " (true", tau2,   ")\n")
  
  cat("\n--- separation ---\n")
  cat("cor(mu_hat, eta_true) =", cor(mu_hat, dat$eta), "\n")
  cat("cor(b_hat,  b_true)   =", cor(b_hat,  dat$b),   "\n")
  cat("cor(mu_hat, b_true)   =", cor(mu_hat, dat$b),
      "  # watch this (confounding)\n")
  
  invisible(list(dat = dat, obj = obj, mu_hat = mu_hat, b_hat = b_hat))
}

out2 <- stage2_run_once()

