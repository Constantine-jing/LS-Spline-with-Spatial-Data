# ============================================================
# Sim4_sweep.R
#   Sim3 setting + add garbage X5, X6
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")

library(mvtnorm)

# --- truth mean (same) ---
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

# ------------------------------------------------------------
# Sim4 generator:
#   X1..X4: spatial GP -> rank -> approx Unif(0,1)  (same as Sim3)
#   X5..X6: iid Unif(0,1) (garbage)
#   eta uses ONLY X1..X4
# ------------------------------------------------------------
simulate_sim4 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p = 6, mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
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
  
  D <- pairdist(coords)
  
  # X1..X4: spatial covariates via latent GP then rank->U(0,1)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X14 <- as.matrix(X14)
  
  # X5..X6: garbage iid Unif(0,1)
  X56 <- cbind(runif(n), runif(n))
  
  X <- cbind(X14, X56)
  colnames(X) <- paste0("X", 1:p)
  
  # eta uses only first 4
  eta <- eta_truth_additive(X[,1:4, drop = FALSE], mu = mu)
  
  # residual spatial GP b(s)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R_b + diag(jitter_b, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps
  
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# --- one run => one summary row (like Sim3) ---
sim4_run_once <- function(n, M = 6,
                          rho_X = 0.10, nu_X = 1.0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                          seed = 42, verbose_fit = FALSE) {
  
  stopifnot(ls_tests())
  
  dat <- simulate_sim4(
    n = n, p = 6, mu = 0,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:6)])
  y      <- dat$Y
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 6), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = verbose_fit
  )
  
  mu_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  b_hat  <- as.numeric(obj$b_hat)
  
  data.frame(
    n = n, seed = seed, M = M,
    rhoX_true = rho_X,
    rho_true = rho, rho_hat = obj$fit$rho,
    sigma2_true = sigma2, sigma2_hat = obj$fit$sigma2,
    tau2_true = tau2, tau2_hat = obj$fit$tau2,
    cor_etab_true = cor(dat$eta, dat$b),
    rmse_eta = rmse(mu_hat, dat$eta),
    cor_eta  = cor(mu_hat, dat$eta),
    cor_b    = cor(b_hat, dat$b),
    cor_eta_btrue = cor(mu_hat, dat$b)
  )
}

# --- sweep driver ---
sim4_sweep <- function(n_vec = c(100, 400, 1000, 10000),
                       M = 6,
                       rho_X = 0.10, nu_X = 1.0,
                       sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                       seed = 42,
                       max_n_dense = 10000) {
  
  out <- list()
  for (n in n_vec) {
    if (n > max_n_dense) {
      message(sprintf("[Sim4] skip n=%d (dense GP O(n^3) too expensive).", n))
      next
    }
    message(sprintf("[Sim4] running n=%d ...", n))
    out[[as.character(n)]] <- sim4_run_once(
      n = n, M = M,
      rho_X = rho_X, nu_X = nu_X,
      sigma2 = sigma2, rho = rho, nu = nu, tau2 = tau2,
      seed = seed, verbose_fit = FALSE
    )
  }
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Optional: one-off marginal plot for X1..X6 (expect X5,X6 flat)
# ------------------------------------------------------------
sim4_plot_marginals_once <- function(n = 1000, M = 6,
                                     rho_X = 0.10, nu_X = 1.0,
                                     sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                                     seed = 42,
                                     png_file = "figs/marginal_sim4.png") {
  
  stopifnot(ls_tests())
  
  dat <- simulate_sim4(
    n = n, p = 6, mu = 0,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:6)])
  y      <- dat$Y
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 6), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = FALSE
  )
  
  # for truth-centering inside marginal_curves (your existing convention)
  obj$X_raw_for_marginal <- X_raw
  
  truth_f_list <- list(
    function(x) 2*sin(pi*x),
    function(x) 1.5*exp(x - 0.5),
    function(x) 0.7*(x^2),
    function(x) 0.5*sin(2*pi*x),
    function(x) 0*x,  # X5 garbage
    function(x) 0*x   # X6 garbage
  )
  
  var_names <- paste0("X", 1:6)
  
  curves <- marginal_curves(
    obj,
    x_grid = seq(0, 1, length.out = 101),
    truth_f_list = truth_f_list,
    clip = TRUE
  )
  
  print(curve_error_table(curves, var_names = var_names))
  
  dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
  png(png_file, width = 1200, height = 800, res = 120)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(2,3))
  plot_marginals(curves, var_names = paste0("X",1:6), show_rmse = TRUE,
                 ylim_shared = TRUE, symmetric = TRUE)
  invisible(list(dat = dat, obj = obj, curves = curves))
}

# ------------------------------------------------------------
# Make 3 marginal plots for n = 100, 400, 1000
# ------------------------------------------------------------
sim4_plot_marginals_triplet <- function(n_vec = c(100, 400, 1000),
                                        M = 6,
                                        rho_X = 0.10, nu_X = 1.0,
                                        sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                                        seed = 42,
                                        out_dir = "figs") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  outs <- vector("list", length(n_vec))
  names(outs) <- as.character(n_vec)
  
  for (k in seq_along(n_vec)) {
    n <- n_vec[k]
    png_file <- file.path(out_dir, sprintf("marginal_sim4_n%d.png", n))
    
    message(sprintf("[Sim4] plotting marginals for n=%d -> %s", n, png_file))
    
    outs[[k]] <- sim4_plot_marginals_once(
      n = n, M = M,
      rho_X = rho_X, nu_X = nu_X,
      sigma2 = sigma2, rho = rho, nu = nu, tau2 = tau2,
      seed = seed,
      png_file = png_file
    )
  }
  
  invisible(outs)
}

# source("Sim4_sweep.R")
# sim4_plot_marginals_triplet(n_vec = c(100, 400, 1000), seed = 42)

# res4 <- sim4_sweep(n_vec = c(100, 400, 1000), seed = 42)
# print(res4)
# sim4_plot_marginals_once(n = 10000, seed = 42, png_file = "figs/marginal_sim4_n10000.png")





