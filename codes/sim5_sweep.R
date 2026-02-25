# ============================================================
# Sim5_sweep.R
#   Sim3 setting + scale up: 10 real covariates, 10 garbage
#   X1..X10:  spatial GP -> rank -> approx Unif(0,1)
#   X11..X20: iid Unif(0,1) (garbage)
#   eta uses ONLY X1..X10
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")

library(mvtnorm)

# --- truth mean: 10-covariate additive ---
eta_truth_additive_10 <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 10)
  x1  <- X[,1];  x2  <- X[,2];  x3  <- X[,3];  x4  <- X[,4]
  x5  <- X[,5];  x6  <- X[,6];  x7  <- X[,7];  x8  <- X[,8]
  x9  <- X[,9];  x10 <- X[,10]
  mu +
    2.0 * sin(pi * x1) +
    1.5 * exp(x2 - 0.5) +
    0.7 * (x3^2) +
    0.5 * sin(2 * pi * x4) +
    1.0 * (x5 - 0.5)^3 +
    0.8 * cos(pi * x6) +
    1.2 * log(x7 + 0.1) +
    0.6 * abs(x8 - 0.5) +
    0.9 * sin(3 * pi * x9) +
    0.4 * (x10^2 - x10)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

# ------------------------------------------------------------
# Sim5 generator:
#   X1..X10:  spatial GP -> rank -> approx Unif(0,1) (same mechanism as Sim3)
#   X11..X20: iid Unif(0,1) (garbage)
#   eta uses ONLY X1..X10
# ------------------------------------------------------------
simulate_sim5 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p_real = 10, p_garbage = 10,
                          mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  p <- p_real + p_garbage

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

  # X1..X10: spatial covariates via latent GP then rank->U(0,1)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)

  X_real <- sapply(1:p_real, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X_real <- as.matrix(X_real)

  # X11..X20: garbage iid Unif(0,1)
  X_garbage <- matrix(runif(n * p_garbage), n, p_garbage)

  X <- cbind(X_real, X_garbage)
  colnames(X) <- paste0("X", 1:p)

  # eta uses only first 10
  eta <- eta_truth_additive_10(X[, 1:p_real, drop = FALSE], mu = mu)

  # residual spatial GP b(s)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R_b + diag(jitter_b, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))

  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps

  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# --- one run => one summary row ---
sim5_run_once <- function(n, M = 6,
                          p_real = 10, p_garbage = 10,
                          rho_X = 0.10, nu_X = 1.0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                          seed = 42, verbose_fit = FALSE) {

  stopifnot(ls_tests())
  p <- p_real + p_garbage

  dat <- simulate_sim5(
    n = n, p_real = p_real, p_garbage = p_garbage, mu = 0,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )

  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y

  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, p), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = verbose_fit
  )

  mu_hat <- as.numeric(obj$X_fix %*% obj$fit$beta)
  b_hat  <- as.numeric(obj$b_hat)

  data.frame(
    n = n, seed = seed, M = M,
    p_real = p_real, p_garbage = p_garbage,
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
sim5_sweep <- function(n_vec = c(100, 400, 1000, 10000),
                       M = 6,
                       p_real = 10, p_garbage = 10,
                       rho_X = 0.10, nu_X = 1.0,
                       sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                       seed = 42,
                       max_n_dense = 10000) {

  out <- list()
  for (n in n_vec) {
    if (n > max_n_dense) {
      message(sprintf("[Sim5] skip n=%d (dense GP O(n^3) too expensive).", n))
      next
    }
    message(sprintf("[Sim5] running n=%d ...", n))
    out[[as.character(n)]] <- sim5_run_once(
      n = n, M = M,
      p_real = p_real, p_garbage = p_garbage,
      rho_X = rho_X, nu_X = nu_X,
      sigma2 = sigma2, rho = rho, nu = nu, tau2 = tau2,
      seed = seed, verbose_fit = FALSE
    )
  }
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Truth function list for marginal plots (10 real + 10 garbage)
# ------------------------------------------------------------
sim5_truth_f_list <- function(p_real = 10, p_garbage = 10) {
  real_fns <- list(
    function(x) 2.0 * sin(pi * x),
    function(x) 1.5 * exp(x - 0.5),
    function(x) 0.7 * (x^2),
    function(x) 0.5 * sin(2 * pi * x),
    function(x) 1.0 * (x - 0.5)^3,
    function(x) 0.8 * cos(pi * x),
    function(x) 1.2 * log(x + 0.1),
    function(x) 0.6 * abs(x - 0.5),
    function(x) 0.9 * sin(3 * pi * x),
    function(x) 0.4 * (x^2 - x)
  )
  garbage_fns <- lapply(1:p_garbage, function(j) function(x) 0 * x)
  c(real_fns, garbage_fns)
}

# ------------------------------------------------------------
# Optional: one-off marginal plot for X1..X20
#   (expect X11..X20 flat)
# ------------------------------------------------------------
sim5_plot_marginals_once <- function(n = 1000, M = 6,
                                     p_real = 10, p_garbage = 10,
                                     rho_X = 0.10, nu_X = 1.0,
                                     sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                                     seed = 42,
                                     png_file = "figs/marginal_sim5.png") {

  stopifnot(ls_tests())
  p <- p_real + p_garbage

  dat <- simulate_sim5(
    n = n, p_real = p_real, p_garbage = p_garbage, mu = 0,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2, seed = seed
  )

  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y

  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, p), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = FALSE
  )

  # for truth-centering inside marginal_curves
  obj$X_raw_for_marginal <- X_raw

  truth_f_list <- sim5_truth_f_list(p_real, p_garbage)
  var_names <- paste0("X", 1:p)

  curves <- marginal_curves(
    obj,
    x_grid = seq(0, 1, length.out = 101),
    truth_f_list = truth_f_list,
    clip = TRUE
  )

  print(curve_error_table(curves, var_names = var_names))

  dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
  png(png_file, width = 1800, height = 1600, res = 120)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mfrow = c(4, 5))  # 20 covariates in a 4x5 grid
  plot_marginals(curves, var_names = var_names, show_rmse = TRUE,
                 ylim_shared = TRUE, symmetric = TRUE)
  invisible(list(dat = dat, obj = obj, curves = curves))
}

# ------------------------------------------------------------
# Make marginal plots for multiple n values
# ------------------------------------------------------------
sim5_plot_marginals_triplet <- function(n_vec = c(100, 400, 1000),
                                        M = 6,
                                        p_real = 10, p_garbage = 10,
                                        rho_X = 0.10, nu_X = 1.0,
                                        sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15,
                                        seed = 42,
                                        out_dir = "figs") {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  outs <- vector("list", length(n_vec))
  names(outs) <- as.character(n_vec)

  for (k in seq_along(n_vec)) {
    n <- n_vec[k]
    png_file <- file.path(out_dir, sprintf("marginal_sim5_n%d.png", n))

    message(sprintf("[Sim5] plotting marginals for n=%d -> %s", n, png_file))

    outs[[k]] <- sim5_plot_marginals_once(
      n = n, M = M,
      p_real = p_real, p_garbage = p_garbage,
      rho_X = rho_X, nu_X = nu_X,
      sigma2 = sigma2, rho = rho, nu = nu, tau2 = tau2,
      seed = seed,
      png_file = png_file
    )
  }

  invisible(outs)
}

# source("sim5_sweep.R")
# sim5_plot_marginals_triplet(n_vec = c(100, 400, 1000), seed = 42)

# res5 <- sim5_sweep(n_vec = c(100, 400, 1000), seed = 42)
# print(res5)
# sim5_plot_marginals_once(n = 10000, seed = 42, png_file = "figs/marginal_sim5_n10000.png")
