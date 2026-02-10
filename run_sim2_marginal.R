source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
library(mvtnorm)

n <- 400
M <- 6
seed <- 42
sigma2 <- 0.8
rho <- 0.2
nu <- 1.5
tau2 <- 0.15

eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

simulate_sim2 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
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
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R + diag(jitter, n)))
  
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps
  
  data.frame(x = coords[,1], y_coord = coords[,2], Y = y, eta = eta, b = b, eps = eps, X)
}

stopifnot(ls_tests())

# ---- add this near the top (replace n <- 400) ----
n_vec <- c(100, 200, 400, 1000)

stopifnot(ls_tests())
dir.create("marginal_out", showWarnings = FALSE)

truth_f_list <- default_truth_f_list()
var_names <- paste0("X", 1:4)

for (n in n_vec) {
  
  dat <- simulate_sim2(n=n, sigma2=sigma2, rho=rho, nu=nu, tau2=tau2, seed=seed)
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = TRUE
  )
  
  imp_tab <- importance_table(obj, var_names = var_names)
  
  # --- consistency key: needed for truth-centering in marginal_curves() ---
  obj$X_raw_for_marginal <- X_raw
  
  curves  <- marginal_curves(
    obj, x_grid = seq(0,1,length.out=101),
    truth_f_list = truth_f_list, clip = TRUE
  )
  err_tab <- curve_error_table(curves, var_names = var_names)
  
  write.csv(imp_tab, sprintf("marginal_out/sim2_importance_n%d.csv", n), row.names = FALSE)
  write.csv(err_tab, sprintf("marginal_out/sim2_curveerr_n%d.csv", n), row.names = FALSE)
  
  png(filename = sprintf("marginal_out/sim2_marginals_n%d.png", n), width = 1200, height = 900)
  par(mfrow = c(2,2), mar = c(4,4,2,1))   # safer margins
  plot_marginals(curves, var_names = var_names, show_rmse = TRUE)
  dev.off()
  
  cat("Saved to marginal_out/ for Sim2, n=", n, "\n")
}
