source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")

# ---- sample sizes ----
n_vec <- c(100, 200, 400, 1000)

M <- 6
seed <- 42
tau2 <- 0.15
nu <- 1.5
rho_init <- 0.2
lambda_init <- 0.1
rho_upper_mult <- 1.0  # Sim1 cap

# ---- simulate (Sim1) ----
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

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
  data.frame(x = coords[,1], y_coord = coords[,2], Y = y, eta = eta, eps = eps, X)
}

stopifnot(ls_tests())
dir.create("marginal_out", showWarnings = FALSE)

truth_f_list <- default_truth_f_list()
var_names <- paste0("X", 1:4)

for (n in n_vec) {
  
  dat <- simulate_sim1(n = n, tau2 = tau2, seed = seed)
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  # ---- Sim1 rho upper cap (MUST recompute for each n) ----
  D <- pairdist(coords)
  rho_upper <- rho_upper_mult * max(D)
  
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4), nu = nu,
    rho_init = rho_init, lambda_init = lambda_init,
    rho_upper = rho_upper,
    verbose = TRUE
  )
  
  imp_tab <- importance_table(obj, var_names = var_names)
  
  # --- consistency key for truth-centering ---
  obj$X_raw_for_marginal <- X_raw
  
  curves  <- marginal_curves(
    obj, x_grid = seq(0,1,length.out=101),
    truth_f_list = truth_f_list, clip = TRUE
  )
  err_tab <- curve_error_table(curves, var_names = var_names)
  
  write.csv(imp_tab, file = sprintf("marginal_out/sim1_importance_n%d.csv", n), row.names = FALSE)
  write.csv(err_tab, file = sprintf("marginal_out/sim1_curveerr_n%d.csv", n), row.names = FALSE)
  
  png(filename = sprintf("marginal_out/sim1_marginals_n%d.png", n), width = 1200, height = 900)
  par(mfrow = c(2,2), mar = c(4,4,2,1))
  plot_marginals(curves, var_names = var_names, show_rmse = TRUE)
  dev.off()
  
  cat("Saved to marginal_out/ for Sim1, n=", n, "\n")
}
