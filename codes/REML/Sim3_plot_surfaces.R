# ============================================================
# Sim3_plot_surfaces.R
# Purpose: Sim3 (spatial X + spatial b) on a GRID and plot
#          2D surfaces: truth vs fitted vs errors (eta and b)
# Dependencies: source("ls_basis.R"), source("spatial_utils.R"),
#               source("fit_spatial_reml.R")
# ============================================================

# ---- load your project code ----
source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")

suppressPackageStartupMessages({
  library(mvtnorm)
})

# ============================================================
# Truth mean function (shared)
# ============================================================
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
  mu + 2*sin(pi*x1) + 1.5*exp(x2 - 0.5) + 0.7*(x3^2) + 0.5*sin(2*pi*x4)
}

# ============================================================
# Sim3 generator (copied here so this file is independent)
#   - spatial covariates X_j(s): latent GP then rank -> Unif(0,1)
#   - spatial residual b(s): Matérn GP
#   - y(s) = eta(s) + b(s) + eps
# ============================================================
simulate_sim3 <- function(n = 400, domain = c(0,1,0,1),
                          design = c("random","grid")[1],
                          p = 4, mu = 0,
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
  
  # spatial covariates via latent GP then rank->U(0,1)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  
  X <- sapply(1:p, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X <- as.matrix(X)
  colnames(X) <- paste0("X", 1:p)
  
  eta <- eta_truth_additive(X, mu = mu)
  
  # residual spatial GP b(s)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  Sigma_b <- sigma2 * R_b + diag(jitter_b, n)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  
  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps
  
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# ============================================================
# Grid reshape helpers (no dots; true surfaces via image())
# ============================================================
.grid_matrix <- function(dat, z) {
  gx <- sort(unique(dat$x))
  gy <- sort(unique(dat$y_coord))
  m  <- length(gx)
  stopifnot(length(gy) == m)
  
  # expand.grid(gx, gy): x varies fastest, then y
  ord <- order(dat$y_coord, dat$x)
  
  # byrow=TRUE matches "x changes fastest" within each y-row
  zmat <- matrix(z[ord], nrow = m, ncol = m, byrow = TRUE)
  list(gx = gx, gy = gy, zmat = zmat)
}

.plot_surface <- function(dat, z, main, add_contour = TRUE) {
  g <- .grid_matrix(dat, z)
  image(g$gx, g$gy, g$zmat, xlab = "x", ylab = "y", main = main, asp = 1)
  if (add_contour) contour(g$gx, g$gy, g$zmat, add = TRUE)
  invisible(g)
}

plot_sim3_recovery_panels <- function(dat, eta_hat, b_hat,
                                      add_contour = TRUE) {
  # Save ONLY what we change (this avoids "pin" issues)
  op <- par(mfrow = c(2,3),
            mar   = c(3,3,2,1),
            mgp   = c(2,0.7,0))
  on.exit(par(op), add = TRUE)
  
  .plot_surface(dat, dat$eta,           "True eta(s)", add_contour)
  .plot_surface(dat, eta_hat,           "Fitted eta_hat(s)", add_contour)
  .plot_surface(dat, eta_hat - dat$eta, "Error: eta_hat - eta", add_contour)
  
  .plot_surface(dat, dat$b,             "True b(s)", add_contour)
  .plot_surface(dat, b_hat,             "Fitted b_hat(s)", add_contour)
  .plot_surface(dat, b_hat - dat$b,     "Error: b_hat - b", add_contour)
}


# ============================================================
# Main: simulate (grid) -> fit -> plot surfaces
# ============================================================
run_sim3_surface_demo <- function(n = 400, M = 6,
                                  seed = 42,
                                  # X(s) params
                                  rho_X = 0.10, nu_X = 1.0,
                                  # b(s) params (truth + fit nu fixed)
                                  sigma2 = 0.8, rho = 0.2, nu = 1.5,
                                  # nugget
                                  tau2 = 0.15,
                                  verbose_fit = FALSE,
                                  add_contour = TRUE) {
  
  stopifnot(ls_tests())
  
  # --- simulate on a GRID ---
  dat <- simulate_sim3(
    n = n, design = "grid", seed = seed,
    rho_X = rho_X, nu_X = nu_X,
    sigma2 = sigma2, rho = rho, nu = nu,
    tau2 = tau2
  )
  
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:4)])
  y      <- dat$Y
  
  # --- fit (your existing REML pipeline) ---
  obj <- fit_ls_spatial(
    y = y, X_raw = X_raw, coords = coords,
    M_vec = rep(M, 4),
    nu = nu,
    rho_init = rho,
    lambda_init = tau2 / sigma2,
    verbose = verbose_fit
  )
  
  # fitted components at observed grid points
  eta_hat <- drop(obj$X_fix %*% obj$fit$beta)
  b_hat   <- drop(obj$b_hat)
  
  # --- plot ---
  plot_sim3_recovery_panels(dat, eta_hat, b_hat, add_contour = add_contour)
  
  invisible(list(dat = dat, obj = obj, eta_hat = eta_hat, b_hat = b_hat))
}

# ============================================================
# Run
# ============================================================
# You can change n to any square-ish number; code will use m=ceil(sqrt(n)) and then n=m^2.
out <- run_sim3_surface_demo(n = 400, M = 6, seed = 42)

# these now exist
dat     <- out$dat
eta_hat <- out$eta_hat
b_hat   <- out$b_hat

# plot again if you want
plot_sim3_recovery_panels(dat, eta_hat, b_hat)

save_sim3_surfaces_png <- function(filename = "sim3_surfaces.png",
                                   width = 1800, height = 1000, res = 200,
                                   dat, eta_hat, b_hat, add_contour = TRUE) {
  png(filename, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  plot_sim3_recovery_panels(dat, eta_hat, b_hat, add_contour = add_contour)
}
save_sim3_surfaces_png("sim3_surfaces.png", dat = dat, eta_hat = eta_hat, b_hat = b_hat)
