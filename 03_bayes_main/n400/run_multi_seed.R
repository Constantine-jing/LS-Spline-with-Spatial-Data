# ============================================================
# run_multi_seed.R
# Multi-seed replication: run full Bayes on Sim2 with 5 seeds
# Shows results aren't seed-dependent
#
# Run: Rscript run_multi_seed.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

M_C <- 20; nu <- 1.5; n <- 400
x_grid <- seq(0, 1, length.out = 101)
truth_f_list <- default_truth_f_list()
var_names <- paste0("X", 1:4)
dir.create("bayes_out", showWarnings = FALSE)

sigma2_true <- 0.8; rho_true <- 0.2; tau2_true <- 0.15

eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim2 <- function(n = 400, p = 4, mu = 0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  X <- matrix(runif(n*p,0,1), n, p); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu = mu)
  D <- pairdist(coords)
  R <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2*R + diag(jitter,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

stopifnot(ls_tests())

seeds <- c(42, 123, 456, 789, 2025)

cat("\n####################################################\n")
cat("# Multi-seed replication: Sim2 full Bayes\n")
cat("####################################################\n")

all_results <- vector("list", length(seeds))

for (k in seq_along(seeds)) {
  s <- seeds[k]
  cat(sprintf("\n========== Seed %d (%d/%d) ==========\n", s, k, length(seeds)))
  
  dat <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                       nu = nu, tau2 = tau2_true, seed = s)
  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X",1:4)])
  y      <- dat$Y
  b_true <- dat$b
  
  des <- ls_additive_build(X_raw, M_vec = rep(M_C, 4))
  H <- cbind(1, des$W)
  D <- pairdist(coords)
  
  # REML init
  obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                         M_vec = rep(6, 4), nu = nu,
                         rho_init = rho_true,
                         lambda_init = tau2_true / sigma2_true,
                         verbose = FALSE)
  
  R_init <- matern_cor(D, rho = obj6$fit$rho, nu = nu)
  Sigma_init <- max(obj6$fit$sigma2, 0.01) * R_init + max(obj6$fit$tau2, 0.01) * diag(n)
  L_init <- chol(Sigma_init + diag(1e-8, n))
  y_w <- forwardsolve(t(L_init), y)
  H_w <- forwardsolve(t(L_init), H)
  eta_init <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H)),
                               crossprod(H_w, y_w)))
  
  init <- list(eta = eta_init,
               sigma2 = max(obj6$fit$sigma2, 0.01),
               tau2 = max(obj6$fit$tau2, 0.01),
               rho = obj6$fit$rho,
               tau2_s = rep(1.0, 4))
  
  t0 <- proc.time()
  gs <- gibbs_full_sampler(
    y = y, H = H, D = D, nu = nu, col_map = des$col_map,
    n_iter = 6000, n_burn = 1500, n_thin = 1,
    kappa2 = 1e6,
    a_sigma = 2, b_sigma = 1,
    a_tau = 2, b_tau = 0.3,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = log(0.2), log_rho_sd = 1.0,
    mh_sd_log_sigma2 = 0.3, mh_sd_log_tau2 = 0.3, mh_sd_log_rho = 0.2,
    init = init, verbose = TRUE
  )
  elapsed <- (proc.time() - t0)[3]
  
  gsm <- gibbs_summary(gs)
  rho_mean <- mean(gs$rho_samples)
  
  obj_tmp <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)
  
  cov_tab <- gibbs_curve_coverage(gs, obj_tmp, x_grid, truth_f_list,
                                   X_raw = X_raw, var_names = var_names)
  
  b_post <- colMeans(gs$b_samples)
  b_cor <- cor(b_true, b_post)
  
  all_results[[k]] <- list(
    seed = s, gsm = gsm, rho_mean = rho_mean,
    cov_tab = cov_tab, b_cor = b_cor, elapsed = elapsed
  )
  
  cat(sprintf("  sigma2=%.4f  tau2=%.4f  rho=%.4f  b_cor=%.4f  time=%.0fs\n",
              gsm$sigma2["mean"], gsm$tau2["mean"], rho_mean, b_cor, elapsed))
  print(cov_tab, digits = 4)
}

# --- Summary ---
cat("\n\n####################################################\n")
cat("# Multi-seed Summary (Sim2, n=400, full Bayes)\n")
cat("####################################################\n\n")

cat(sprintf("%-6s  %8s  %8s  %8s  %6s  %6s  %6s  %6s  %6s\n",
            "seed", "sigma2", "tau2", "rho", "cov_X1", "cov_X2", "cov_X3", "cov_X4", "b_cor"))
cat(paste(rep("-", 85), collapse=""), "\n")

for (r in all_results) {
  cat(sprintf("%-6d  %8.4f  %8.4f  %8.4f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n",
              r$seed,
              r$gsm$sigma2["mean"], r$gsm$tau2["mean"], r$rho_mean,
              r$cov_tab$band_coverage[1], r$cov_tab$band_coverage[2],
              r$cov_tab$band_coverage[3], r$cov_tab$band_coverage[4],
              r$b_cor))
}
cat(sprintf("\nTrue:   %8.4f  %8.4f  %8.4f\n", sigma2_true, tau2_true, rho_true))

# Averages
avg_sigma2 <- mean(sapply(all_results, function(r) r$gsm$sigma2["mean"]))
avg_tau2   <- mean(sapply(all_results, function(r) r$gsm$tau2["mean"]))
avg_rho    <- mean(sapply(all_results, function(r) r$rho_mean))
avg_cov    <- colMeans(do.call(rbind, lapply(all_results, function(r) r$cov_tab$band_coverage)))
cat(sprintf("Avg:    %8.4f  %8.4f  %8.4f  %6.3f  %6.3f  %6.3f  %6.3f\n",
            avg_sigma2, avg_tau2, avg_rho, avg_cov[1], avg_cov[2], avg_cov[3], avg_cov[4]))
