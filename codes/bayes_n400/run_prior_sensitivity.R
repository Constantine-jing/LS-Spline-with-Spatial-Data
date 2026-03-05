# ============================================================
# run_prior_sensitivity.R
# Prior sensitivity: run Sim2 with 3 prior bundles
# to show results are robust to prior choice
#
# Bundle 1: Default (what we've been using)
# Bundle 2: More informative (tighter IG)
# Bundle 3: Less informative (wider IG)
#
# Run: Rscript run_prior_sensitivity.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

M_C <- 20; seed <- 42; nu <- 1.5; n <- 400
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

# --- Simulate data ONCE ---
dat <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                     nu = nu, tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat[, c("x","y_coord")])
X_raw  <- as.matrix(dat[, paste0("X",1:4)])
y      <- dat$Y
b_true <- dat$b

# Build basis
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

base_init <- list(eta = eta_init,
                  sigma2 = max(obj6$fit$sigma2, 0.01),
                  tau2 = max(obj6$fit$tau2, 0.01),
                  rho = obj6$fit$rho,
                  tau2_s = rep(1.0, 4))

# --- Define 3 prior bundles ---
bundles <- list(
  default = list(
    label = "Default: IG(2,1), IG(2,0.3), IG(1,0.005)",
    a_sigma = 2, b_sigma = 1,
    a_tau = 2, b_tau = 0.3,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = log(0.2), log_rho_sd = 1.0
  ),
  informative = list(
    label = "Informative: IG(3,2), IG(3,0.3), IG(2,0.01)",
    a_sigma = 3, b_sigma = 2,
    a_tau = 3, b_tau = 0.3,
    a_smooth = 2, b_smooth = 0.01,
    log_rho_mu = log(0.2), log_rho_sd = 0.5
  ),
  diffuse = list(
    label = "Diffuse: IG(1,0.5), IG(1,0.1), IG(0.5,0.001)",
    a_sigma = 1, b_sigma = 0.5,
    a_tau = 1, b_tau = 0.1,
    a_smooth = 0.5, b_smooth = 0.001,
    log_rho_mu = log(0.2), log_rho_sd = 1.5
  )
)

# --- Run each bundle ---
results <- list()

for (bname in names(bundles)) {
  bndl <- bundles[[bname]]
  cat(sprintf("\n####################################################\n"))
  cat(sprintf("# Prior bundle: %s\n", bndl$label))
  cat(sprintf("####################################################\n"))
  
  t0 <- proc.time()
  gs <- gibbs_full_sampler(
    y = y, H = H, D = D, nu = nu, col_map = des$col_map,
    n_iter = 6000, n_burn = 1500, n_thin = 1,
    kappa2 = 1e6,
    a_sigma = bndl$a_sigma, b_sigma = bndl$b_sigma,
    a_tau = bndl$a_tau, b_tau = bndl$b_tau,
    a_smooth = bndl$a_smooth, b_smooth = bndl$b_smooth,
    log_rho_mu = bndl$log_rho_mu, log_rho_sd = bndl$log_rho_sd,
    mh_sd_log_sigma2 = 0.3, mh_sd_log_tau2 = 0.3, mh_sd_log_rho = 0.2,
    init = base_init, verbose = TRUE
  )
  elapsed <- (proc.time() - t0)[3]
  
  gsm <- gibbs_summary(gs)
  rho_sum <- c(mean = mean(gs$rho_samples),
               quantile(gs$rho_samples, c(0.025, 0.975)))
  
  obj_tmp <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)
  
  cov_tab <- gibbs_curve_coverage(gs, obj_tmp, x_grid, truth_f_list,
                                   X_raw = X_raw, var_names = var_names)
  
  cat(sprintf("\n  sigma2=%.4f [%.4f, %.4f]  tau2=%.4f [%.4f, %.4f]  rho=%.4f [%.4f, %.4f]\n",
              gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"],
              gsm$tau2["mean"], gsm$tau2["2.5%"], gsm$tau2["97.5%"],
              rho_sum["mean"], rho_sum["2.5%"], rho_sum["97.5%"]))
  print(cov_tab, digits = 4)
  cat(sprintf("  Time: %.0f seconds\n", elapsed))
  
  results[[bname]] <- list(gsm = gsm, rho_sum = rho_sum,
                            cov_tab = cov_tab, elapsed = elapsed)
}

# --- Summary table ---
cat("\n\n####################################################\n")
cat("# Prior Sensitivity Summary (Sim2, n=400)\n")
cat("####################################################\n\n")

cat(sprintf("%-15s  %8s  %8s  %8s  %6s  %6s  %6s  %6s\n",
            "Bundle", "sigma2", "tau2", "rho", "cov_X1", "cov_X2", "cov_X3", "cov_X4"))
cat(paste(rep("-", 80), collapse=""), "\n")

for (bname in names(results)) {
  r <- results[[bname]]
  cat(sprintf("%-15s  %8.4f  %8.4f  %8.4f  %6.3f  %6.3f  %6.3f  %6.3f\n",
              bname,
              r$gsm$sigma2["mean"], r$gsm$tau2["mean"], r$rho_sum["mean"],
              r$cov_tab$band_coverage[1], r$cov_tab$band_coverage[2],
              r$cov_tab$band_coverage[3], r$cov_tab$band_coverage[4]))
}
cat(sprintf("\nTrue values:     %8.4f  %8.4f  %8.4f\n", sigma2_true, tau2_true, rho_true))
