# ============================================================
# run_full_bayes_all.R
# Full Bayesian model: rho sampled + RW2 + M=20
# Tests on Sim1-4
#
# Run: Rscript run_full_bayes_all.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

M_C    <- 20
seed   <- 42
nu     <- 1.5
n      <- 400
x_grid <- seq(0, 1, length.out = 101)

truth_f_list_4 <- default_truth_f_list()
truth_f_list_6 <- c(truth_f_list_4,
                     list(function(x) rep(0, length(x)),
                          function(x) rep(0, length(x))))

dir.create("bayes_out", showWarnings = FALSE)

sigma2_true <- 0.8; rho_true <- 0.2; tau2_true <- 0.15
rho_X <- 0.10; nu_X <- 1.0

# ---- Simulation functions ----
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim1 <- function(n = 400, p = 4, mu = 0, tau2 = 0.15, seed = 42) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  X <- matrix(runif(n*p,0,1), n, p); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu = mu)
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+eps,
             eta=eta, b=rep(0,n), eps=eps, X)
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

simulate_sim3 <- function(n = 400, p = 4, mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  D <- pairdist(coords)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X <- sapply(1:p, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X <- as.matrix(X); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2*R_b + diag(jitter_b,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

simulate_sim4 <- function(n = 400, p = 6, mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n,0,1), y = runif(n,0,1))
  D <- pairdist(coords)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X56 <- cbind(runif(n), runif(n))
  X <- cbind(X14, X56); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X[,1:4,drop=FALSE], mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2*R_b + diag(jitter_b,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

stopifnot(ls_tests())


# ============================================================
# Generic runner for one simulation
# ============================================================
run_full_one <- function(sim_name, dat, coords, X_raw, y, b_true,
                          sigma2_t, tau2_t, rho_t,
                          truth_f_list, var_names,
                          n_iter = 6000, n_burn = 1500) {
  
  cat(sprintf("\n========== Full Bayes: %s, n=%d, p=%d, M=%d ==========\n",
              sim_name, length(y), ncol(X_raw), M_C))
  
  p_covs <- ncol(X_raw)
  n <- length(y)
  
  # --- Build LS basis ---
  des <- ls_additive_build(X_raw, M_vec = rep(M_C, p_covs))
  H <- cbind(1, des$W)
  D <- pairdist(coords)
  
  cat(sprintf("  H: %d x %d\n", nrow(H), ncol(H)))
  
  # --- Initialize via REML ---
  # Quick REML with M=6 for starting values
  obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                         M_vec = rep(6, p_covs), nu = nu,
                         rho_init = rho_t,
                         lambda_init = tau2_t / max(sigma2_t, 0.01),
                         verbose = FALSE)
  
  # GLS init for M=20 eta
  R_init <- matern_cor(D, rho = obj6$fit$rho, nu = nu)
  Sigma_init <- max(obj6$fit$sigma2, 0.01) * R_init + max(obj6$fit$tau2, 0.01) * diag(n)
  L_init <- chol(Sigma_init + diag(1e-8, n))
  y_w <- forwardsolve(t(L_init), y)
  H_w <- forwardsolve(t(L_init), H)
  eta_init <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H)),
                               crossprod(H_w, y_w)))
  
  init <- list(
    eta    = eta_init,
    sigma2 = max(obj6$fit$sigma2, 0.01),
    tau2   = max(obj6$fit$tau2, 0.01),
    rho    = obj6$fit$rho,
    tau2_s = rep(1.0, p_covs)
  )
  
  cat(sprintf("  REML init: rho=%.4f sigma2=%.4f tau2=%.4f\n",
              init$rho, init$sigma2, init$tau2))
  
  # --- Run full Gibbs ---
  t0 <- proc.time()
  gs <- gibbs_full_sampler(
    y = y, H = H, D = D, nu = nu, col_map = des$col_map,
    n_iter = n_iter, n_burn = n_burn, n_thin = 1,
    kappa2 = 1e6,
    a_sigma = 2, b_sigma = 1,
    a_tau   = 2, b_tau   = 0.3,
    a_smooth = 1, b_smooth = 0.005,
    log_rho_mu = log(0.2),   # center prior near typical range
    log_rho_sd = 1.0,        # diffuse: covers ~[0.03, 1.5]
    mh_sd_log_sigma2 = 0.3,
    mh_sd_log_tau2   = 0.3,
    mh_sd_log_rho    = 0.2,
    init = init, verbose = TRUE
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  Done in %.1f seconds.\n", elapsed))
  
  # --- Summaries ---
  gsm <- gibbs_summary(gs)
  
  rho_summary <- c(mean = mean(gs$rho_samples), sd = sd(gs$rho_samples),
                    quantile(gs$rho_samples, c(0.025, 0.5, 0.975)))
  
  cat(sprintf("  sigma2: mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
              gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"], sigma2_t))
  cat(sprintf("  tau2:   mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
              gsm$tau2["mean"], gsm$tau2["2.5%"], gsm$tau2["97.5%"], tau2_t))
  cat(sprintf("  rho:    mean=%.4f  95%%CI=[%.4f, %.4f]  (true=%.4f)\n",
              rho_summary["mean"], rho_summary["2.5%"], rho_summary["97.5%"], rho_t))
  
  cat("  tau2_s: ")
  for (j in 1:p_covs) cat(sprintf("%s=%.4f ", var_names[j], mean(gs$tau2_s_samples[,j])))
  cat("\n")
  
  # --- Build obj structure ---
  obj <- list(des = des, X_fix = H,
              fit = list(beta = colMeans(gs$eta_samples)),
              X_raw_for_marginal = X_raw)
  
  # --- Curve coverage ---
  cov_tab <- gibbs_curve_coverage(gs, obj, x_grid, truth_f_list,
                                   X_raw = X_raw, var_names = var_names)
  cat("  Curve metrics:\n")
  print(cov_tab, digits = 4)
  
  # --- b recovery ---
  if (sigma2_t > 0) {
    b_post <- colMeans(gs$b_samples)
    cat(sprintf("  b: cor=%.4f  RMSE=%.4f\n", cor(b_true, b_post),
                sqrt(mean((b_true - b_post)^2))))
  }
  
  # --- Plots ---
  pdf(sprintf("bayes_out/full_%s.pdf", sim_name), width = 12, height = 9)
  par(mfrow = c(ceiling(p_covs/2), 2))
  plot_gibbs_marginals(gs, obj, x_grid = x_grid,
                        truth_f_list = truth_f_list,
                        var_names = var_names, X_raw = X_raw)
  dev.off()
  
  pdf(sprintf("bayes_out/full_%s_trace.pdf", sim_name), width = 10, height = 10)
  plot_gibbs_trace_full(gs, true_sigma2 = sigma2_t, true_tau2 = tau2_t,
                         true_rho = rho_t)
  dev.off()
  
  write.csv(cov_tab, sprintf("bayes_out/full_%s_coverage.csv", sim_name),
            row.names = FALSE)
  
  list(gs = gs, gsm = gsm, rho_summary = rho_summary,
       obj = obj, cov_tab = cov_tab, elapsed = elapsed)
}


# ============================================================
# Run all 4 sims
# ============================================================

cat("\n####################################################\n")
cat("# FULL BAYESIAN MODEL (rho sampled)\n")
cat("####################################################\n")

# --- Sim1 ---
cat("\n# Sim1: no spatial b\n")
dat1 <- simulate_sim1(n = n, tau2 = tau2_true, seed = seed)
res1 <- run_full_one("sim1", dat1,
  coords = as.matrix(dat1[, c("x","y_coord")]),
  X_raw  = as.matrix(dat1[, paste0("X",1:4)]),
  y = dat1$Y, b_true = dat1$b,
  sigma2_t = 0, tau2_t = tau2_true, rho_t = rho_true,
  truth_f_list = truth_f_list_4, var_names = paste0("X",1:4))

# --- Sim2 ---
cat("\n# Sim2: spatial b, iid X\n")
dat2 <- simulate_sim2(n = n, sigma2 = sigma2_true, rho = rho_true,
                      nu = nu, tau2 = tau2_true, seed = seed)
res2 <- run_full_one("sim2", dat2,
  coords = as.matrix(dat2[, c("x","y_coord")]),
  X_raw  = as.matrix(dat2[, paste0("X",1:4)]),
  y = dat2$Y, b_true = dat2$b,
  sigma2_t = sigma2_true, tau2_t = tau2_true, rho_t = rho_true,
  truth_f_list = truth_f_list_4, var_names = paste0("X",1:4))

# --- Sim3 ---
cat("\n# Sim3: spatial b + spatial X\n")
dat3 <- simulate_sim3(n = n, rho_X = rho_X, nu_X = nu_X,
                      sigma2 = sigma2_true, rho = rho_true,
                      nu = nu, tau2 = tau2_true, seed = seed)
res3 <- run_full_one("sim3", dat3,
  coords = as.matrix(dat3[, c("x","y_coord")]),
  X_raw  = as.matrix(dat3[, paste0("X",1:4)]),
  y = dat3$Y, b_true = dat3$b,
  sigma2_t = sigma2_true, tau2_t = tau2_true, rho_t = rho_true,
  truth_f_list = truth_f_list_4, var_names = paste0("X",1:4))

# --- Sim4 ---
cat("\n# Sim4: spatial confounding + garbage X5,X6\n")
dat4 <- simulate_sim4(n = n, rho_X = rho_X, nu_X = nu_X,
                      sigma2 = sigma2_true, rho = rho_true,
                      nu = nu, tau2 = tau2_true, seed = seed)
res4 <- run_full_one("sim4", dat4,
  coords = as.matrix(dat4[, c("x","y_coord")]),
  X_raw  = as.matrix(dat4[, paste0("X",1:6)]),
  y = dat4$Y, b_true = dat4$b,
  sigma2_t = sigma2_true, tau2_t = tau2_true, rho_t = rho_true,
  truth_f_list = truth_f_list_6, var_names = paste0("X",1:6))


# ============================================================
# Grand summary
# ============================================================
cat("\n\n####################################################\n")
cat("# GRAND SUMMARY: Full Bayesian model\n")
cat("####################################################\n\n")

for (nm in c("sim1", "sim2", "sim3", "sim4")) {
  res <- get(paste0("res", substr(nm, 4, 4)))
  sigma2_t <- if (nm == "sim1") 0 else sigma2_true
  
  cat(sprintf("--- %s ---\n", nm))
  print(res$cov_tab, digits = 4)
  cat(sprintf("  sigma2=%.4f (true=%.2f)  tau2=%.4f (true=%.2f)  rho=%.4f (true=%.2f)\n",
              res$gsm$sigma2["mean"], sigma2_t,
              res$gsm$tau2["mean"], tau2_true,
              res$rho_summary["mean"], rho_true))
  cat(sprintf("  MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              res$gs$accept_rate["sigma2"],
              res$gs$accept_rate["tau2"],
              res$gs$accept_rate["rho"]))
  cat("\n")
}

# Sim4 smoothing variances
cat("--- Sim4: Smoothing variances ---\n")
tau2s_means <- colMeans(res4$gs$tau2_s_samples)
names(tau2s_means) <- paste0("X", 1:6)
print(round(tau2s_means, 5))

cat(sprintf("\nTotal time: Sim1=%.0fs  Sim2=%.0fs  Sim3=%.0fs  Sim4=%.0fs\n",
            res1$elapsed, res2$elapsed, res3$elapsed, res4$elapsed))
cat("\nAll plots in bayes_out/full_sim*.pdf\n")
