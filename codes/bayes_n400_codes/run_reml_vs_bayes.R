# ============================================================
# run_reml_vs_bayes.R
# Head-to-head comparison: REML (M=6) vs Full Bayes (M=20+RW2)
# on Sim2 and Sim3
#
# Run: Rscript run_reml_vs_bayes.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

M_C <- 20; M_reml <- 6; seed <- 42; nu <- 1.5; n <- 400
x_grid <- seq(0, 1, length.out = 101)
truth_f_list <- default_truth_f_list()
var_names <- paste0("X", 1:4)
dir.create("bayes_out", showWarnings = FALSE)

sigma2_true <- 0.8; rho_true <- 0.2; tau2_true <- 0.15
rho_X <- 0.10; nu_X <- 1.0

eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim2 <- function(n=400, p=4, mu=0, sigma2=0.8, rho=0.2, nu=1.5,
                          tau2=0.15, seed=42, jitter=1e-8) {
  set.seed(seed)
  coords <- cbind(x=runif(n,0,1), y=runif(n,0,1))
  X <- matrix(runif(n*p,0,1), n, p); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu=mu)
  D <- pairdist(coords); R <- matern_cor(D, rho=rho, nu=nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma=sigma2*R+diag(jitter,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

simulate_sim3 <- function(n=400, p=4, mu=0, rho_X=0.10, nu_X=1.0, jitter_X=1e-8,
                          sigma2=0.8, rho=0.2, nu=1.5, tau2=0.15, seed=42, jitter_b=1e-8) {
  set.seed(seed)
  coords <- cbind(x=runif(n,0,1), y=runif(n,0,1))
  D <- pairdist(coords)
  R_X <- matern_cor(D, rho=rho_X, nu=nu_X); Sigma_X <- R_X+diag(jitter_X,n)
  X <- sapply(1:p, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma=Sigma_X))
    (rank(z, ties.method="average")-0.5)/n
  })
  X <- as.matrix(X); colnames(X) <- paste0("X",1:p)
  eta <- eta_truth_additive(X, mu=mu)
  R_b <- matern_cor(D, rho=rho, nu=nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma=sigma2*R_b+diag(jitter_b,n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x=coords[,1], y_coord=coords[,2], Y=eta+b+eps,
             eta=eta, b=b, eps=eps, X)
}

stopifnot(ls_tests())

# ============================================================
# Helper: run REML + compute curve RMSE
# ============================================================
run_reml <- function(y, X_raw, coords, M, truth_f_list, var_names) {
  obj <- fit_ls_spatial(y=y, X_raw=X_raw, coords=coords,
                        M_vec=rep(M, ncol(X_raw)), nu=nu,
                        rho_init=rho_true,
                        lambda_init=tau2_true/sigma2_true,
                        verbose=FALSE)
  obj$X_raw_for_marginal <- X_raw
  curves <- marginal_curves(obj, x_grid=x_grid, truth_f_list=truth_f_list, clip=TRUE)
  err <- curve_error_table(curves, var_names=var_names)
  list(obj=obj, rmse=err$rmse_curve)
}

# ============================================================
# Helper: run full Bayes + compute curve metrics
# ============================================================
run_bayes <- function(y, X_raw, coords, b_true, truth_f_list, var_names) {
  des <- ls_additive_build(X_raw, M_vec=rep(M_C, ncol(X_raw)))
  H <- cbind(1, des$W)
  D <- pairdist(coords)
  
  obj6 <- fit_ls_spatial(y=y, X_raw=X_raw, coords=coords,
                         M_vec=rep(6, ncol(X_raw)), nu=nu,
                         rho_init=rho_true,
                         lambda_init=tau2_true/sigma2_true,
                         verbose=FALSE)
  
  R_init <- matern_cor(D, rho=obj6$fit$rho, nu=nu)
  Sigma_init <- max(obj6$fit$sigma2, 0.01)*R_init + max(obj6$fit$tau2, 0.01)*diag(length(y))
  L_init <- chol(Sigma_init + diag(1e-8, length(y)))
  y_w <- forwardsolve(t(L_init), y)
  H_w <- forwardsolve(t(L_init), H)
  eta_init <- as.vector(solve(crossprod(H_w)+diag(1e-4, ncol(H)), crossprod(H_w, y_w)))
  
  init <- list(eta=eta_init, sigma2=max(obj6$fit$sigma2, 0.01),
               tau2=max(obj6$fit$tau2, 0.01), rho=obj6$fit$rho, tau2_s=rep(1.0, ncol(X_raw)))
  
  gs <- gibbs_full_sampler(
    y=y, H=H, D=D, nu=nu, col_map=des$col_map,
    n_iter=6000, n_burn=1500, n_thin=1,
    kappa2=1e6, a_sigma=2, b_sigma=1, a_tau=2, b_tau=0.3,
    a_smooth=1, b_smooth=0.005,
    log_rho_mu=log(0.2), log_rho_sd=1.0,
    mh_sd_log_sigma2=0.3, mh_sd_log_tau2=0.3, mh_sd_log_rho=0.2,
    init=init, verbose=TRUE
  )
  
  gsm <- gibbs_summary(gs)
  obj_tmp <- list(des=des, X_fix=H,
                  fit=list(beta=colMeans(gs$eta_samples)),
                  X_raw_for_marginal=X_raw)
  
  cov_tab <- gibbs_curve_coverage(gs, obj_tmp, x_grid, truth_f_list,
                                   X_raw=X_raw, var_names=var_names)
  b_post <- colMeans(gs$b_samples)
  
  list(gs=gs, gsm=gsm, cov_tab=cov_tab,
       rho_mean=mean(gs$rho_samples),
       b_cor=cor(b_true, b_post),
       rmse=cov_tab$rmse_curve)
}


# ============================================================
# Sim2
# ============================================================
cat("\n####################################################\n")
cat("# REML vs Full Bayes: Sim2\n")
cat("####################################################\n")

dat2 <- simulate_sim2(n=n, sigma2=sigma2_true, rho=rho_true, nu=nu, tau2=tau2_true, seed=seed)
coords2 <- as.matrix(dat2[, c("x","y_coord")])
X_raw2  <- as.matrix(dat2[, paste0("X",1:4)])

cat("\nRunning REML (M=6)...\n")
t0 <- proc.time()
reml2 <- run_reml(dat2$Y, X_raw2, coords2, M_reml, truth_f_list, var_names)
t_reml2 <- (proc.time()-t0)[3]

cat("\nRunning Full Bayes (M=20)...\n")
t0 <- proc.time()
bayes2 <- run_bayes(dat2$Y, X_raw2, coords2, dat2$b, truth_f_list, var_names)
t_bayes2 <- (proc.time()-t0)[3]


# ============================================================
# Sim3
# ============================================================
cat("\n####################################################\n")
cat("# REML vs Full Bayes: Sim3\n")
cat("####################################################\n")

dat3 <- simulate_sim3(n=n, rho_X=rho_X, nu_X=nu_X, sigma2=sigma2_true,
                      rho=rho_true, nu=nu, tau2=tau2_true, seed=seed)
coords3 <- as.matrix(dat3[, c("x","y_coord")])
X_raw3  <- as.matrix(dat3[, paste0("X",1:4)])

cat("\nRunning REML (M=6)...\n")
t0 <- proc.time()
reml3 <- run_reml(dat3$Y, X_raw3, coords3, M_reml, truth_f_list, var_names)
t_reml3 <- (proc.time()-t0)[3]

cat("\nRunning Full Bayes (M=20)...\n")
t0 <- proc.time()
bayes3 <- run_bayes(dat3$Y, X_raw3, coords3, dat3$b, truth_f_list, var_names)
t_bayes3 <- (proc.time()-t0)[3]


# ============================================================
# Comparison tables
# ============================================================
cat("\n\n####################################################\n")
cat("# REML vs Full Bayes: Comparison\n")
cat("####################################################\n\n")

cat("--- Sim2: Curve RMSE ---\n")
cat(sprintf("%-8s  %8s  %8s  %8s  %8s\n", "", "X1", "X2", "X3", "X4"))
cat(sprintf("%-8s  %8.4f  %8.4f  %8.4f  %8.4f\n", "REML",
            reml2$rmse[1], reml2$rmse[2], reml2$rmse[3], reml2$rmse[4]))
cat(sprintf("%-8s  %8.4f  %8.4f  %8.4f  %8.4f\n", "Bayes",
            bayes2$rmse[1], bayes2$rmse[2], bayes2$rmse[3], bayes2$rmse[4]))

cat("\n--- Sim2: Bayes coverage ---\n")
print(bayes2$cov_tab, digits = 4)

cat(sprintf("\n--- Sim2: Variance estimates ---\n"))
cat(sprintf("  REML:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml2$obj$fit$sigma2, reml2$obj$fit$tau2, reml2$obj$fit$rho))
cat(sprintf("  Bayes: sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            bayes2$gsm$sigma2["mean"], bayes2$gsm$tau2["mean"], bayes2$rho_mean))
cat(sprintf("  True:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            sigma2_true, tau2_true, rho_true))

cat(sprintf("\n--- Sim2: Time ---\n"))
cat(sprintf("  REML: %.1fs  Bayes: %.1fs\n", t_reml2, t_bayes2))

cat("\n\n--- Sim3: Curve RMSE ---\n")
cat(sprintf("%-8s  %8s  %8s  %8s  %8s\n", "", "X1", "X2", "X3", "X4"))
cat(sprintf("%-8s  %8.4f  %8.4f  %8.4f  %8.4f\n", "REML",
            reml3$rmse[1], reml3$rmse[2], reml3$rmse[3], reml3$rmse[4]))
cat(sprintf("%-8s  %8.4f  %8.4f  %8.4f  %8.4f\n", "Bayes",
            bayes3$rmse[1], bayes3$rmse[2], bayes3$rmse[3], bayes3$rmse[4]))

cat("\n--- Sim3: Bayes coverage ---\n")
print(bayes3$cov_tab, digits = 4)

cat(sprintf("\n--- Sim3: Variance estimates ---\n"))
cat(sprintf("  REML:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml3$obj$fit$sigma2, reml3$obj$fit$tau2, reml3$obj$fit$rho))
cat(sprintf("  Bayes: sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            bayes3$gsm$sigma2["mean"], bayes3$gsm$tau2["mean"], bayes3$rho_mean))
cat(sprintf("  True:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            sigma2_true, tau2_true, rho_true))

cat(sprintf("\n--- Sim3: Time ---\n"))
cat(sprintf("  REML: %.1fs  Bayes: %.1fs\n", t_reml3, t_bayes3))

# Save
comp_table <- data.frame(
  sim = rep(c("sim2","sim3"), each=2),
  method = rep(c("REML","Bayes"), 2),
  rmse_X1 = c(reml2$rmse[1], bayes2$rmse[1], reml3$rmse[1], bayes3$rmse[1]),
  rmse_X2 = c(reml2$rmse[2], bayes2$rmse[2], reml3$rmse[2], bayes3$rmse[2]),
  rmse_X3 = c(reml2$rmse[3], bayes2$rmse[3], reml3$rmse[3], bayes3$rmse[3]),
  rmse_X4 = c(reml2$rmse[4], bayes2$rmse[4], reml3$rmse[4], bayes3$rmse[4]),
  sigma2 = c(reml2$obj$fit$sigma2, bayes2$gsm$sigma2["mean"],
             reml3$obj$fit$sigma2, bayes3$gsm$sigma2["mean"]),
  tau2 = c(reml2$obj$fit$tau2, bayes2$gsm$tau2["mean"],
           reml3$obj$fit$tau2, bayes3$gsm$tau2["mean"]),
  rho = c(reml2$obj$fit$rho, bayes2$rho_mean,
          reml3$obj$fit$rho, bayes3$rho_mean)
)
write.csv(comp_table, "bayes_out/reml_vs_bayes.csv", row.names=FALSE)
cat("\nSaved to bayes_out/reml_vs_bayes.csv\n")
