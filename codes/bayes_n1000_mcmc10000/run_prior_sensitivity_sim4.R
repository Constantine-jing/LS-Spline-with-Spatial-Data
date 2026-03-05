# ============================================================
# run_prior_sensitivity_sim4.R
#
# Prior sensitivity analysis for Sim4 (n=1000, p=6, 4 real + 2 garbage)
# Same data, 5 prior bundles, fixed seed=42
#
# Scientific goal: check whether conclusions about
#   (a) garbage covariate shrinkage (tau2_s for X5,X6)
#   (b) variance recovery (sigma2, tau2, rho)
#   (c) curve coverage for real covariates X1-X4
# are stable across plausible prior choices.
#
# Bundle design: vary ONE group at a time to isolate effects
#
#   B1 (Default)   : what we've been using throughout
#   B2 (Smooth-tight) : tighter smoothing prior -> stronger shrinkage on X5,X6
#   B3 (Smooth-loose) : looser smoothing prior -> weaker shrinkage, harder test
#   B4 (Var-diffuse)  : more diffuse IG on sigma2/tau2
#   B5 (Rho-misspec)  : rho prior centered away from truth (log-mean at 0.5)
#                       tests robustness to prior misspecification on range
#
# Run: Rscript run_prior_sensitivity_sim4.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

# ============================================================
# Settings (match run_sim4_n1000.R)
# ============================================================
M_C    <- 20
n      <- 1000
n_iter <- 10000
n_burn <- 2500
n_thin <- 1
seed   <- 42
nu     <- 1.5
p      <- 6

x_grid    <- seq(0, 1, length.out = 101)
var_names <- paste0("X", 1:6)

sigma2_true <- 0.8
rho_true    <- 0.2
tau2_true   <- 0.15
rho_X       <- 0.10
nu_X        <- 1.0

truth_f_list_6 <- c(
  default_truth_f_list(),
  list(function(x) rep(0, length(x)),
       function(x) rep(0, length(x)))
)

dir.create("bayes_out", showWarnings = FALSE)
stopifnot(ls_tests())

# ============================================================
# Simulate Sim4 ONCE (same data for all bundles)
# ============================================================
eta_truth_additive <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 4)
  mu + 2*sin(pi*X[,1]) + 1.5*exp(X[,2]-0.5) + 0.7*(X[,3]^2) + 0.5*sin(2*pi*X[,4])
}

simulate_sim4 <- function(n = 1000, p = 6, mu = 0,
                          rho_X = 0.10, nu_X = 1.0, jitter_X = 1e-8,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42, jitter_b = 1e-8) {
  set.seed(seed)
  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D      <- pairdist(coords)
  R_X    <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)
  X14 <- sapply(1:4, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X56 <- cbind(runif(n), runif(n))
  X   <- cbind(X14, X56); colnames(X) <- paste0("X", 1:p)
  eta <- eta_truth_additive(X[, 1:4, drop = FALSE], mu = mu)
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b   <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))
  eps <- rnorm(n, 0, sqrt(tau2))
  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = eta + b + eps, eta = eta, b = b, eps = eps, X)
}

cat(sprintf("\n=== Simulating Sim4: n=%d, p=%d (fixed seed=%d) ===\n", n, p, seed))
dat    <- simulate_sim4(n = n, p = p, rho_X = rho_X, nu_X = nu_X,
                        sigma2 = sigma2_true, rho = rho_true, nu = nu,
                        tau2 = tau2_true, seed = seed)
coords <- as.matrix(dat[, c("x", "y_coord")])
X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
y      <- dat$Y
b_true <- dat$b
D      <- pairdist(coords)

# ============================================================
# Build basis and shared REML initialisation
# ============================================================
des <- ls_additive_build(X_raw, M_vec = rep(M_C, p))
H   <- cbind(1, des$W)
cat(sprintf("  H: %d x %d\n", nrow(H), ncol(H)))

obj6 <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                       M_vec = rep(6, p), nu = nu,
                       rho_init = rho_true,
                       lambda_init = tau2_true / max(sigma2_true, 0.01),
                       verbose = FALSE)

R_init     <- matern_cor(D, rho = obj6$fit$rho, nu = nu)
Sigma_init <- max(obj6$fit$sigma2, 0.01) * R_init + max(obj6$fit$tau2, 0.01) * diag(n)
L_init     <- chol(Sigma_init + diag(1e-8, n))
y_w        <- forwardsolve(t(L_init), y)
H_w        <- forwardsolve(t(L_init), H)
eta_init   <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H)),
                               crossprod(H_w, y_w)))

base_init <- list(
  eta    = eta_init,
  sigma2 = max(obj6$fit$sigma2, 0.01),
  tau2   = max(obj6$fit$tau2,   0.01),
  rho    = obj6$fit$rho,
  tau2_s = rep(1.0, p)
)
cat(sprintf("  REML init (shared): rho=%.4f  sigma2=%.4f  tau2=%.4f\n",
            base_init$rho, base_init$sigma2, base_init$tau2))

# ============================================================
# Prior bundles
#
# Each parameter group and what it controls:
#
#  sigma2 ~ IG(a_sigma, b_sigma) : spatial variance
#    IG mean = b/(a-1), mode = b/(a+1)
#    Default IG(2,1): mean=1, mode=0.5   -- truth=0.8
#
#  tau2 ~ IG(a_tau, b_tau) : nugget variance
#    Default IG(2,0.3): mean=0.3, mode=0.15  -- truth=0.15
#
#  tau2_s_j ~ IG(a_smooth, b_smooth) : per-covariate smoothing var (RW2)
#    Default IG(1,0.005): very diffuse, allows both rough and flat curves
#    Tight   IG(2,0.01) : pushes toward small tau2_s -> MORE smoothing
#    Loose   IG(0.5,0.001): even more diffuse -> less regularisation
#
#  log(rho) ~ N(log_rho_mu, log_rho_sd^2)
#    Default: centered at log(0.2), sd=1.0  -- truth rho=0.2
#    Misspec: centered at log(0.5), sd=1.0  -- prior mean too large
# ============================================================

bundles <- list(

  # B1: Default -- unchanged from main analysis
  B1_default = list(
    label      = "B1 Default  [IG(2,1) | IG(2,0.3) | IG(1,0.005) | rho~logN(log0.2,1)]",
    a_sigma    = 2,    b_sigma  = 1,
    a_tau      = 2,    b_tau    = 0.3,
    a_smooth   = 1,    b_smooth = 0.005,
    log_rho_mu = log(0.2), log_rho_sd = 1.0
  ),

  # B2: Tighter smoothing prior -> stronger shrinkage on garbage covariates
  # IG(2, 0.01): mode = 0.01/3 ~ 0.003, mean = 0.01  (pulls tau2_s toward small)
  # Expect: X5/X6 tau2_s shrinks more, curves flatter
  B2_smooth_tight = list(
    label      = "B2 Smooth-tight  [IG(2,0.01) on tau2_s]",
    a_sigma    = 2,    b_sigma  = 1,
    a_tau      = 2,    b_tau    = 0.3,
    a_smooth   = 2,    b_smooth = 0.01,
    log_rho_mu = log(0.2), log_rho_sd = 1.0
  ),

  # B3: Looser smoothing prior -> weaker shrinkage, harder test for garbage detection
  # IG(0.5, 0.001): very heavy-tailed, allows large tau2_s freely
  # Expect: X5/X6 tau2_s less controlled, noisier garbage curve
  B3_smooth_loose = list(
    label      = "B3 Smooth-loose  [IG(0.5,0.001) on tau2_s]",
    a_sigma    = 2,    b_sigma  = 1,
    a_tau      = 2,    b_tau    = 0.3,
    a_smooth   = 0.5,  b_smooth = 0.001,
    log_rho_mu = log(0.2), log_rho_sd = 1.0
  ),

  # B4: Diffuse variance priors -- test robustness of sigma2/tau2 recovery
  # IG(1,0.5) for sigma2: mean=0.5, heavy tails
  # IG(1,0.1) for tau2: mean=0.1
  # Smoothing unchanged
  B4_var_diffuse = list(
    label      = "B4 Var-diffuse  [IG(1,0.5) sigma2 | IG(1,0.1) tau2]",
    a_sigma    = 1,    b_sigma  = 0.5,
    a_tau      = 1,    b_tau    = 0.1,
    a_smooth   = 1,    b_smooth = 0.005,
    log_rho_mu = log(0.2), log_rho_sd = 1.0
  ),

  # B5: Rho prior misspecified -- centered at 0.5 instead of 0.2
  # Tests: does the data (n=1000) correct a wrong prior on range?
  # Expect: posterior should pull back toward rho~0.2 with enough data
  B5_rho_misspec = list(
    label      = "B5 Rho-misspec  [log(rho)~N(log0.5,1), prior mean rho=0.5]",
    a_sigma    = 2,    b_sigma  = 1,
    a_tau      = 2,    b_tau    = 0.3,
    a_smooth   = 1,    b_smooth = 0.005,
    log_rho_mu = log(0.5), log_rho_sd = 1.0
  )
)

# ============================================================
# Run all bundles
# ============================================================
results <- list()

for (bname in names(bundles)) {
  bndl <- bundles[[bname]]
  cat(sprintf("\n####################################################\n"))
  cat(sprintf("# %s\n", bndl$label))
  cat(sprintf("####################################################\n"))

  t0 <- proc.time()
  gs <- gibbs_full_sampler(
    y        = y,
    H        = H,
    D        = D,
    nu       = nu,
    col_map  = des$col_map,
    n_iter   = n_iter,
    n_burn   = n_burn,
    n_thin   = n_thin,
    kappa2           = 1e6,
    a_sigma          = bndl$a_sigma,   b_sigma  = bndl$b_sigma,
    a_tau            = bndl$a_tau,     b_tau    = bndl$b_tau,
    a_smooth         = bndl$a_smooth,  b_smooth = bndl$b_smooth,
    log_rho_mu       = bndl$log_rho_mu,
    log_rho_sd       = bndl$log_rho_sd,
    mh_sd_log_sigma2 = 0.3,
    mh_sd_log_tau2   = 0.3,
    mh_sd_log_rho    = 0.2,
    init    = base_init,
    verbose = TRUE
  )
  elapsed <- (proc.time() - t0)[3]

  gsm         <- gibbs_summary(gs)
  rho_summary <- c(mean = mean(gs$rho_samples),
                   sd   = sd(gs$rho_samples),
                   quantile(gs$rho_samples, c(0.025, 0.5, 0.975)))
  tau2s_means <- colMeans(gs$tau2_s_samples)
  names(tau2s_means) <- var_names

  obj_tmp <- list(des = des, X_fix = H,
                  fit = list(beta = colMeans(gs$eta_samples)),
                  X_raw_for_marginal = X_raw)

  cov_tab <- gibbs_curve_coverage(gs, obj_tmp, x_grid, truth_f_list_6,
                                   X_raw = X_raw, var_names = var_names)

  b_post <- colMeans(gs$b_samples)
  b_cor  <- cor(b_true, b_post)

  cat(sprintf("  sigma2=%.4f [%.4f,%.4f]  true=%.2f\n",
              gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"], sigma2_true))
  cat(sprintf("  tau2  =%.4f [%.4f,%.4f]  true=%.2f\n",
              gsm$tau2["mean"],   gsm$tau2["2.5%"],   gsm$tau2["97.5%"],   tau2_true))
  cat(sprintf("  rho   =%.4f [%.4f,%.4f]  true=%.2f\n",
              rho_summary["mean"], rho_summary["2.5%"], rho_summary["97.5%"], rho_true))
  cat(sprintf("  b_cor =%.4f  time=%.0fs\n", b_cor, elapsed))
  cat(sprintf("  MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
              gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))
  cat("  tau2_s:\n")
  print(round(tau2s_means, 6))
  cat("  Curve coverage:\n")
  print(cov_tab, digits = 4)

  # Marginal curves plot
  pdf(sprintf("bayes_out/priorsens_%s_marginals.pdf", bname), width = 12, height = 9)
  par(mfrow = c(3, 2))
  plot_gibbs_marginals(gs, obj_tmp, x_grid = x_grid,
                       truth_f_list = truth_f_list_6,
                       var_names = var_names, X_raw = X_raw)
  dev.off()

  results[[bname]] <- list(
    gsm = gsm, rho_summary = rho_summary,
    tau2s_means = tau2s_means, cov_tab = cov_tab,
    b_cor = b_cor, elapsed = elapsed, gs = gs
  )
}

# ============================================================
# Grand summary table
# ============================================================
cat("\n\n####################################################\n")
cat("# PRIOR SENSITIVITY SUMMARY  (Sim4, n=1000)\n")
cat("####################################################\n\n")

# --- Variance recovery ---
cat("--- Variance recovery (posterior means) ---\n")
cat(sprintf("%-20s  %8s  %8s  %8s  %8s\n", "Bundle", "sigma2", "tau2", "rho", "b_cor"))
cat(paste(rep("-", 65), collapse=""), "\n")
for (bname in names(results)) {
  r <- results[[bname]]
  cat(sprintf("%-20s  %8.4f  %8.4f  %8.4f  %8.4f\n",
              bname,
              r$gsm$sigma2["mean"], r$gsm$tau2["mean"],
              r$rho_summary["mean"], r$b_cor))
}
cat(sprintf("%-20s  %8.4f  %8.4f  %8.4f\n", "TRUE",
            sigma2_true, tau2_true, rho_true))

# --- Smoothing variances (key: X5,X6 should be near 0) ---
cat("\n--- Smoothing variances tau2_s (mean over posterior) ---\n")
cat(sprintf("%-20s", "Bundle"))
for (v in var_names) cat(sprintf("  %8s", v))
cat("\n")
cat(paste(rep("-", 20 + 10*p), collapse=""), "\n")
for (bname in names(results)) {
  r <- results[[bname]]
  cat(sprintf("%-20s", bname))
  for (v in r$tau2s_means) cat(sprintf("  %8.5f", v))
  cat("\n")
}
cat(sprintf("  (X1-X4 = real, X5-X6 = garbage, truth f5=f6=0)\n")  )

# --- Curve RMSE ---
cat("\n--- Curve RMSE ---\n")
cat(sprintf("%-20s", "Bundle"))
for (v in var_names) cat(sprintf("  %6s", v))
cat("\n")
cat(paste(rep("-", 20 + 8*p), collapse=""), "\n")
for (bname in names(results)) {
  r <- results[[bname]]
  cat(sprintf("%-20s", bname))
  for (v in r$cov_tab$rmse_curve) cat(sprintf("  %6.4f", v))
  cat("\n")
}

# --- 95% band coverage ---
cat("\n--- 95% band coverage (nominal = 0.95) ---\n")
cat(sprintf("%-20s", "Bundle"))
for (v in var_names) cat(sprintf("  %6s", v))
cat("\n")
cat(paste(rep("-", 20 + 8*p), collapse=""), "\n")
for (bname in names(results)) {
  r <- results[[bname]]
  cat(sprintf("%-20s", bname))
  for (v in r$cov_tab$band_coverage) cat(sprintf("  %6.3f", v))
  cat("\n")
}

# --- Timing ---
cat("\n--- Timing ---\n")
for (bname in names(results)) {
  cat(sprintf("  %-20s  %.0f seconds\n", bname, results[[bname]]$elapsed))
}

# ============================================================
# Combined trace plot for rho across bundles (diagnostic)
# ============================================================
pdf("bayes_out/priorsens_rho_comparison.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))
for (bname in names(results)) {
  rho_chain <- results[[bname]]$gs$rho_samples
  plot(rho_chain, type = "l", col = "steelblue",
       main = bname, xlab = "Iteration (post-burnin)",
       ylab = expression(rho), ylim = c(0, 1))
  abline(h = rho_true, col = "red", lty = 2, lwd = 2)
  abline(h = mean(rho_chain), col = "darkgreen", lty = 1, lwd = 1.5)
  legend("topright", legend = c("chain", "true", "post.mean"),
         col = c("steelblue","red","darkgreen"), lty = c(1,2,1), bty = "n", cex = 0.8)
}
dev.off()
cat("\nRho comparison trace: bayes_out/priorsens_rho_comparison.pdf\n")

# ============================================================
# Save summary CSV
# ============================================================
summary_df <- do.call(rbind, lapply(names(results), function(bname) {
  r <- results[[bname]]
  data.frame(
    bundle   = bname,
    sigma2   = r$gsm$sigma2["mean"],
    tau2     = r$gsm$tau2["mean"],
    rho      = r$rho_summary["mean"],
    b_cor    = r$b_cor,
    tau2s_X1 = r$tau2s_means[1], tau2s_X2 = r$tau2s_means[2],
    tau2s_X3 = r$tau2s_means[3], tau2s_X4 = r$tau2s_means[4],
    tau2s_X5 = r$tau2s_means[5], tau2s_X6 = r$tau2s_means[6],
    rmse_X1  = r$cov_tab$rmse_curve[1], rmse_X2 = r$cov_tab$rmse_curve[2],
    rmse_X3  = r$cov_tab$rmse_curve[3], rmse_X4 = r$cov_tab$rmse_curve[4],
    rmse_X5  = r$cov_tab$rmse_curve[5], rmse_X6 = r$cov_tab$rmse_curve[6],
    cov_X1   = r$cov_tab$band_coverage[1], cov_X2 = r$cov_tab$band_coverage[2],
    cov_X3   = r$cov_tab$band_coverage[3], cov_X4 = r$cov_tab$band_coverage[4],
    cov_X5   = r$cov_tab$band_coverage[5], cov_X6 = r$cov_tab$band_coverage[6],
    time_sec = r$elapsed,
    row.names = NULL
  )
}))
write.csv(summary_df, "bayes_out/priorsens_sim4_n1000_summary.csv", row.names = FALSE)
cat("Saved: bayes_out/priorsens_sim4_n1000_summary.csv\n")

cat("\nAll outputs in bayes_out/:\n")
cat("  priorsens_B*_marginals.pdf  (one per bundle)\n")
cat("  priorsens_rho_comparison.pdf\n")
cat("  priorsens_sim4_n1000_summary.csv\n")
