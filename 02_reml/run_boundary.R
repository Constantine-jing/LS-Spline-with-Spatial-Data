# ============================================================
# run_boundary.R
#
# Main runner for boundary experiments on Hellbender.
# Sources all dependencies and runs the three experiments.
#
# Exp 1 (failed — SIM_FAIL), Exp 2 passed, Exp 3 passed
#
# Output goes to boundary_out/ directory:
#   - exp1_matrix_boundary.csv
#   - exp2_stat_boundary.csv
#   - exp3_signal_boundary.csv
#   - boundary_log.txt  (console output captured by SLURM)
# ============================================================

cat("=== run_boundary.R started ===\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Node:", Sys.info()["nodename"], "\n")
cat("R version:", R.version.string, "\n\n")

# ---- Source dependencies ----
source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("variable_selection.R")

library(mvtnorm)

# ---- Sanity check ----
stopifnot(ls_tests())
cat("ls_tests() passed.\n\n")

# ---- Output directory ----
dir.create("boundary_out", showWarnings = FALSE)

# ============================================================
# Shared tools (same as boundary_experiments.R)
# ============================================================

rmse <- function(a, b) sqrt(mean((a - b)^2))

# "Strong" truth (same as Sim5)
eta_strong <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 10)
  mu +
    2.0 * sin(pi * X[,1]) +
    1.5 * exp(X[,2] - 0.5) +
    0.7 * (X[,3]^2) +
    0.5 * sin(2 * pi * X[,4]) +
    1.0 * (X[,5] - 0.5)^3 +
    0.8 * cos(pi * X[,6]) +
    1.2 * log(X[,7] + 0.1) +
    0.6 * abs(X[,8] - 0.5) +
    0.9 * sin(3 * pi * X[,9]) +
    0.4 * (X[,10]^2 - X[,10])
}

# "Weak" truth: same shapes, amplitudes ~3x smaller
eta_weak <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 10)
  mu +
    0.6 * sin(pi * X[,1]) +
    0.5 * exp(X[,2] - 0.5) +
    0.25 * (X[,3]^2) +
    0.15 * sin(2 * pi * X[,4]) +
    0.3 * (X[,5] - 0.5)^3 +
    0.25 * cos(pi * X[,6]) +
    0.4 * log(X[,7] + 0.1) +
    0.2 * abs(X[,8] - 0.5) +
    0.3 * sin(3 * pi * X[,9]) +
    0.12 * (X[,10]^2 - X[,10])
}

# Generic data generator
simulate_generic <- function(n, p_real, p_garbage,
                             eta_fun = eta_strong,
                             rho_X = 0.10, nu_X = 1.0,
                             sigma2 = 0.8, rho = 0.2, nu = 1.5,
                             tau2 = 0.15, seed = 42,
                             jitter_X = 1e-8, jitter_b = 1e-8) {
  set.seed(seed)
  p <- p_real + p_garbage

  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D <- pairdist(coords)

  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)

  X_real <- sapply(1:p_real, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X_real <- as.matrix(X_real)

  X_garbage <- matrix(runif(n * p_garbage), n, p_garbage)

  X <- cbind(X_real, X_garbage)
  colnames(X) <- paste0("X", 1:p)

  eta <- eta_fun(X[, 1:p_real, drop = FALSE])

  R_b <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))

  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps

  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# Fit + Wald test + one-row summary
run_one_boundary <- function(n, p_real, p_garbage, M = 6,
                              eta_fun = eta_strong,
                              seed = 42, alpha = 0.05,
                              verbose_fit = FALSE) {
  p <- p_real + p_garbage
  k <- 1 + p * (M - 1)

  if (k >= n) {
    return(data.frame(
      n=n, p_real=p_real, p_garbage=p_garbage, p_total=p, k=k,
      df_resid=n-k, status="INFEASIBLE",
      power=NA, FPR=NA, n_missed=NA, missed="",
      rho_hat=NA, sigma2_hat=NA, tau2_hat=NA,
      stringsAsFactors=FALSE))
  }

  df_resid <- n - k
  status <- if (df_resid < 20) "LOW_DF" else "OK"

  dat <- tryCatch(
    simulate_generic(n=n, p_real=p_real, p_garbage=p_garbage,
                     eta_fun=eta_fun, seed=seed),
    error = function(e) NULL)
  if (is.null(dat)) {
    return(data.frame(
      n=n, p_real=p_real, p_garbage=p_garbage, p_total=p, k=k,
      df_resid=df_resid, status="SIM_FAIL",
      power=NA, FPR=NA, n_missed=NA, missed="",
      rho_hat=NA, sigma2_hat=NA, tau2_hat=NA,
      stringsAsFactors=FALSE))
  }

  coords <- as.matrix(dat[, c("x","y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y

  obj <- tryCatch(
    fit_ls_spatial(y=y, X_raw=X_raw, coords=coords,
                   M_vec=rep(M, p), nu=1.5,
                   rho_init=0.2, lambda_init=0.15/0.8,
                   verbose=verbose_fit),
    error = function(e) NULL)

  if (is.null(obj)) {
    return(data.frame(
      n=n, p_real=p_real, p_garbage=p_garbage, p_total=p, k=k,
      df_resid=df_resid, status="FIT_FAIL",
      power=NA, FPR=NA, n_missed=NA, missed="",
      rho_hat=NA, sigma2_hat=NA, tau2_hat=NA,
      stringsAsFactors=FALSE))
  }

  sel <- tryCatch(
    select_covariates(obj, var_names=paste0("X",1:p),
                       alpha=alpha, method="BH"),
    error = function(e) NULL)

  if (is.null(sel)) {
    return(data.frame(
      n=n, p_real=p_real, p_garbage=p_garbage, p_total=p, k=k,
      df_resid=df_resid, status="WALD_FAIL",
      power=NA, FPR=NA, n_missed=NA, missed="",
      rho_hat=obj$fit$rho, sigma2_hat=obj$fit$sigma2,
      tau2_hat=obj$fit$tau2,
      stringsAsFactors=FALSE))
  }

  sel$truth <- c(rep("real", p_real), rep("garbage", p_garbage))
  power  <- mean(sel$selected[sel$truth == "real"])
  fpr    <- if (p_garbage > 0) mean(sel$selected[sel$truth == "garbage"]) else 0
  missed <- sel$var[sel$truth == "real" & !sel$selected]

  data.frame(
    n=n, p_real=p_real, p_garbage=p_garbage, p_total=p, k=k,
    df_resid=df_resid, status=status,
    power=power, FPR=fpr,
    n_missed=length(missed),
    missed=paste(missed, collapse=","),
    rho_hat=obj$fit$rho, sigma2_hat=obj$fit$sigma2,
    tau2_hat=obj$fit$tau2,
    stringsAsFactors=FALSE)
}


# ############################################################
# EXPERIMENT 1: Matrix Boundary
# ############################################################

cat("\n============================================================\n")
cat("EXPERIMENT 1: Matrix Boundary\n")
cat("============================================================\n\n")

exp1 <- data.frame()
M <- 6

for (n in c(100, 200, 400)) {
  p_max <- floor((n - 1) / 5)
  p_real <- 5
  frac_vec <- c(0.3, 0.5, 0.7, 0.85, 0.95, 1.0)
  pg_vec <- unique(sort(c(
    pmax(0, round(frac_vec * p_max) - p_real),
    p_max - p_real + 1   # one past boundary
  )))

  for (pg in pg_vec) {
    ptot <- p_real + pg; k <- 1 + ptot * 5
    cat(sprintf("  n=%d p_garbage=%d (k=%d df=%d) ... ", n, pg, k, n-k))
    t0 <- proc.time()
    row <- run_one_boundary(n=n, p_real=p_real, p_garbage=pg,
                             eta_fun=eta_strong, seed=42)
    dt <- (proc.time() - t0)[3]
    cat(sprintf("[%s] power=%.2f FPR=%.2f (%.1fs)\n",
                row$status,
                ifelse(is.na(row$power), -1, row$power),
                ifelse(is.na(row$FPR), -1, row$FPR), dt))
    exp1 <- rbind(exp1, row)
  }
  cat("\n")
}

write.csv(exp1, "boundary_out/exp1_matrix_boundary.csv", row.names=FALSE)
cat("Saved boundary_out/exp1_matrix_boundary.csv\n\n")
print(exp1[, c("n","p_garbage","p_total","k","df_resid","status",
               "power","FPR","sigma2_hat","tau2_hat")],
      row.names=FALSE, digits=3)


# ############################################################
# EXPERIMENT 2: Statistical Boundary
# ############################################################

cat("\n\n============================================================\n")
cat("EXPERIMENT 2: Statistical Boundary (garbage sweep)\n")
cat("============================================================\n\n")

exp2 <- data.frame()

for (n in c(1000, 2000)) {
  p_real <- 10
  pg_max <- floor((n - 1) / 5) - p_real
  pg_candidates <- c(10, 20, 40, 60, 80, 100, 130, 160)
  pg_vec <- pg_candidates[pg_candidates <= pg_max]

  cat(sprintf("n=%d (max feasible p_garbage=%d)\n", n, pg_max))

  for (pg in pg_vec) {
    k <- 1 + (10 + pg) * 5
    cat(sprintf("  p_garbage=%3d (k=%d df=%d) ... ", pg, k, n-k))
    t0 <- proc.time()
    row <- run_one_boundary(n=n, p_real=10, p_garbage=pg,
                             eta_fun=eta_strong, seed=42)
    dt <- (proc.time() - t0)[3]
    cat(sprintf("[%s] power=%.2f FPR=%.3f (%.1fs)", row$status, row$power, row$FPR, dt))
    if (nchar(row$missed) > 0) cat(sprintf("  missed: %s", row$missed))
    cat("\n")
    exp2 <- rbind(exp2, row)
  }
  cat("\n")
}

write.csv(exp2, "boundary_out/exp2_stat_boundary.csv", row.names=FALSE)
cat("Saved boundary_out/exp2_stat_boundary.csv\n\n")
print(exp2[, c("n","p_garbage","p_total","df_resid","power","FPR","n_missed","missed")],
      row.names=FALSE, digits=3)


# ############################################################
# EXPERIMENT 3: Signal Strength Boundary
# ############################################################

cat("\n\n============================================================\n")
cat("EXPERIMENT 3: Signal Strength Boundary\n")
cat("============================================================\n\n")

exp3 <- data.frame()
n <- 1000
p_real <- 10
pg_vec <- c(10, 20, 40, 60, 80, 100, 130)

for (strength in c("strong", "weak")) {
  eta_fun <- if (strength == "strong") eta_strong else eta_weak
  cat(sprintf("--- %s signals (n=%d) ---\n", toupper(strength), n))

  for (pg in pg_vec) {
    k <- 1 + (10 + pg) * 5
    if (k >= n) { cat(sprintf("  p_garbage=%3d: SKIP (k=%d >= n)\n", pg, k)); next }

    cat(sprintf("  p_garbage=%3d (k=%d df=%d) ... ", pg, k, n-k))
    t0 <- proc.time()
    row <- run_one_boundary(n=n, p_real=10, p_garbage=pg,
                             eta_fun=eta_fun, seed=42)
    row$strength <- strength
    dt <- (proc.time() - t0)[3]
    cat(sprintf("[%s] power=%.2f FPR=%.3f (%.1fs)", row$status, row$power, row$FPR, dt))
    if (nchar(row$missed) > 0) cat(sprintf("  missed: %s", row$missed))
    cat("\n")
    exp3 <- rbind(exp3, row)
  }
  cat("\n")
}

write.csv(exp3, "boundary_out/exp3_signal_boundary.csv", row.names=FALSE)
cat("Saved boundary_out/exp3_signal_boundary.csv\n\n")

# Side-by-side
cat("--- Strong vs Weak comparison ---\n")
for (pg in pg_vec) {
  s <- exp3[exp3$strength == "strong" & exp3$p_garbage == pg, ]
  w <- exp3[exp3$strength == "weak"   & exp3$p_garbage == pg, ]
  if (nrow(s) == 0 || nrow(w) == 0) next
  cat(sprintf("  p_garbage=%3d:  strong power=%.2f FPR=%.3f  |  weak power=%.2f FPR=%.3f\n",
              pg, s$power, s$FPR, w$power, w$FPR))
}


# ############################################################
# Done
# ############################################################
cat("\n\n=== run_boundary.R finished ===\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
