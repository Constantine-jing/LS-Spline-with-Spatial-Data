# ============================================================
# run_boundary_exp1.R
#
# EXPERIMENT 1 ONLY: Matrix Boundary
# Fixed: uses a 5-covariate truth function (not 10)
#
# Tests: fix p_real=5, push p_garbage toward k=n boundary
#        n in {100, 200, 400}
#
# Usage:  Rscript run_boundary_exp1.R
# Output: boundary_out/exp1_matrix_boundary.csv
# ============================================================

cat("=== run_boundary_exp1.R started ===\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("variable_selection.R")

library(mvtnorm)

stopifnot(ls_tests())
cat("ls_tests() passed.\n\n")

dir.create("boundary_out", showWarnings = FALSE)

rmse <- function(a, b) sqrt(mean((a - b)^2))

# ---- 5-covariate truth (for Experiment 1) ----
eta_5cov <- function(X, mu = 0) {
  X <- as.matrix(X); stopifnot(ncol(X) >= 5)
  mu +
    2.0 * sin(pi * X[,1]) +
    1.5 * exp(X[,2] - 0.5) +
    0.7 * (X[,3]^2) +
    0.5 * sin(2 * pi * X[,4]) +
    0.8 * cos(pi * X[,5])
}

# ---- Data generator for Exp1 ----
simulate_exp1 <- function(n, p_real = 5, p_garbage,
                          rho_X = 0.10, nu_X = 1.0,
                          sigma2 = 0.8, rho = 0.2, nu = 1.5,
                          tau2 = 0.15, seed = 42,
                          jitter_X = 1e-8, jitter_b = 1e-8) {
  set.seed(seed)
  p <- p_real + p_garbage

  coords <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1))
  D <- pairdist(coords)

  # Real covariates: spatial GP -> rank -> Unif(0,1)
  R_X <- matern_cor(D, rho = rho_X, nu = nu_X)
  Sigma_X <- R_X + diag(jitter_X, n)

  X_real <- sapply(1:p_real, function(j) {
    z <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_X))
    (rank(z, ties.method = "average") - 0.5) / n
  })
  X_real <- as.matrix(X_real)

  # Garbage covariates: iid Unif(0,1)
  X_garbage <- matrix(runif(n * p_garbage), n, p_garbage)

  X <- cbind(X_real, X_garbage)
  colnames(X) <- paste0("X", 1:p)

  eta <- eta_5cov(X[, 1:p_real, drop = FALSE])

  # Spatial GP residual
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))

  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps

  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# ---- Fit + Wald + summarize ----
run_one_exp1 <- function(n, p_real = 5, p_garbage, M = 6,
                         seed = 42, alpha = 0.05, verbose_fit = FALSE) {
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
    simulate_exp1(n=n, p_real=p_real, p_garbage=p_garbage, seed=seed),
    error = function(e) { cat("  SIM error:", e$message, "\n"); NULL })
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
    error = function(e) { cat("  FIT error:", e$message, "\n"); NULL })

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
    error = function(e) { cat("  WALD error:", e$message, "\n"); NULL })

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
#
# M=6, each covariate = 5 df, k = 1 + p_total*5
# p_real = 5 (all strong signals)
#
# n=100  => max p_total = 19  => max p_garbage = 14
# n=200  => max p_total = 39  => max p_garbage = 34
# n=400  => max p_total = 79  => max p_garbage = 74

cat("============================================================\n")
cat("EXPERIMENT 1: Matrix Boundary (p_real=5, M=6)\n")
cat("============================================================\n\n")

M <- 6
exp1 <- data.frame()

for (n in c(100, 200, 400)) {
  p_real <- 5
  p_max <- floor((n - 1) / 5)  # max total covariates
  max_garbage <- p_max - p_real

  # Test at several fractions of capacity + one past boundary
  frac_vec <- c(0.3, 0.5, 0.7, 0.85, 0.95, 1.0)
  pg_vec <- unique(sort(c(
    pmax(0, round(frac_vec * p_max) - p_real),
    max_garbage + 1   # past boundary
  )))

  cat(sprintf("n = %d (max p_total=%d, max p_garbage=%d)\n", n, p_max, max_garbage))

  for (pg in pg_vec) {
    ptot <- p_real + pg
    k <- 1 + ptot * 5
    cat(sprintf("  p_garbage=%2d (p_total=%d, k=%d, df=%d) ... ",
                pg, ptot, k, n - k))
    t0 <- proc.time()
    row <- run_one_exp1(n=n, p_real=p_real, p_garbage=pg, seed=42)
    dt <- (proc.time() - t0)[3]

    cat(sprintf("[%s]", row$status))
    if (!is.na(row$power)) {
      cat(sprintf(" power=%.2f FPR=%.3f", row$power, row$FPR))
      if (nchar(row$missed) > 0) cat(sprintf("  missed: %s", row$missed))
      cat(sprintf("  sigma2=%.3f tau2=%.3f",
                  row$sigma2_hat, row$tau2_hat))
    }
    cat(sprintf(" (%.1fs)\n", dt))

    exp1 <- rbind(exp1, row)
  }
  cat("\n")
}

write.csv(exp1, "boundary_out/exp1_matrix_boundary.csv", row.names=FALSE)
cat("\nSaved boundary_out/exp1_matrix_boundary.csv\n\n")

cat("--- Full table ---\n")
print(exp1[, c("n","p_real","p_garbage","p_total","k","df_resid",
               "status","power","FPR","n_missed","sigma2_hat","tau2_hat")],
      row.names=FALSE, digits=3)

cat("\n\n=== run_boundary_exp1.R finished ===\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
