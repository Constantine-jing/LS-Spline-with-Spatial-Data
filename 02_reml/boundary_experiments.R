# ============================================================
# boundary_experiments.R
#
# is the local/interactive version
# 
# Three boundary experiments for variable selection in the
# LS-spline + spatial GP additive model.
#
# SETUP:
#   - Block Wald F-test per covariate (from variable_selection.R)
#   - BH-corrected p-values at alpha = 0.05
#   - M = 6 knots => each covariate uses (M-1) = 5 df
#   - Total fixed-effect columns: k = 1 + p_total * 5
#   - Hard requirement: n > k  (otherwise GLS is impossible)
#
# EXPERIMENT 1 (Matrix boundary):
#   Fix p_real = 5, vary p_garbage so that k approaches n.
#   Small n values: n in {100, 200, 400}.
#   Goal: find where the fit breaks or becomes unreliable
#   as p_total * 5 + 1 -> n.
#
# EXPERIMENT 2 (Statistical boundary):
#   Fix n = 1000 (or 2000), p_real = 10 (the Sim5 truth).
#   Increase p_garbage: 10, 30, 50, 80, 100, 150.
#   Goal: track power and FPR as garbage grows (within the
#   matrix-feasible region).
#
# EXPERIMENT 3 (Signal strength boundary):
#   Fix n = 1000, compare two covariate sets:
#     "strong 10": amplitudes roughly 0.5 to 2.0 (our Sim5 set)
#     "weak 10":   amplitudes scaled down by ~3x (0.15 to 0.6)
#   Sweep p_garbage for each set.
#   Goal: show that the garbage boundary is earlier for weak signals.
#
# THEORETICAL CONTEXT:
# --------------------
# Matrix boundary:
#   The GLS estimator requires X'Sigma0^{-1}X to be invertible,
#   which needs n > k = 1 + p*(M-1). When n - k is small, the
#   residual df (n-k) for sigma2_hat = RSS/(n-k) becomes tiny,
#   making variance estimates unstable. Classical result:
#   Hurvich & Tsai (1989) on AIC correction for small n-k/n.
#
# Statistical boundary:
#   Under H0, the Wald F-stat ~ F(M-1, n-k). As p_garbage grows
#   (k grows), df2 = n-k shrinks, the F distribution becomes
#   heavier-tailed, and more garbage covariates produce large
#   F-stats by chance. With BH correction, the threshold adapts
#   but eventually the weakest real signals become undetectable.
#   Related: Portnoy (1984, 1985) on asymptotics of regression
#   when p grows with n; Fan & Li (2001) on variable selection
#   in high-dimensional regression.
#
# Signal strength boundary:
#   Detection depends on the noncentrality parameter of the
#   Wald statistic: delta_j ~ beta_j' V_j^{-1} beta_j, which
#   scales with ||beta_j||^2. Weaker signals have smaller
#   noncentrality and need either larger n or fewer garbage
#   covariates to be detectable. This connects to the concept
#   of "minimum detectable effect size" in power analysis.
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("variable_selection.R")

library(mvtnorm)

# ============================================================
# Shared simulation generator (flexible signal set)
# ============================================================

# "Strong" truth (same as Sim5)
eta_strong <- function(X, mu = 0) {
  X <- as.matrix(X)
  p <- ncol(X)
  stopifnot(p >= 10)
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

# "Weak" truth: same shapes, amplitudes scaled down by ~3x
eta_weak <- function(X, mu = 0) {
  X <- as.matrix(X)
  p <- ncol(X)
  stopifnot(p >= 10)
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

# Generic data generator with pluggable eta function
simulate_generic <- function(n, p_real, p_garbage,
                             eta_fun = eta_strong,
                             rho_X = 0.10, nu_X = 1.0,
                             sigma2 = 0.8, rho = 0.2, nu = 1.5,
                             tau2 = 0.15, seed = 42,
                             jitter_X = 1e-8, jitter_b = 1e-8) {
  set.seed(seed)
  p <- p_real + p_garbage

  coords <- cbind(
    x = runif(n, 0, 1),
    y = runif(n, 0, 1)
  )

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

  eta <- eta_fun(X[, 1:p_real, drop = FALSE])

  # Spatial GP residual
  R_b <- matern_cor(D, rho = rho, nu = nu)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = sigma2 * R_b + diag(jitter_b, n)))

  eps <- rnorm(n, 0, sqrt(tau2))
  y <- eta + b + eps

  data.frame(x = coords[,1], y_coord = coords[,2],
             Y = y, eta = eta, b = b, eps = eps, X)
}

# ============================================================
# Helper: fit + Wald test + summarize
# ============================================================
run_one_boundary <- function(n, p_real, p_garbage, M = 6,
                              eta_fun = eta_strong,
                              seed = 42, alpha = 0.05,
                              verbose_fit = FALSE) {

  p <- p_real + p_garbage
  k <- 1 + p * (M - 1)

  # Check matrix boundary
  if (k >= n) {
    return(data.frame(
      n = n, p_real = p_real, p_garbage = p_garbage,
      p_total = p, k = k, df_resid = n - k,
      status = "INFEASIBLE",
      power = NA, FPR = NA, n_missed = NA, missed = "",
      rho_hat = NA, sigma2_hat = NA, tau2_hat = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Warn if df is dangerously low
  df_resid <- n - k
  status <- if (df_resid < 20) "LOW_DF" else "OK"

  dat <- tryCatch(
    simulate_generic(n = n, p_real = p_real, p_garbage = p_garbage,
                     eta_fun = eta_fun, seed = seed),
    error = function(e) NULL
  )
  if (is.null(dat)) {
    return(data.frame(
      n = n, p_real = p_real, p_garbage = p_garbage,
      p_total = p, k = k, df_resid = df_resid,
      status = "SIM_FAIL",
      power = NA, FPR = NA, n_missed = NA, missed = "",
      rho_hat = NA, sigma2_hat = NA, tau2_hat = NA,
      stringsAsFactors = FALSE
    ))
  }

  coords <- as.matrix(dat[, c("x", "y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y

  obj <- tryCatch(
    fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                   M_vec = rep(M, p), nu = 1.5,
                   rho_init = 0.2, lambda_init = 0.15/0.8,
                   verbose = verbose_fit),
    error = function(e) NULL
  )

  if (is.null(obj)) {
    return(data.frame(
      n = n, p_real = p_real, p_garbage = p_garbage,
      p_total = p, k = k, df_resid = df_resid,
      status = "FIT_FAIL",
      power = NA, FPR = NA, n_missed = NA, missed = "",
      rho_hat = NA, sigma2_hat = NA, tau2_hat = NA,
      stringsAsFactors = FALSE
    ))
  }

  sel <- tryCatch(
    select_covariates(obj, var_names = paste0("X", 1:p),
                       alpha = alpha, method = "BH"),
    error = function(e) NULL
  )

  if (is.null(sel)) {
    return(data.frame(
      n = n, p_real = p_real, p_garbage = p_garbage,
      p_total = p, k = k, df_resid = df_resid,
      status = "WALD_FAIL",
      power = NA, FPR = NA, n_missed = NA, missed = "",
      rho_hat = obj$fit$rho, sigma2_hat = obj$fit$sigma2,
      tau2_hat = obj$fit$tau2,
      stringsAsFactors = FALSE
    ))
  }

  sel$truth <- c(rep("real", p_real), rep("garbage", p_garbage))
  power <- mean(sel$selected[sel$truth == "real"])
  fpr   <- if (p_garbage > 0) mean(sel$selected[sel$truth == "garbage"]) else 0
  missed <- sel$var[sel$truth == "real" & !sel$selected]

  data.frame(
    n = n, p_real = p_real, p_garbage = p_garbage,
    p_total = p, k = k, df_resid = df_resid,
    status = status,
    power = power, FPR = fpr,
    n_missed = length(missed),
    missed = paste(missed, collapse = ","),
    rho_hat = obj$fit$rho,
    sigma2_hat = obj$fit$sigma2,
    tau2_hat = obj$fit$tau2,
    stringsAsFactors = FALSE
  )
}


# ############################################################
# EXPERIMENT 1: Matrix Boundary
# ############################################################
#
# With M=6, each covariate costs 5 df.
# k = 1 + p_total * 5
# Need n > k, i.e. p_total < (n-1)/5
#
# n=100  => p_total < 19.8  => max ~19
# n=200  => p_total < 39.8  => max ~39
# n=400  => p_total < 79.8  => max ~79
#
# We fix p_real=5 (to keep it small) and push p_garbage
# toward the boundary.

cat("\n")
cat("############################################################\n")
cat("# EXPERIMENT 1: Matrix Boundary                            #\n")
cat("############################################################\n\n")

cat("With M=6: k = 1 + p_total * 5.  Need n > k.\n")
cat("  n=100 => max p_total = 19\n")
cat("  n=200 => max p_total = 39\n")
cat("  n=400 => max p_total = 79\n\n")

exp1_results <- data.frame()

for (n in c(100, 200, 400)) {
  p_max <- floor((n - 1) / 5)   # theoretical max p_total
  p_real <- 5

  # Test at several fractions of capacity
  # going from comfortable to right at the edge
  frac_vec <- c(0.3, 0.5, 0.7, 0.85, 0.95, 1.0)
  p_garbage_vec <- unique(pmax(0, round(frac_vec * p_max) - p_real))
  # also add one past the boundary
  p_garbage_vec <- c(p_garbage_vec, p_max - p_real + 1)
  p_garbage_vec <- sort(unique(p_garbage_vec))

  for (pg in p_garbage_vec) {
    p_total <- p_real + pg
    k <- 1 + p_total * 5
    cat(sprintf("  n=%d, p_real=%d, p_garbage=%d (p_total=%d, k=%d, df=%d) ... ",
                n, p_real, pg, p_total, k, n - k))

    row <- run_one_boundary(n = n, p_real = p_real, p_garbage = pg,
                             eta_fun = eta_strong, seed = 42)
    cat(sprintf("[%s] power=%.2f FPR=%.2f\n",
                row$status,
                ifelse(is.na(row$power), -1, row$power),
                ifelse(is.na(row$FPR), -1, row$FPR)))

    exp1_results <- rbind(exp1_results, row)
  }
  cat("\n")
}

cat("\n--- Experiment 1: Matrix Boundary Summary ---\n")
print(exp1_results[, c("n","p_real","p_garbage","p_total","k","df_resid",
                        "status","power","FPR","sigma2_hat","tau2_hat")],
      row.names = FALSE, digits = 3)


# ############################################################
# EXPERIMENT 2: Statistical Boundary (garbage sweep)
# ############################################################
#
# Fix p_real = 10 (the Sim5 strong signals), fix n.
# Increase p_garbage.
# n=1000 => max p_total = 199 => max p_garbage = 189
# n=2000 => max p_total = 399 => max p_garbage = 389

cat("\n\n")
cat("############################################################\n")
cat("# EXPERIMENT 2: Statistical Boundary (garbage sweep)       #\n")
cat("############################################################\n\n")

exp2_results <- data.frame()

for (n in c(1000, 2000)) {
  p_real <- 10
  p_max_garbage <- floor((n - 1) / 5) - p_real

  # Sweep: start small, go up to ~80% of max
  pg_candidates <- c(10, 20, 40, 60, 80, 100, 130, 160)
  pg_vec <- pg_candidates[pg_candidates <= p_max_garbage]

  cat(sprintf("n = %d (max p_garbage = %d)\n", n, p_max_garbage))

  for (pg in pg_vec) {
    p_total <- p_real + pg
    k <- 1 + p_total * 5
    cat(sprintf("  p_garbage=%3d (k=%d, df=%d) ... ", pg, k, n - k))

    row <- run_one_boundary(n = n, p_real = 10, p_garbage = pg,
                             eta_fun = eta_strong, seed = 42)
    cat(sprintf("[%s] power=%.2f FPR=%.3f",
                row$status, row$power, row$FPR))
    if (nchar(row$missed) > 0) cat(sprintf("  missed: %s", row$missed))
    cat("\n")

    exp2_results <- rbind(exp2_results, row)
  }
  cat("\n")
}

cat("\n--- Experiment 2: Statistical Boundary Summary ---\n")
print(exp2_results[, c("n","p_garbage","p_total","k","df_resid",
                        "status","power","FPR","n_missed","missed")],
      row.names = FALSE, digits = 3)


# ############################################################
# EXPERIMENT 3: Signal Strength Boundary
# ############################################################
#
# Fix n = 1000, p_real = 10.
# Compare "strong" vs "weak" signal sets.
# Sweep p_garbage for each.
# Show that the garbage boundary shifts earlier for weak signals.

cat("\n\n")
cat("############################################################\n")
cat("# EXPERIMENT 3: Signal Strength Boundary                   #\n")
cat("############################################################\n\n")

# First, show the amplitude comparison
cat("Signal amplitudes (approx range on [0,1]):\n")
cat("  Strong set:  X1=2.0, X2=1.5, X3=0.7, X4=0.5, X5=1.0, X6=0.8, X7=1.2, X8=0.6, X9=0.9, X10=0.4\n")
cat("  Weak set:    X1=0.6, X2=0.5, X3=0.25, X4=0.15, X5=0.3, X6=0.25, X7=0.4, X8=0.2, X9=0.3, X10=0.12\n\n")

exp3_results <- data.frame()

n <- 1000
p_real <- 10
pg_vec <- c(10, 20, 40, 60, 80, 100, 130)

for (strength in c("strong", "weak")) {
  eta_fun <- if (strength == "strong") eta_strong else eta_weak

  cat(sprintf("--- %s signals (n=%d) ---\n", toupper(strength), n))

  for (pg in pg_vec) {
    p_total <- p_real + pg
    k <- 1 + p_total * 5
    if (k >= n) { cat(sprintf("  p_garbage=%3d: SKIP (k=%d >= n)\n", pg, k)); next }

    cat(sprintf("  p_garbage=%3d (k=%d, df=%d) ... ", pg, k, n - k))

    row <- run_one_boundary(n = n, p_real = 10, p_garbage = pg,
                             eta_fun = eta_fun, seed = 42)
    row$strength <- strength
    cat(sprintf("[%s] power=%.2f FPR=%.3f",
                row$status, row$power, row$FPR))
    if (nchar(row$missed) > 0) cat(sprintf("  missed: %s", row$missed))
    cat("\n")

    exp3_results <- rbind(exp3_results, row)
  }
  cat("\n")
}

cat("\n--- Experiment 3: Signal Strength Comparison ---\n")
print(exp3_results[, c("strength","n","p_garbage","p_total","df_resid",
                        "power","FPR","n_missed","missed")],
      row.names = FALSE, digits = 3)

# Side-by-side comparison
cat("\n--- Strong vs Weak at each p_garbage ---\n")
for (pg in pg_vec) {
  s_row <- exp3_results[exp3_results$strength == "strong" & exp3_results$p_garbage == pg, ]
  w_row <- exp3_results[exp3_results$strength == "weak"   & exp3_results$p_garbage == pg, ]
  if (nrow(s_row) == 0 || nrow(w_row) == 0) next
  cat(sprintf("  p_garbage=%3d:  strong power=%.2f FPR=%.3f  |  weak power=%.2f FPR=%.3f\n",
              pg, s_row$power, s_row$FPR, w_row$power, w_row$FPR))
}


# ############################################################
# Summary
# ############################################################
cat("\n\n")
cat("############################################################\n")
cat("# SUMMARY                                                  #\n")
cat("############################################################\n\n")
cat("1. MATRIX BOUNDARY: k = 1 + p_total*(M-1) must be < n.\n")
cat("   With M=6: p_total < (n-1)/5.\n")
cat("   As df_resid = n-k approaches 0, sigma2_hat and Wald\n")
cat("   tests become unreliable even before exact singularity.\n\n")
cat("2. STATISTICAL BOUNDARY: even within the feasible region,\n")
cat("   power to detect weak real signals drops as p_garbage\n")
cat("   grows, because (a) df_resid shrinks and (b) BH\n")
cat("   correction becomes more conservative with more tests.\n\n")
cat("3. SIGNAL STRENGTH: the garbage boundary (where power\n")
cat("   drops below, say, 0.8) occurs at a smaller p_garbage\n")
cat("   for weak signals than for strong signals.\n\n")

cat("Done.\n")
