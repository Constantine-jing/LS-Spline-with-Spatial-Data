# ============================================================
# variable_selection.R
# Block Wald test for covariate selection in the LS-spline
# spatial additive model fitted by fit_ls_spatial()
#
# WHY THIS APPROACH:
# ------------------
# Our model is  y = X_fix * beta + b + eps
# where X_fix = [1 | W_1 | ... | W_p], and we estimate beta
# by GLS:  beta_hat = (X' Sigma^{-1} X)^{-1} X' Sigma^{-1} y
#
# The GLS covariance of beta_hat is:
#   Cov(beta_hat) = sigma2 * (X' Sigma0^{-1} X)^{-1}
#
# where Sigma0 = R + lambda*I  (the "normalized" covariance).
#
# For covariate j, the identified LS-spline block has (M-1)
# coefficients beta_j. To test H0: beta_j = 0, we use the
# block Wald statistic:
#
#   W_j = beta_j' * [Cov(beta_j)]^{-1} * beta_j
#
# Under H0 (and Gaussian errors), W_j / (M-1) ~ F(M-1, n-k)
# approximately, where k = ncol(X_fix).
#
# This is:
#   - Exact (no simulation/permutation needed)
#   - Free (computed from quantities already in the fit)
#   - Handles the spatial correlation correctly (via GLS)
#   - Gives a p-value per covariate block
#   - Natural threshold: reject H0 at level alpha
#
# As garbage count grows, real covariates should have tiny
# p-values while garbage covariates should be uniform on [0,1].
# We can also apply Bonferroni or BH correction for multiple
# testing across p covariates.
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")

# ------------------------------------------------------------
# block_wald_test(obj)
#
# INPUT:  obj = output from fit_ls_spatial()
# OUTPUT: data.frame with one row per covariate:
#         var, block_size, Wald_stat, F_stat, df1, df2, p_value
# ------------------------------------------------------------
block_wald_test <- function(obj, var_names = NULL) {

  beta_all <- as.numeric(obj$fit$beta)
  col_map  <- obj$des$col_map
  p <- length(col_map)
  n <- nrow(obj$X_fix)
  k <- ncol(obj$X_fix)   # total fixed-effect columns (intercept + all spline blocks)

  if (is.null(var_names)) var_names <- paste0("X", 1:p)

  # --- GLS covariance of beta_hat ---
  # Sigma0 = R + lambda * I  (already computed in fit)
  # Cov(beta_hat) = sigma2 * (X_fix' Sigma0^{-1} X_fix)^{-1}
  #
  # We need Sigma0^{-1} X_fix.  Use Cholesky of Sigma0.
  Sigma0 <- obj$fit$Sigma0
  jitter <- 1e-8
  U0 <- chol(Sigma0 + diag(jitter, n))

  # Whiten: X_t = L^{-1} X_fix
  X_t <- forwardsolve(t(U0), obj$X_fix)
  XtX <- crossprod(X_t)                     # = X_fix' Sigma0^{-1} X_fix
  V_beta <- obj$fit$sigma2 * solve(XtX)     # full covariance of beta_hat

  # --- Block Wald test for each covariate j ---
  results <- data.frame(
    var        = var_names,
    block_size = integer(p),
    Wald_stat  = numeric(p),
    F_stat     = numeric(p),
    df1        = integer(p),
    df2        = integer(p),
    p_value    = numeric(p),
    stringsAsFactors = FALSE
  )

  for (j in 1:p) {
    # Column indices in beta_all: +1 for intercept
    idx_in_beta <- 1 + col_map[[j]]
    beta_j <- beta_all[idx_in_beta]
    q_j <- length(beta_j)    # = M_j - 1

    # Extract the q_j x q_j sub-block of V_beta
    V_j <- V_beta[idx_in_beta, idx_in_beta, drop = FALSE]

    # Wald statistic: beta_j' V_j^{-1} beta_j
    W_j <- as.numeric(t(beta_j) %*% solve(V_j, beta_j))

    # F statistic: W_j / q_j  ~  F(q_j, n - k) under H0
    F_j <- W_j / q_j
    df1 <- q_j
    df2 <- n - k
    pval <- pf(F_j, df1, df2, lower.tail = FALSE)

    results$block_size[j] <- q_j
    results$Wald_stat[j]  <- W_j
    results$F_stat[j]     <- F_j
    results$df1[j]        <- df1
    results$df2[j]        <- df2
    results$p_value[j]    <- pval
  }

  results
}

# ------------------------------------------------------------
# select_covariates(obj, alpha = 0.05, method = "BH")
#
# Applies block Wald test + multiple testing correction.
# method: "none", "bonferroni", "BH" (Benjamini-Hochberg)
#
# Returns the Wald table augmented with:
#   p_adjusted, selected (TRUE/FALSE)
# ------------------------------------------------------------
select_covariates <- function(obj, var_names = NULL,
                               alpha = 0.05,
                               method = c("BH", "bonferroni", "none")) {
  method <- match.arg(method)

  tab <- block_wald_test(obj, var_names = var_names)

  if (method == "none") {
    tab$p_adjusted <- tab$p_value
  } else if (method == "bonferroni") {
    tab$p_adjusted <- pmin(tab$p_value * nrow(tab), 1)
  } else if (method == "BH") {
    tab$p_adjusted <- p.adjust(tab$p_value, method = "BH")
  }

  tab$selected <- tab$p_adjusted < alpha
  tab
}


# ============================================================
# TESTING: run on Sim5 with known real/garbage split
# ============================================================

if (FALSE) {
  # Source sim5
  source("sim5_sweep.R")

  n <- 1000; M <- 6; seed <- 42
  p_real <- 10; p_garbage <- 10; p <- p_real + p_garbage

  dat <- simulate_sim5(n = n, p_real = p_real, p_garbage = p_garbage,
                       sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15, seed = seed)

  coords <- as.matrix(dat[, c("x", "y_coord")])
  X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
  y      <- dat$Y

  obj <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                        M_vec = rep(M, p), nu = 1.5,
                        rho_init = 0.2, lambda_init = 0.15/0.8, verbose = TRUE)

  var_names <- paste0("X", 1:p)
  sel <- select_covariates(obj, var_names = var_names, alpha = 0.05, method = "BH")
  print(sel)
}
