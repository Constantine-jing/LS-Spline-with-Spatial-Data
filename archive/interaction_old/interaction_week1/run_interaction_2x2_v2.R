# ============================================================
# run_interaction_2x2_v2.R
#
# UNIFIED 2x2 EXPERIMENT (centered Sim4 DGP, p in {2, 3, 4})
#
# DGP (centered Sim4-style):
#   f1(x) = 2.0 * sin(pi * x)
#   f2(x) = 1.5 * (exp(x - 0.5) - C2),    C2 = exp(0.5) - exp(-0.5)
#   f3(x) = 0.7 * (x^2 - 1/3)
#   f4(x) = 0.5 * sin(2 * pi * x)
#   f12(x1,x2) = c_int * sin(pi*x1) * (exp(x2-0.5) - C2)
#   All other pairwise interactions = 0.
#   b ~ Matern(rho=0.3, nu=1.5), sigma2=0.5; eps ~ N(0, tau2=0.25)
#   X_j ~ Uniform(0,1) i.i.d.; locations ~ Uniform([0,1]^2)
#
# Two-phase rollout for each p in {2,3,4}:
#   Phase 1: M1 includes only X1xX2 interaction
#   Phase 2: M1 includes all p*(p-1)/2 pairwise interactions
#
# 2x2 design per (p, phase):
#   {DGP-NULL (c_int=0), DGP-INT (c_int=1.5)} x {Model M0, Model M1}
#
# Output:
#   interaction_2x2_v2_summary.csv          (one row per cell, 24 rows)
#   interaction_2x2_v2_results.rds          (full posterior summaries)
#   interaction_2x2_v2_p{p}_phase{ph}_marginals.pdf
#   interaction_2x2_v2_p{p}_phase{ph}_surfaces.pdf
# ============================================================


# ============================================================
# CONSTANTS for centered DGP
# ============================================================
.C2_F2 <- exp(0.5) - exp(-0.5)   # makes f2 centered around 0 over [0,1]


# ============================================================
# True functions (centered)
# ============================================================
true_f_unified <- function(j, x) {
  switch(j,
    `1` = 2.0 * (sin(pi * x) - 2/pi),
    `2` = 1.5 * (exp(x - 0.5) - .C2_F2),
    `3` = 0.7 * (x^2 - 1/3),
    `4` = 0.5 * sin(2 * pi * x),
    stop("j must be in 1..4")
  )
}

true_f12_unified <- function(x1, x2, c_int) {
  c_int * (sin(pi * x1) - 2/pi) * (exp(x2 - 0.5) - .C2_F2)
}


# ============================================================
# WAIC computation (relative metric, for M0-vs-M1 comparison)
# ============================================================
compute_waic <- function(logLik_mat) {
  lpd  <- sum(log(colMeans(exp(logLik_mat))))
  pwaic <- sum(apply(logLik_mat, 2, var))
  list(waic = -2 * (lpd - pwaic), lpd = lpd, p_waic = pwaic)
}

pointwise_loglik <- function(y, fitted_mat, sigma2_vec, tau2_vec) {
  n_keep <- nrow(fitted_mat); n <- length(y)
  ll_mat <- matrix(NA, n_keep, n)
  for (k in seq_len(n_keep)) {
    sd_k <- sqrt(sigma2_vec[k] + tau2_vec[k])
    ll_mat[k, ] <- dnorm(y, mean = fitted_mat[k, ], sd = sd_k, log = TRUE)
  }
  ll_mat
}


# ============================================================
# simulate_unified(): generate data for given p
# ============================================================
simulate_unified <- function(n, p, c_int, data_seed,
                              sigma2 = 0.5, rho = 0.3, tau2 = 0.25, nu = 1.5) {
  stopifnot(p %in% c(2, 3, 4))
  set.seed(data_seed)
  locs <- matrix(runif(2 * n), n, 2)
  D <- as.matrix(dist(locs))

  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("X", seq_len(p))

  # Main effects
  f_true <- vector("list", p)
  for (j in seq_len(p)) f_true[[j]] <- true_f_unified(j, X[, j])
  names(f_true) <- paste0("true_f", seq_len(p))

  # Only X1xX2 interaction is true-nonzero
  f12_true <- true_f12_unified(X[, 1], X[, 2], c_int)

  # Spatial GP
  R_true <- matern_cor(D, rho = rho, nu = nu)
  L_R <- chol(R_true + diag(1e-8, n))
  b <- as.vector(t(L_R) %*% rnorm(n)) * sqrt(sigma2)
  eps <- rnorm(n, 0, sqrt(tau2))

  y <- Reduce("+", f_true) + f12_true + b + eps

  out <- list(y = y, D = D, locs = locs, X = X, b_true = b,
              true_f12 = f12_true,
              true_sigma2 = sigma2, true_tau2 = tau2,
              true_rho = rho, true_nu = nu, c_int = c_int, p = p)
  for (j in seq_len(p)) out[[paste0("X", j)]] <- X[, j]
  for (j in seq_len(p)) out[[paste0("true_f", j)]] <- f_true[[j]]
  out
}


# ============================================================
# build_design_unified(): unified design builder
#
# include_int:
#   "none"      -> M0 (no interactions)
#   "x1x2_only" -> M1 with only X1xX2 interaction
#   "all_pairs" -> M1 with all p*(p-1)/2 pairwise interactions
# ============================================================
build_design_unified <- function(X, M, include_int = c("none", "x1x2_only", "all_pairs")) {
  include_int <- match.arg(include_int)
  p <- ncol(X)

  # Build per-covariate LS objects
  objs <- vector("list", p)
  K_main_list <- vector("list", p)
  W_list <- vector("list", p)
  d_main <- integer(p)

  K_raw <- build_rw2_penalty_1d(M)
  for (j in seq_len(p)) {
    objs[[j]] <- ls_build_one_full(X[, j], M = M)
    W_list[[j]] <- objs[[j]]$W
    K_main_list[[j]] <- t(objs[[j]]$T) %*% K_raw %*% objs[[j]]$T
    d_main[j] <- ncol(W_list[[j]])
  }

  # Determine interaction pairs
  if (include_int == "none") {
    pairs <- list()
  } else if (include_int == "x1x2_only") {
    pairs <- list(c(1, 2))
  } else {
    pairs <- list()
    for (a in 1:(p - 1)) for (b in (a + 1):p) pairs[[length(pairs) + 1]] <- c(a, b)
  }

  # Build interaction blocks
  int_recipes <- list()
  K_int_list <- list()
  W_int_list <- list()
  d_int <- integer(0)
  for (k in seq_along(pairs)) {
    a <- pairs[[k]][1]; b <- pairs[[k]][2]
    int_ab <- ls_build_interaction(objs[[a]], objs[[b]])
    int_recipes[[k]] <- int_ab
    K_int_list[[k]] <- int_ab$K_uv
    W_int_list[[k]] <- int_ab$W_uv
    d_int[k] <- ncol(int_ab$W_uv)
  }

  # Assemble H = [1 | W1 ... Wp | W_int_1 ... W_int_K]
  H <- do.call(cbind, c(list(1), W_list, W_int_list))

  # 1-indexed col_maps
  col_map_main <- vector("list", p)
  s <- 1L  # last filled column (intercept = col 1)
  for (j in seq_len(p)) {
    col_map_main[[j]] <- seq(s + 1, s + d_main[j])
    s <- s + d_main[j]
  }
  col_map_int <- vector("list", length(pairs))
  for (k in seq_along(pairs)) {
    col_map_int[[k]] <- seq(s + 1, s + d_int[k])
    s <- s + d_int[k]
  }

  # Convert to 0-indexed (sampler convention)
  col_map_main <- lapply(col_map_main, function(x) x - 1L)
  col_map_int  <- lapply(col_map_int,  function(x) x - 1L)

  list(H = H,
       col_map_main = col_map_main, col_map_int = col_map_int,
       K_main_list  = K_main_list,  K_int_list  = K_int_list,
       objs = objs, int_recipes = int_recipes,
       d_main = d_main, d_int = d_int,
       pairs = pairs, include_int = include_int)
}


# ============================================================
# fit_one_cell(): run sampler for one cell with timing
# ============================================================
fit_one_cell <- function(y, design, D, n_iter, n_burn, mcmc_seed,
                          true_nu = 1.5, verbose = FALSE) {
  set.seed(mcmc_seed)
  t0 <- proc.time()
  gs <- gibbs_interaction_sampler(
    y = y, H = design$H, D = D, nu = true_nu,
    col_map_main = design$col_map_main,
    K_main_list  = design$K_main_list,
    col_map_int  = design$col_map_int,
    K_int_list   = design$K_int_list,
    n_iter = n_iter, n_burn = n_burn, verbose = verbose
  )
  list(gs = gs, elapsed = as.numeric((proc.time() - t0)["elapsed"]))
}


# ============================================================
# summarise_cell_unified(): RMSE + WAIC + grid posteriors
# ============================================================
summarise_cell_unified <- function(fit, design, sim,
                                    grid_main = NULL, grid_int = NULL) {
  gs <- fit$gs
  H  <- design$H
  y  <- sim$y
  p  <- sim$p
  d_main <- design$d_main
  d_int  <- design$d_int
  has_int <- length(design$col_map_int) > 0
  pairs   <- design$pairs

  # Eta indices in H
  starts <- cumsum(c(2, d_main))[1:p]
  ends   <- starts + d_main - 1
  idx_main_list <- mapply(seq, starts, ends, SIMPLIFY = FALSE)

  if (has_int) {
    s <- max(unlist(idx_main_list))
    idx_int_list <- vector("list", length(pairs))
    for (k in seq_along(pairs)) {
      idx_int_list[[k]] <- seq(s + 1, s + d_int[k])
      s <- s + d_int[k]
    }
  } else {
    idx_int_list <- list()
  }

  eta_mean <- colMeans(gs$eta_samples)

  # Main-effect fitted at training X
  rmse_main <- numeric(p)
  for (j in seq_len(p)) {
    Wj <- design$objs[[j]]$W
    f_hat <- as.vector(Wj %*% eta_mean[idx_main_list[[j]]])
    rmse_main[j] <- sqrt(mean((f_hat - sim[[paste0("true_f", j)]])^2))
  }
  names(rmse_main) <- paste0("RMSE_f", seq_len(p))

  # Marginal grid posteriors
  marginal <- NULL
  if (!is.null(grid_main)) {
    marginal <- vector("list", p)
    for (j in seq_len(p)) {
      x_grid <- grid_main[[j]]
      Wj_grid <- design$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)
      beta_j_post <- gs$eta_samples[, idx_main_list[[j]], drop = FALSE]
      f_grid_samp <- beta_j_post %*% t(Wj_grid)
      marginal[[j]] <- list(
        x         = x_grid,
        post_mean = colMeans(f_grid_samp),
        post_lo   = apply(f_grid_samp, 2, quantile, 0.025),
        post_hi   = apply(f_grid_samp, 2, quantile, 0.975)
      )
    }
    names(marginal) <- paste0("f", seq_len(p))
  }

  # Interaction (X1xX2 only — first interaction block by convention)
  rmse_f12 <- NA; cov_f12 <- NA
  surface12 <- NULL
  f12_train_hat <- NULL; f12_train_lo <- NULL; f12_train_hi <- NULL
  if (has_int) {
    # Find which block is X1xX2 (in any phase, it's pairs[[1]] = c(1,2))
    k12 <- which(sapply(pairs, function(p_pair) all(p_pair == c(1, 2))))
    if (length(k12) == 1) {
      beta12_post <- gs$eta_samples[, idx_int_list[[k12]], drop = FALSE]
      W12 <- design$int_recipes[[k12]]$W_uv
      f12_samp <- beta12_post %*% t(W12)
      f12_train_hat <- colMeans(f12_samp)
      f12_train_lo  <- apply(f12_samp, 2, quantile, 0.025)
      f12_train_hi  <- apply(f12_samp, 2, quantile, 0.975)
      rmse_f12 <- sqrt(mean((f12_train_hat - sim$true_f12)^2))
      cov_f12  <- mean(f12_train_lo <= sim$true_f12 & sim$true_f12 <= f12_train_hi)

      # Surface on regular grid
      if (!is.null(grid_int)) {
        W12_grid <- ls_interaction_design_new(grid_int$X1_long, grid_int$X2_long,
                                              design$int_recipes[[k12]]$recipe,
                                              clip = TRUE)
        surf_samp <- beta12_post %*% t(W12_grid)
        surface12 <- list(
          x1 = grid_int$x1, x2 = grid_int$x2,
          post_mean = matrix(colMeans(surf_samp),
                             length(grid_int$x1), length(grid_int$x2)),
          post_lo   = matrix(apply(surf_samp, 2, quantile, 0.025),
                             length(grid_int$x1), length(grid_int$x2)),
          post_hi   = matrix(apply(surf_samp, 2, quantile, 0.975),
                             length(grid_int$x1), length(grid_int$x2))
        )
        surface12$band_width <- surface12$post_hi - surface12$post_lo
      }
    }
  }

  fitted_mat <- gs$eta_samples %*% t(H)
  ll_mat <- pointwise_loglik(y, fitted_mat, gs$sigma2_samples, gs$tau2_samples)
  waic_obj <- compute_waic(ll_mat)
  rmse_y <- sqrt(mean((colMeans(fitted_mat) - y)^2))

  list(
    rmse_main   = rmse_main,
    rmse_f12    = rmse_f12, cov_f12 = cov_f12,
    rmse_y      = rmse_y,
    waic        = waic_obj$waic, p_waic = waic_obj$p_waic,
    sigma2_mean = mean(gs$sigma2_samples),
    tau2_mean   = mean(gs$tau2_samples),
    rho_mean    = mean(gs$rho_samples),
    accept      = if (!is.null(gs$accept_rate)) gs$accept_rate
                  else c(sigma2 = NA, tau2 = NA, rho = NA),
    elapsed     = fit$elapsed,
    n_keep      = nrow(gs$eta_samples),
    marginal    = marginal,
    surface12   = surface12,
    f12_train_hat = f12_train_hat,
    f12_train_lo  = f12_train_lo,
    f12_train_hi  = f12_train_hi,
    n_int_blocks  = length(design$col_map_int),
    # Raw MCMC chains for trace/density/ACF plots
    chains = list(
      sigma2         = gs$sigma2_samples,
      tau2           = gs$tau2_samples,
      rho            = gs$rho_samples,
      tau2_s_main    = gs$tau2_s_main_samp,
      tau2_s_int     = gs$tau2_s_int_samp
    )
  )
}


# ============================================================
# print_2x2_table_v2(): pretty-print the 2x2 summary
# ============================================================
print_2x2_table_v2 <- function(results, p, phase) {
  cells <- list(A = results$null_M0, B = results$null_M1,
                C = results$int_M0,  D = results$int_M1)
  cat("\n========================================================\n")
  cat(sprintf("  RESULTS: p = %d  |  Phase %d  (%s)\n", p, phase,
              ifelse(phase == 1, "M1 = X1xX2 only", "M1 = all pairwise")))
  cat("========================================================\n")
  cat(sprintf("%-12s | %-22s | %-22s\n", "", "M0 (no interaction)", "M1 (with interaction)"))
  cat(strrep("-", 64), "\n", sep = "")
  for (dgp in c("NULL", "INT")) {
    cM0 <- if (dgp == "NULL") cells$A else cells$C
    cM1 <- if (dgp == "NULL") cells$B else cells$D
    rm0 <- paste(sprintf("%.3f", cM0$rmse_main), collapse = ",")
    rm1 <- paste(sprintf("%.3f", cM1$rmse_main), collapse = ",")
    cat(sprintf("%-12s | RMSE_main=%-12s | RMSE_main=%-12s\n",
                paste0("DGP-", dgp), rm0, rm1))
    cat(sprintf("%-12s | RMSE_y   =%-12.4f | RMSE_y   =%-12.4f\n",
                "", cM0$rmse_y, cM1$rmse_y))
    cat(sprintf("%-12s | sigma2   =%-12.3f | sigma2   =%-12.3f\n",
                "", cM0$sigma2_mean, cM1$sigma2_mean))
    cat(sprintf("%-12s | tau2     =%-12.3f | tau2     =%-12.3f\n",
                "", cM0$tau2_mean, cM1$tau2_mean))
    cat(sprintf("%-12s | rho      =%-12.3f | rho      =%-12.3f\n",
                "", cM0$rho_mean, cM1$rho_mean))
    cat(sprintf("%-12s | WAIC     =%-12.1f | WAIC     =%-12.1f\n",
                "", cM0$waic, cM1$waic))
    if (!is.na(cM1$rmse_f12)) {
      cat(sprintf("%-12s | f12      = (n/a)        | RMSE_f12 =%-12.4f\n",
                  "", cM1$rmse_f12))
      cat(sprintf("%-12s |                         | cov_f12  =%-12.3f\n",
                  "", cM1$cov_f12))
    }
    cat(sprintf("%-12s | time(s)  =%-12.1f | time(s)  =%-12.1f\n",
                "", cM0$elapsed, cM1$elapsed))
    cat(strrep("-", 64), "\n", sep = "")
  }
  cat("\nKey comparisons:\n")
  cat(sprintf("  WAIC delta (NULL: M1-M0) = %+.1f  [>0 means M0 preferred]\n",
              cells$B$waic - cells$A$waic))
  cat(sprintf("  WAIC delta (INT:  M1-M0) = %+.1f  [<0 means M1 preferred]\n",
              cells$D$waic - cells$C$waic))
  cat(sprintf("  Time mult M1/M0 (INT)   = %.2fx\n",
              cells$D$elapsed / cells$C$elapsed))
}


# ============================================================
# run_2x2(): unified runner for any (p, phase)
# ============================================================
run_2x2 <- function(p, phase, n = 300, M = 8,
                    n_iter = 3000, n_burn = 1000, c_int = 1.5,
                    data_seed_null = NULL, data_seed_int = NULL,
                    mcmc_seed = NULL, verbose = FALSE) {
  stopifnot(p %in% c(2, 3, 4), phase %in% c(1, 2))

  # Default seeds depend on p, phase for reproducibility
  if (is.null(data_seed_null)) data_seed_null <- 100L * p + phase
  if (is.null(data_seed_int))  data_seed_int  <- 200L * p + phase
  if (is.null(mcmc_seed))      mcmc_seed      <- 300L * p + phase

  cat(sprintf("\n##########  p = %d  PHASE %d  ##########\n", p, phase))
  cat(sprintf("n=%d  M=%d  n_iter=%d  n_burn=%d  c_int=%.2f\n",
              n, M, n_iter, n_burn, c_int))
  inc_M1 <- if (phase == 1) "x1x2_only" else "all_pairs"
  cat(sprintf("M1 includes interactions: %s\n", inc_M1))

  # Simulate two DGPs
  sim_null <- simulate_unified(n, p, c_int = 0,     data_seed = data_seed_null)
  sim_int  <- simulate_unified(n, p, c_int = c_int, data_seed = data_seed_int)

  # Build designs (M0 same for both DGPs of given p; M1 same too, just X varies)
  des_M0_null <- build_design_unified(sim_null$X, M, include_int = "none")
  des_M1_null <- build_design_unified(sim_null$X, M, include_int = inc_M1)
  des_M0_int  <- build_design_unified(sim_int$X,  M, include_int = "none")
  des_M1_int  <- build_design_unified(sim_int$X,  M, include_int = inc_M1)

  cat(sprintf("\n  H dims: M0 = %d x %d  | M1 = %d x %d\n",
              nrow(des_M0_null$H), ncol(des_M0_null$H),
              nrow(des_M1_null$H), ncol(des_M1_null$H)))

  cat("\nFitting Cell A (NULL + M0)...\n")
  fit_A <- fit_one_cell(sim_null$y, des_M0_null, sim_null$D,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_A$elapsed))

  cat("Fitting Cell B (NULL + M1)...\n")
  fit_B <- fit_one_cell(sim_null$y, des_M1_null, sim_null$D,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_B$elapsed))

  cat("Fitting Cell C (INT + M0)...\n")
  fit_C <- fit_one_cell(sim_int$y, des_M0_int, sim_int$D,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_C$elapsed))

  cat("Fitting Cell D (INT + M1)...\n")
  fit_D <- fit_one_cell(sim_int$y, des_M1_int, sim_int$D,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_D$elapsed))

  # Grids for plotting
  g1d <- seq(0, 1, length.out = 200)
  grid_main <- replicate(p, g1d, simplify = FALSE)
  g2d <- seq(0, 1, length.out = 30)
  gg2d <- expand.grid(x1 = g2d, x2 = g2d)
  grid_int <- list(x1 = g2d, x2 = g2d, X1_long = gg2d$x1, X2_long = gg2d$x2)

  truth_null <- list()
  truth_int  <- list()
  for (j in seq_len(p)) {
    truth_null[[paste0("f", j)]] <- true_f_unified(j, g1d)
    truth_int [[paste0("f", j)]] <- true_f_unified(j, g1d)
  }
  truth_null$f12 <- matrix(0, length(g2d), length(g2d))
  truth_int$f12  <- matrix(true_f12_unified(gg2d$x1, gg2d$x2, c_int),
                           length(g2d), length(g2d))

  results <- list(
    null_M0 = summarise_cell_unified(fit_A, des_M0_null, sim_null,
                                      grid_main = grid_main, grid_int = NULL),
    null_M1 = summarise_cell_unified(fit_B, des_M1_null, sim_null,
                                      grid_main = grid_main, grid_int = grid_int),
    int_M0  = summarise_cell_unified(fit_C, des_M0_int,  sim_int,
                                      grid_main = grid_main, grid_int = NULL),
    int_M1  = summarise_cell_unified(fit_D, des_M1_int,  sim_int,
                                      grid_main = grid_main, grid_int = grid_int)
  )

  print_2x2_table_v2(results, p = p, phase = phase)

  invisible(list(results = results,
                 sim_null = sim_null, sim_int = sim_int,
                 settings = list(n = n, M = M, p = p, phase = phase,
                                 n_iter = n_iter, n_burn = n_burn,
                                 c_int = c_int),
                 truth_grid_null = truth_null,
                 truth_grid_int  = truth_int))
}


# ============================================================
# results_to_dataframe_v2(): tidy CSV-ready frame
# ============================================================
results_to_dataframe_v2 <- function(res_obj) {
  s <- res_obj$settings
  results <- res_obj$results
  cells <- list(
    list(name = "A", dgp = "NULL", model = "M0", obj = results$null_M0),
    list(name = "B", dgp = "NULL", model = "M1", obj = results$null_M1),
    list(name = "C", dgp = "INT",  model = "M0", obj = results$int_M0),
    list(name = "D", dgp = "INT",  model = "M1", obj = results$int_M1)
  )
  rows <- lapply(cells, function(cell) {
    rm <- cell$obj$rmse_main
    data.frame(
      p = s$p, phase = s$phase, n = s$n, M = s$M,
      n_iter = s$n_iter, n_burn = s$n_burn, c_int_true = s$c_int,
      cell = cell$name, DGP = cell$dgp, Model = cell$model,
      n_int_blocks = cell$obj$n_int_blocks,
      RMSE_f1 = rm[1],
      RMSE_f2 = if (length(rm) >= 2) rm[2] else NA,
      RMSE_f3 = if (length(rm) >= 3) rm[3] else NA,
      RMSE_f4 = if (length(rm) >= 4) rm[4] else NA,
      RMSE_f12 = cell$obj$rmse_f12,
      cov_f12  = cell$obj$cov_f12,
      RMSE_y   = cell$obj$rmse_y,
      sigma2_mean = cell$obj$sigma2_mean,
      tau2_mean   = cell$obj$tau2_mean,
      rho_mean    = cell$obj$rho_mean,
      WAIC = cell$obj$waic, p_WAIC = cell$obj$p_waic,
      accept_sigma2 = cell$obj$accept[1],
      accept_tau2   = cell$obj$accept[2],
      accept_rho    = cell$obj$accept[3],
      elapsed_sec   = cell$obj$elapsed,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows); rownames(df) <- NULL
  df
}


# ============================================================
# plot_marginals_v2(): one full page per main effect, 2x2 cell layout
# ============================================================
plot_marginals_v2 <- function(res_obj, outfile) {
  s <- res_obj$settings
  results <- res_obj$results
  truth_null <- res_obj$truth_grid_null
  truth_int  <- res_obj$truth_grid_int

  cells <- list(
    list(name = "NULL + M0 (Cell A)", obj = results$null_M0, truth = truth_null),
    list(name = "NULL + M1 (Cell B)", obj = results$null_M1, truth = truth_null),
    list(name = "INT  + M0 (Cell C)", obj = results$int_M0,  truth = truth_int),
    list(name = "INT  + M1 (Cell D)", obj = results$int_M1,  truth = truth_int)
  )

  pdf(outfile, width = 11, height = 8.5)
  for (j in seq_len(s$p)) {
    fname <- paste0("f", j)
    y_all <- c()
    for (cell in cells) {
      m <- cell$obj$marginal[[fname]]
      if (!is.null(m)) y_all <- c(y_all, m$post_lo, m$post_hi, cell$truth[[fname]])
    }
    ylim <- range(y_all, na.rm = TRUE) * 1.05

    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3.5, 1.5),
        oma = c(0, 0, 3, 0), cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)

    for (cell in cells) {
      m <- cell$obj$marginal[[fname]]
      if (is.null(m)) { plot.new(); next }
      x <- m$x
      true_y <- cell$truth[[fname]]
      plot(x, true_y, type = "n", ylim = ylim,
           xlab = paste0("X", j), ylab = paste0("f", j, "(X", j, ")"),
           main = sprintf("%s   RMSE = %.3f", cell$name, cell$obj$rmse_main[j]))
      polygon(c(x, rev(x)), c(m$post_lo, rev(m$post_hi)),
              col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
      lines(x, m$post_mean, col = "red",   lwd = 2.5, lty = 2)
      lines(x, true_y,      col = "black", lwd = 2.5)
      legend("topright", c("truth", "posterior mean", "95% CI"),
             col = c("black", "red", adjustcolor("steelblue", alpha.f = 0.5)),
             lwd = c(2.5, 2.5, 8), lty = c(1, 2, 1), bty = "n", cex = 0.95)
      grid()
    }
    mtext(sprintf("Marginal recovery: f%d  (p = %d, phase %d)", j, s$p, s$phase),
          outer = TRUE, cex = 1.5, font = 2)
  }
  dev.off()
  cat(sprintf("  Marginal plots: %s\n", outfile))
}


# ============================================================
# plot_surfaces_v2(): one heatmap per page for f12
# ============================================================
plot_surfaces_v2 <- function(res_obj, outfile) {
  s <- res_obj$settings
  results <- res_obj$results
  truth_int <- res_obj$truth_grid_int$f12

  surfD <- results$int_M1$surface12
  surfB <- results$null_M1$surface12

  pdf(outfile, width = 9, height = 8)
  if (is.null(surfD)) {
    plot.new(); text(0.5, 0.5, "No interaction surface stored")
    dev.off(); return(invisible(NULL))
  }
  x1 <- surfD$x1; x2 <- surfD$x2
  z_all <- c(as.vector(truth_int), as.vector(surfD$post_mean),
             if (!is.null(surfB)) as.vector(surfB$post_mean) else numeric(0))
  z_lim <- max(abs(z_all))
  z_breaks <- seq(-z_lim, z_lim, length.out = 51)
  pal <- hcl.colors(50, "Blue-Red 2")

  par(mar = c(5, 5, 4, 6), cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)

  image(x1, x2, truth_int, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("TRUE f_12  (p=%d, phase %d)\nrange [%.2f, %.2f]",
                       s$p, s$phase, min(truth_int), max(truth_int)))
  contour(x1, x2, truth_int, add = TRUE, lwd = 1.2, labcex = 0.9)

  image(x1, x2, surfD$post_mean, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("POSTERIOR MEAN f_12  (Cell D: INT + M1)\nRMSE = %.3f, cov of truth = %.2f",
                       results$int_M1$rmse_f12, results$int_M1$cov_f12))
  contour(x1, x2, surfD$post_mean, add = TRUE, lwd = 1.2, labcex = 0.9)

  pal_unc <- hcl.colors(50, "YlOrRd", rev = TRUE)
  image(x1, x2, surfD$band_width, col = pal_unc,
        xlab = "X1", ylab = "X2",
        main = sprintf("UNCERTAINTY: 95%% band width  (Cell D)\nmean = %.3f, max = %.3f",
                       mean(surfD$band_width), max(surfD$band_width)))
  contour(x1, x2, surfD$band_width, add = TRUE, lwd = 1.2, labcex = 0.9)

  if (!is.null(surfB)) {
    image(x1, x2, surfB$post_mean, col = pal, breaks = z_breaks,
          xlab = "X1", ylab = "X2",
          main = sprintf("POSTERIOR MEAN f_12  (Cell B: NULL + M1, truth=0)\nRMSE = %.4f, cov of 0 = %.2f",
                         results$null_M1$rmse_f12, results$null_M1$cov_f12))
    contour(x1, x2, surfB$post_mean, add = TRUE, lwd = 1.2, labcex = 0.9)
  }

  resid_surf <- surfD$post_mean - truth_int
  r_lim <- max(abs(resid_surf))
  r_breaks <- seq(-r_lim, r_lim, length.out = 51)
  image(x1, x2, resid_surf, col = pal, breaks = r_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("RESIDUAL: posterior mean - truth  (Cell D)\nrange [%.2f, %.2f]",
                       min(resid_surf), max(resid_surf)))
  contour(x1, x2, resid_surf, add = TRUE, lwd = 1.2, labcex = 0.9)
  dev.off()
  cat(sprintf("  Surface plots:  %s\n", outfile))
}


# ============================================================
# plot_surfaces_3d_v2(): Fahrmeir/Kneib-style 3D wireframe surfaces
#
# Produces 3D perspective plots (persp) for:
#   Page 1: TRUE f_12 surface (INT DGP)
#   Page 2: Posterior mean f_12 (Cell D: INT + M1)
#   Page 3: Posterior mean + upper/lower 95% CI bands (3 surfaces)
#   Page 4: Posterior mean f_12 (Cell B: NULL + M1) — near-zero
#   Page 5: Residual (fitted − truth)
# ============================================================
plot_surfaces_3d_v2 <- function(res_obj, outfile) {
  s <- res_obj$settings
  results <- res_obj$results
  truth_int <- res_obj$truth_grid_int$f12

  surfD <- results$int_M1$surface12
  surfB <- results$null_M1$surface12

  pdf(outfile, width = 9, height = 7.5)

  if (is.null(surfD)) {
    plot.new(); text(0.5, 0.5, "No interaction surface stored")
    dev.off(); return(invisible(NULL))
  }

  x1 <- surfD$x1; x2 <- surfD$x2

  # Common z-range for truth + fitted
  z_all <- c(as.vector(truth_int), as.vector(surfD$post_mean))
  zlim <- c(min(z_all), max(z_all)) * 1.1

  # Viewing angle (matches Fahrmeir/Kneib style)
  theta <- -30; phi <- 25

  # Color function for wireframe
  make_facet_colors <- function(z_mat, pal = "Blue-Red 2") {
    nrz <- nrow(z_mat); ncz <- ncol(z_mat)
    # Average of 4 corners per facet
    z_facet <- (z_mat[-1, -1] + z_mat[-1, -ncz] +
                  z_mat[-nrz, -1] + z_mat[-nrz, -ncz]) / 4
    n_col <- 50
    cols <- hcl.colors(n_col, pal)
    z_range <- range(z_facet, na.rm = TRUE)
    if (diff(z_range) == 0) return(rep(cols[n_col %/% 2], length(z_facet)))
    idx <- findInterval(z_facet, seq(z_range[1], z_range[2], length.out = n_col + 1),
                        all.inside = TRUE)
    cols[idx]
  }

  par(mar = c(1, 1, 3, 1))

  # Page 1: TRUTH
  persp(x1, x2, truth_int, zlim = zlim,
        theta = theta, phi = phi, expand = 0.6,
        col = make_facet_colors(truth_int),
        shade = 0.3, border = NA, ltheta = 120,
        xlab = "X1", ylab = "X2", zlab = "f_12",
        main = sprintf("TRUE f_12  (p=%d, phase %d)\nrange [%.2f, %.2f]",
                       s$p, s$phase, min(truth_int), max(truth_int)),
        ticktype = "detailed", cex.main = 1.3)

  # Page 2: Posterior mean (Cell D)
  persp(x1, x2, surfD$post_mean, zlim = zlim,
        theta = theta, phi = phi, expand = 0.6,
        col = make_facet_colors(surfD$post_mean),
        shade = 0.3, border = NA, ltheta = 120,
        xlab = "X1", ylab = "X2", zlab = "f_12",
        main = sprintf("POSTERIOR MEAN f_12  (Cell D: INT + M1)\nRMSE = %.3f, cov = %.2f",
                       results$int_M1$rmse_f12, results$int_M1$cov_f12),
        ticktype = "detailed", cex.main = 1.3)

  # Page 3: Posterior mean + CI bands (wireframe with transparency)
  # Mean as solid colored surface, lo/hi as gray wireframes
  pmat <- persp(x1, x2, surfD$post_mean, zlim = zlim,
                theta = theta, phi = phi, expand = 0.6,
                col = make_facet_colors(surfD$post_mean),
                shade = 0.3, border = NA, ltheta = 120,
                xlab = "X1", ylab = "X2", zlab = "f_12",
                main = sprintf("POSTERIOR MEAN + 95%% CI  (Cell D)\nmean width = %.3f",
                               mean(surfD$band_width)),
                ticktype = "detailed", cex.main = 1.3)
  # Overlay upper band wireframe
  par(new = TRUE)
  persp(x1, x2, surfD$post_hi, zlim = zlim,
        theta = theta, phi = phi, expand = 0.6,
        col = NA, border = adjustcolor("gray50", alpha.f = 0.3),
        lwd = 0.3, shade = NA,
        xlab = "", ylab = "", zlab = "",
        axes = FALSE, box = FALSE)
  par(new = TRUE)
  persp(x1, x2, surfD$post_lo, zlim = zlim,
        theta = theta, phi = phi, expand = 0.6,
        col = NA, border = adjustcolor("gray50", alpha.f = 0.3),
        lwd = 0.3, shade = NA,
        xlab = "", ylab = "", zlab = "",
        axes = FALSE, box = FALSE)

  # Page 4: NULL + M1 (Cell B) — should be near-flat
  if (!is.null(surfB)) {
    zlim_null <- range(c(as.vector(surfB$post_lo), as.vector(surfB$post_hi))) * 1.2
    if (diff(zlim_null) < 0.01) zlim_null <- c(-0.1, 0.1)
    persp(x1, x2, surfB$post_mean, zlim = zlim_null,
          theta = theta, phi = phi, expand = 0.6,
          col = make_facet_colors(surfB$post_mean),
          shade = 0.3, border = NA, ltheta = 120,
          xlab = "X1", ylab = "X2", zlab = "f_12",
          main = sprintf("NULL + M1 (Cell B, truth=0)\nRMSE = %.4f, cov of 0 = %.2f",
                         results$null_M1$rmse_f12, results$null_M1$cov_f12),
          ticktype = "detailed", cex.main = 1.3)
  }

  # Page 5: Residual
  resid_surf <- surfD$post_mean - truth_int
  zlim_res <- max(abs(resid_surf)) * c(-1, 1) * 1.1
  persp(x1, x2, resid_surf, zlim = zlim_res,
        theta = theta, phi = phi, expand = 0.6,
        col = make_facet_colors(resid_surf),
        shade = 0.3, border = NA, ltheta = 120,
        xlab = "X1", ylab = "X2", zlab = "Residual",
        main = sprintf("RESIDUAL: fitted - truth  (Cell D)\nrange [%.2f, %.2f]",
                       min(resid_surf), max(resid_surf)),
        ticktype = "detailed", cex.main = 1.3)

  dev.off()
  cat(sprintf("  3D surfaces:    %s\n", outfile))
}


# ============================================================
# plot_diagnostics_v2(): trace plots, posterior densities, ACFs
#
# For each cell in the 2x2: one multi-page PDF with:
#   Page 1: Trace plots (sigma2, tau2, rho, tau2_s_main_1..p, tau2_s_int_1..K)
#   Page 2: Posterior densities with truth lines
#   Page 3: ACF plots (for ESS)
# ============================================================
plot_diagnostics_v2 <- function(res_obj, outfile) {
  s <- res_obj$settings
  results <- res_obj$results
  sim_null <- res_obj$sim_null
  sim_int  <- res_obj$sim_int

  cell_list <- list(
    list(name = "Cell A: NULL + M0", obj = results$null_M0, sim = sim_null),
    list(name = "Cell B: NULL + M1", obj = results$null_M1, sim = sim_null),
    list(name = "Cell C: INT  + M0", obj = results$int_M0,  sim = sim_int),
    list(name = "Cell D: INT  + M1", obj = results$int_M1,  sim = sim_int)
  )

  true_sigma2 <- sim_null$true_sigma2
  true_tau2   <- sim_null$true_tau2
  true_rho    <- sim_null$true_rho

  pdf(outfile, width = 11, height = 8.5)

  for (cell in cell_list) {
    ch <- cell$obj$chains
    if (is.null(ch)) next

    p_ <- s$p
    n_main <- ncol(ch$tau2_s_main)
    n_int  <- if (!is.null(ch$tau2_s_int)) ncol(ch$tau2_s_int) else 0

    # Build parameter list: name, samples, truth (NA if unknown)
    params <- list(
      list(name = "sigma2", samp = ch$sigma2, truth = true_sigma2),
      list(name = "tau2",   samp = ch$tau2,   truth = true_tau2),
      list(name = "rho",    samp = ch$rho,    truth = true_rho)
    )
    for (j in seq_len(n_main)) {
      params[[length(params) + 1]] <- list(
        name = sprintf("tau2_s_%d", j),
        samp = ch$tau2_s_main[, j],
        truth = NA    # no single truth for smoothing variance
      )
    }
    for (k in seq_len(n_int)) {
      params[[length(params) + 1]] <- list(
        name = sprintf("tau2_s_int_%d", k),
        samp = ch$tau2_s_int[, k],
        truth = NA
      )
    }

    np <- length(params)

    # --- Page: Trace plots ---
    nr <- ceiling(np / 2); nc <- 2
    par(mfrow = c(min(nr, 4), nc), mar = c(3.5, 4, 3, 1),
        oma = c(0, 0, 3, 0), cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)

    page_count <- 0
    for (i in seq_along(params)) {
      if (page_count >= 8) {
        mtext(sprintf("Trace plots — %s  (p=%d, phase %d)", cell$name, s$p, s$phase),
              outer = TRUE, cex = 1.3, font = 2)
        par(mfrow = c(min(4, ceiling((np - i + 1) / 2)), 2), mar = c(3.5, 4, 3, 1),
            oma = c(0, 0, 3, 0))
        page_count <- 0
      }
      pm <- params[[i]]
      plot(pm$samp, type = "l", main = pm$name, ylab = pm$name, xlab = "Iteration",
           col = "gray30", lwd = 0.5)
      if (!is.na(pm$truth)) abline(h = pm$truth, col = "red", lty = 2, lwd = 2)
      abline(h = mean(pm$samp), col = "steelblue", lty = 3, lwd = 1.5)
      page_count <- page_count + 1
    }
    mtext(sprintf("Trace plots — %s  (p=%d, phase %d)", cell$name, s$p, s$phase),
          outer = TRUE, cex = 1.3, font = 2)

    # --- Page: Posterior densities ---
    par(mfrow = c(min(nr, 4), nc), mar = c(3.5, 4, 3, 1),
        oma = c(0, 0, 3, 0))
    page_count <- 0
    for (i in seq_along(params)) {
      if (page_count >= 8) {
        mtext(sprintf("Posterior densities — %s  (p=%d, phase %d)",
                      cell$name, s$p, s$phase),
              outer = TRUE, cex = 1.3, font = 2)
        par(mfrow = c(min(4, ceiling((np - i + 1) / 2)), 2), mar = c(3.5, 4, 3, 1),
            oma = c(0, 0, 3, 0))
        page_count <- 0
      }
      pm <- params[[i]]
      d_post <- density(pm$samp)
      plot(d_post, main = pm$name, xlab = pm$name, lwd = 2, col = "steelblue")
      if (!is.na(pm$truth)) {
        abline(v = pm$truth, col = "red", lty = 2, lwd = 2)
        legend("topright", c("posterior", "truth"),
               col = c("steelblue", "red"), lty = c(1, 2), lwd = 2,
               bty = "n", cex = 0.85)
      }
      # Posterior mean line
      abline(v = mean(pm$samp), col = "gray40", lty = 3, lwd = 1.5)
      # Credible interval
      ci <- quantile(pm$samp, c(0.025, 0.975))
      abline(v = ci, col = "gray60", lty = 3, lwd = 1)
      page_count <- page_count + 1
    }
    mtext(sprintf("Posterior densities — %s  (p=%d, phase %d)",
                  cell$name, s$p, s$phase),
          outer = TRUE, cex = 1.3, font = 2)

    # --- Page: ACF plots (for ESS) ---
    par(mfrow = c(min(nr, 4), nc), mar = c(3.5, 4, 3, 1),
        oma = c(0, 0, 3, 0))
    page_count <- 0
    for (i in seq_along(params)) {
      if (page_count >= 8) {
        mtext(sprintf("ACF — %s  (p=%d, phase %d)", cell$name, s$p, s$phase),
              outer = TRUE, cex = 1.3, font = 2)
        par(mfrow = c(min(4, ceiling((np - i + 1) / 2)), 2), mar = c(3.5, 4, 3, 1),
            oma = c(0, 0, 3, 0))
        page_count <- 0
      }
      pm <- params[[i]]
      acf(pm$samp, main = sprintf("ACF: %s  (ESS~%.0f)",
                                    pm$name,
                                    length(pm$samp) / (1 + 2 * sum(pmax(0, acf(pm$samp, plot = FALSE)$acf[-1])))),
          lag.max = 50, col = "steelblue", lwd = 2)
      page_count <- page_count + 1
    }
    mtext(sprintf("ACF — %s  (p=%d, phase %d)", cell$name, s$p, s$phase),
          outer = TRUE, cex = 1.3, font = 2)
  }

  dev.off()
  cat(sprintf("  Diagnostics:    %s\n", outfile))
}


# ============================================================
# Driver
# ============================================================
if (!interactive()) {
  source("spatial_utils.R")
  source("ls_basis.R")
  source("ls_interaction.R")
  source("gibbs_interaction.R")

  cat("==============================================================\n")
  cat(" UNIFIED 2x2 EXPERIMENT (centered Sim4 DGP)\n")
  cat(" p in {2,3,4}  x  phase in {1=X1xX2 only, 2=all pairwise}\n")
  cat("==============================================================\n")

  N_OBS  <- 300
  M_KNOTS <- 8
  N_ITER <- 3000
  N_BURN <- 1000
  C_INT  <- 1.5

  t_total <- proc.time()

  all_res <- list()
  for (p_ in c(2, 3, 4)) {
    for (phase_ in c(1, 2)) {
      tag <- sprintf("p%d_phase%d", p_, phase_)
      cat(sprintf("\n=================  %s  =================\n", tag))
      res <- run_2x2(p = p_, phase = phase_,
                     n = N_OBS, M = M_KNOTS,
                     n_iter = N_ITER, n_burn = N_BURN, c_int = C_INT)
      all_res[[tag]] <- res

      plot_marginals_v2(res,
        outfile = sprintf("interaction_2x2_v2_%s_marginals.pdf", tag))
      plot_surfaces_v2(res,
        outfile = sprintf("interaction_2x2_v2_%s_surfaces.pdf", tag))
      plot_diagnostics_v2(res,
        outfile = sprintf("interaction_2x2_v2_%s_diagnostics.pdf", tag))
    }
  }

  total_elapsed <- as.numeric((proc.time() - t_total)["elapsed"])
  cat(sprintf("\n##########  TOTAL: %.1f sec (%.1f min)  ##########\n",
              total_elapsed, total_elapsed / 60))

  saveRDS(list(all = all_res, total_sec = total_elapsed),
          file = "interaction_2x2_v2_results.rds")

  df_all <- do.call(rbind, lapply(all_res, results_to_dataframe_v2))
  rownames(df_all) <- NULL
  write.csv(df_all, "interaction_2x2_v2_summary.csv", row.names = FALSE)
  cat("CSV saved to interaction_2x2_v2_summary.csv\n")
  cat("RDS saved to interaction_2x2_v2_results.rds\n")
}
