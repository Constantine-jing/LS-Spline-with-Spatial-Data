# ============================================================
# run_interaction_2x2.R
#
# 2x2 FACTORIAL COMPARISON: {DGP-NULL, DGP-INT} x {Model M0, Model M1}
#
# Same simulation setting in both rows; same sampler code path in both columns.
# DGP varies via `c_int` (true interaction amplitude); Model varies via
# whether interaction blocks are included.
#
# Runs both p=2 (one interaction) and p=3 (three interactions) experiments.
#
# All four cells per experiment use the same data seed (so M0 vs M1 see
# identical y) and same MCMC seed (so RNG noise is matched).
#
# Output:
#   - Console: per-cell summaries + 2x2 comparison tables
#   - PDF:     traces and key diagnostics
#   - RDS:     full results object
# ============================================================

# Required sources (called in driver block at the bottom)

# ============================================================
# Helpers
# ============================================================

# WAIC computation from per-iter log-likelihood matrix
# logLik_mat: n_keep x n
compute_waic <- function(logLik_mat) {
  lpd  <- sum(log(colMeans(exp(logLik_mat))))   # log pointwise predictive density
  pwaic <- sum(apply(logLik_mat, 2, var))       # effective number of parameters
  waic <- -2 * (lpd - pwaic)
  list(waic = waic, lpd = lpd, p_waic = pwaic)
}

# Per-iteration Gaussian log-likelihood given residuals and (sigma2, tau2, rho)
# This is the MARGINAL likelihood after integrating out b, evaluated pointwise.
# For pointwise WAIC we use the conditional likelihood: y_i | eta, b ~ N(H_i eta + b_i, tau2)
# But b is integrated out by the collapsed sampler, so we use the marginal:
#   y | eta ~ N(H eta, sigma2*R + tau2*I)  -- joint, not pointwise.
# For pointwise contributions, the cleanest is: y_i | eta, sigma2, tau2 ~ N((H eta)_i, sigma2 + tau2)
# treating b_i marginally. This is approximate but standard for WAIC under collapsed samplers.
pointwise_loglik <- function(y, fitted_mat, sigma2_vec, tau2_vec) {
  # y: n; fitted_mat: n_keep x n; sigma2_vec, tau2_vec: n_keep
  n_keep <- nrow(fitted_mat)
  n      <- length(y)
  ll_mat <- matrix(NA, n_keep, n)
  for (k in seq_len(n_keep)) {
    sd_k <- sqrt(sigma2_vec[k] + tau2_vec[k])
    ll_mat[k, ] <- dnorm(y, mean = fitted_mat[k, ], sd = sd_k, log = TRUE)
  }
  ll_mat
}


# ============================================================
# fit_one_cell(): run sampler for a single (DGP, Model) cell
# ============================================================
fit_one_cell <- function(y, H, D, col_map_main, K_main_list,
                          col_map_int, K_int_list,
                          n_iter, n_burn, mcmc_seed,
                          true_nu = 1.5, verbose = FALSE) {
  set.seed(mcmc_seed)
  t0 <- proc.time()
  gs <- gibbs_interaction_sampler(
    y = y, H = H, D = D, nu = true_nu,
    col_map_main = col_map_main,
    K_main_list  = K_main_list,
    col_map_int  = col_map_int,
    K_int_list   = K_int_list,
    n_iter = n_iter, n_burn = n_burn,
    verbose = verbose
  )
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  list(gs = gs, elapsed = elapsed)
}


# ============================================================
# build_design_p2(): build H + col_maps for p=2 model
# include_int = FALSE gives M0; TRUE gives M1
# ============================================================
build_design_p2 <- function(X1, X2, M, include_int = TRUE) {
  obj1 <- ls_build_one_full(X1, M = M)
  obj2 <- ls_build_one_full(X2, M = M)
  W1 <- obj1$W; W2 <- obj2$W
  T1 <- obj1$T; T2 <- obj2$T

  K1_raw <- build_rw2_penalty_1d(M)
  K2_raw <- build_rw2_penalty_1d(M)
  K1_beta <- t(T1) %*% K1_raw %*% T1
  K2_beta <- t(T2) %*% K2_raw %*% T2

  d1 <- ncol(W1); d2 <- ncol(W2)

  if (include_int) {
    int12 <- ls_build_interaction(obj1, obj2)
    W12 <- int12$W_uv; K12 <- int12$K_uv
    d12 <- ncol(W12)
    H <- cbind(1, W1, W2, W12)
    col_map_main <- list(seq(2, d1 + 1), seq(d1 + 2, d1 + d2 + 1))
    col_map_int  <- list(seq(d1 + d2 + 2, d1 + d2 + d12 + 1))
    K_main_list <- list(K1_beta, K2_beta)
    K_int_list  <- list(K12)
    int_recipes <- list(int12)
  } else {
    H <- cbind(1, W1, W2)
    col_map_main <- list(seq(2, d1 + 1), seq(d1 + 2, d1 + d2 + 1))
    col_map_int  <- list()
    K_main_list <- list(K1_beta, K2_beta)
    K_int_list  <- list()
    int_recipes <- list()
  }

  col_map_main <- lapply(col_map_main, function(x) x - 1L)
  col_map_int  <- lapply(col_map_int,  function(x) x - 1L)

  list(H = H, col_map_main = col_map_main, col_map_int = col_map_int,
       K_main_list = K_main_list, K_int_list = K_int_list,
       objs = list(obj1, obj2), int_recipes = int_recipes,
       d_main = c(d1, d2))
}


# ============================================================
# build_design_p3(): build H + col_maps for p=3 model
# 3 main effects, optionally 3 pairwise interactions (1-2, 1-3, 2-3)
# ============================================================
build_design_p3 <- function(X1, X2, X3, M, include_int = TRUE) {
  obj1 <- ls_build_one_full(X1, M = M)
  obj2 <- ls_build_one_full(X2, M = M)
  obj3 <- ls_build_one_full(X3, M = M)
  W1 <- obj1$W; W2 <- obj2$W; W3 <- obj3$W
  T1 <- obj1$T; T2 <- obj2$T; T3 <- obj3$T

  K_raw <- build_rw2_penalty_1d(M)
  K1_beta <- t(T1) %*% K_raw %*% T1
  K2_beta <- t(T2) %*% K_raw %*% T2
  K3_beta <- t(T3) %*% K_raw %*% T3

  d1 <- ncol(W1); d2 <- ncol(W2); d3 <- ncol(W3)

  if (include_int) {
    int12 <- ls_build_interaction(obj1, obj2)
    int13 <- ls_build_interaction(obj1, obj3)
    int23 <- ls_build_interaction(obj2, obj3)
    W12 <- int12$W_uv; W13 <- int13$W_uv; W23 <- int23$W_uv
    K12 <- int12$K_uv; K13 <- int13$K_uv; K23 <- int23$K_uv
    d12 <- ncol(W12); d13 <- ncol(W13); d23 <- ncol(W23)

    H <- cbind(1, W1, W2, W3, W12, W13, W23)

    # 1-indexed column ranges
    s <- 1L  # last filled column (intercept = col 1)
    r_W1  <- seq(s + 1, s + d1);  s <- s + d1
    r_W2  <- seq(s + 1, s + d2);  s <- s + d2
    r_W3  <- seq(s + 1, s + d3);  s <- s + d3
    r_W12 <- seq(s + 1, s + d12); s <- s + d12
    r_W13 <- seq(s + 1, s + d13); s <- s + d13
    r_W23 <- seq(s + 1, s + d23); s <- s + d23

    col_map_main <- list(r_W1, r_W2, r_W3)
    col_map_int  <- list(r_W12, r_W13, r_W23)
    K_main_list  <- list(K1_beta, K2_beta, K3_beta)
    K_int_list   <- list(K12, K13, K23)
    int_recipes  <- list(int12, int13, int23)
  } else {
    H <- cbind(1, W1, W2, W3)
    s <- 1L
    r_W1 <- seq(s + 1, s + d1); s <- s + d1
    r_W2 <- seq(s + 1, s + d2); s <- s + d2
    r_W3 <- seq(s + 1, s + d3); s <- s + d3

    col_map_main <- list(r_W1, r_W2, r_W3)
    col_map_int  <- list()
    K_main_list  <- list(K1_beta, K2_beta, K3_beta)
    K_int_list   <- list()
    int_recipes  <- list()
  }

  col_map_main <- lapply(col_map_main, function(x) x - 1L)
  col_map_int  <- lapply(col_map_int,  function(x) x - 1L)

  list(H = H, col_map_main = col_map_main, col_map_int = col_map_int,
       K_main_list = K_main_list, K_int_list = K_int_list,
       objs = list(obj1, obj2, obj3), int_recipes = int_recipes,
       d_main = c(d1, d2, d3))
}


# ============================================================
# simulate_p2(): generate data with controllable interaction amplitude
# ============================================================
simulate_p2 <- function(n, c_int, data_seed,
                        sigma2 = 0.5, rho = 0.3, tau2 = 0.25, nu = 1.5) {
  set.seed(data_seed)
  locs <- matrix(runif(2 * n), n, 2)
  D <- as.matrix(dist(locs))
  X1 <- runif(n); X2 <- runif(n)

  true_f1  <- sin(2 * pi * X1)
  true_f2  <- (2 * X2 - 1)^2
  true_f12 <- c_int * sin(pi * X1) * cos(pi * X2)

  R_true <- matern_cor(D, rho = rho, nu = nu)
  L_R <- chol(R_true + diag(1e-8, n))
  b <- as.vector(t(L_R) %*% rnorm(n)) * sqrt(sigma2)
  eps <- rnorm(n, 0, sqrt(tau2))

  y <- true_f1 + true_f2 + true_f12 + b + eps
  list(y = y, D = D, X1 = X1, X2 = X2, locs = locs,
       true_f1 = true_f1, true_f2 = true_f2, true_f12 = true_f12,
       true_b = b, true_sigma2 = sigma2, true_tau2 = tau2,
       true_rho = rho, true_nu = nu, c_int = c_int)
}


# ============================================================
# simulate_p3(): p=3 version with one pairwise interaction (1-2 only)
# (X3 has only a main effect; interactions 1-3 and 2-3 are zero)
# ============================================================
simulate_p3 <- function(n, c_int, data_seed,
                        sigma2 = 0.5, rho = 0.3, tau2 = 0.25, nu = 1.5) {
  set.seed(data_seed)
  locs <- matrix(runif(2 * n), n, 2)
  D <- as.matrix(dist(locs))
  X1 <- runif(n); X2 <- runif(n); X3 <- runif(n)

  true_f1  <- sin(2 * pi * X1)
  true_f2  <- (2 * X2 - 1)^2
  true_f3  <- 0.8 * (X3 - 0.5)             # mild linear-ish effect
  true_f12 <- c_int * sin(pi * X1) * cos(pi * X2)
  true_f13 <- rep(0, n)
  true_f23 <- rep(0, n)

  R_true <- matern_cor(D, rho = rho, nu = nu)
  L_R <- chol(R_true + diag(1e-8, n))
  b <- as.vector(t(L_R) %*% rnorm(n)) * sqrt(sigma2)
  eps <- rnorm(n, 0, sqrt(tau2))

  y <- true_f1 + true_f2 + true_f3 + true_f12 + true_f13 + true_f23 + b + eps
  list(y = y, D = D, X1 = X1, X2 = X2, X3 = X3, locs = locs,
       true_f1 = true_f1, true_f2 = true_f2, true_f3 = true_f3,
       true_f12 = true_f12, true_f13 = true_f13, true_f23 = true_f23,
       true_b = b, true_sigma2 = sigma2, true_tau2 = tau2,
       true_rho = rho, true_nu = nu, c_int = c_int)
}


# ============================================================
# summarise_cell(): compute RMSEs, fitted values, WAIC for one fitted cell
#
# If `grid_main` is supplied (length-p list of x-grids), also computes
# posterior samples of f_j on the grid and returns mean + 95% bands.
# If `grid_int` is supplied (list with $X1, $X2), computes f12 on the
# tensor grid (length(X1) x length(X2)).
# ============================================================
summarise_cell <- function(fit, design, sim, p,
                            grid_main = NULL, grid_int = NULL) {
  gs <- fit$gs
  H  <- design$H
  y  <- sim$y
  d_main <- design$d_main
  has_int <- length(design$col_map_int) > 0

  # eta column indices (1-indexed in eta vector, same as H)
  idx_mu <- 1
  starts <- cumsum(c(2, d_main))[1:p]   # start of each main block in eta
  ends   <- starts + d_main - 1
  idx_main_list <- mapply(seq, starts, ends, SIMPLIFY = FALSE)

  # Posterior mean eta
  eta_mean <- colMeans(gs$eta_samples)

  # Recover main-effect fitted values at training X
  f_hat <- list()
  for (j in seq_len(p)) {
    Wj <- design$objs[[j]]$W
    f_hat[[j]] <- as.vector(Wj %*% eta_mean[idx_main_list[[j]]])
  }
  names(f_hat) <- paste0("f", seq_len(p), "_hat")

  # Truth main effects from sim
  f_true <- list()
  for (j in seq_len(p)) {
    f_true[[j]] <- sim[[paste0("true_f", j)]]
  }
  rmse_main <- sapply(seq_len(p), function(j) sqrt(mean((f_hat[[j]] - f_true[[j]])^2)))
  names(rmse_main) <- paste0("RMSE_f", seq_len(p))

  # ---- Marginal grid posteriors (NEW) ----
  marginal <- NULL
  if (!is.null(grid_main)) {
    marginal <- vector("list", p)
    for (j in seq_len(p)) {
      x_grid <- grid_main[[j]]
      Wj_grid <- design$objs[[j]]$design_new(x_grid, type = "W", clip = TRUE)
      beta_j_post <- gs$eta_samples[, idx_main_list[[j]], drop = FALSE]
      f_grid_samp <- beta_j_post %*% t(Wj_grid)   # n_keep x length(x_grid)
      marginal[[j]] <- list(
        x        = x_grid,
        post_mean = colMeans(f_grid_samp),
        post_lo   = apply(f_grid_samp, 2, quantile, 0.025),
        post_hi   = apply(f_grid_samp, 2, quantile, 0.975)
      )
    }
    names(marginal) <- paste0("f", seq_len(p))
  }

  # Interaction fitted (if model includes it) -- training data
  f12_hat <- NULL; f12_lo <- NULL; f12_hi <- NULL; rmse_f12 <- NA; cov_f12 <- NA
  surface12 <- NULL
  if (has_int) {
    # First interaction block is always X1xX2 by our design convention
    int_idx <- design$col_map_int[[1]] + 1L   # 0-indexed -> 1-indexed
    beta12_post <- gs$eta_samples[, int_idx, drop = FALSE]
    W12 <- design$int_recipes[[1]]$W_uv
    f12_samp <- beta12_post %*% t(W12)
    f12_hat <- colMeans(f12_samp)
    f12_lo  <- apply(f12_samp, 2, quantile, 0.025)
    f12_hi  <- apply(f12_samp, 2, quantile, 0.975)
    rmse_f12 <- sqrt(mean((f12_hat - sim$true_f12)^2))
    # If true_f12 is exactly 0 (NULL DGP), report band coverage of 0
    cov_f12 <- mean(f12_lo <= sim$true_f12 & sim$true_f12 <= f12_hi)

    # ---- Interaction surface on regular grid (NEW) ----
    if (!is.null(grid_int)) {
      W12_grid <- ls_interaction_design_new(grid_int$X1_long, grid_int$X2_long,
                                            design$int_recipes[[1]]$recipe,
                                            clip = TRUE)
      surf_samp <- beta12_post %*% t(W12_grid)   # n_keep x (g x g)
      surface12 <- list(
        x1 = grid_int$x1,
        x2 = grid_int$x2,
        post_mean = matrix(colMeans(surf_samp),     length(grid_int$x1), length(grid_int$x2)),
        post_lo   = matrix(apply(surf_samp, 2, quantile, 0.025),
                           length(grid_int$x1), length(grid_int$x2)),
        post_hi   = matrix(apply(surf_samp, 2, quantile, 0.975),
                           length(grid_int$x1), length(grid_int$x2))
      )
      surface12$band_width <- surface12$post_hi - surface12$post_lo
    }
  }

  # Fitted y at each kept iteration (for WAIC)
  fitted_mat <- gs$eta_samples %*% t(H)   # n_keep x n
  ll_mat <- pointwise_loglik(y, fitted_mat, gs$sigma2_samples, gs$tau2_samples)
  waic_obj <- compute_waic(ll_mat)

  # Predictive RMSE (in-sample)
  y_hat <- colMeans(fitted_mat)
  rmse_y <- sqrt(mean((y_hat - y)^2))

  list(
    rmse_main = rmse_main,
    rmse_f12  = rmse_f12,
    cov_f12   = cov_f12,
    rmse_y    = rmse_y,
    waic      = waic_obj$waic,
    p_waic    = waic_obj$p_waic,
    sigma2_mean = mean(gs$sigma2_samples),
    tau2_mean   = mean(gs$tau2_samples),
    rho_mean    = mean(gs$rho_samples),
    accept = if (!is.null(gs$accept_rate)) gs$accept_rate else c(sigma2 = NA, tau2 = NA, rho = NA),
    elapsed = fit$elapsed,
    n_keep = nrow(gs$eta_samples),
    f12_hat = f12_hat, f12_lo = f12_lo, f12_hi = f12_hi,
    marginal  = marginal,    # NEW: list of f_j grid posteriors
    surface12 = surface12    # NEW: f12 grid posterior surface (M1 cells only)
  )
}


# ============================================================
# results_to_dataframe(): flatten 4-cell results into a tidy data.frame
# ============================================================
results_to_dataframe <- function(results, p, settings) {
  cells <- list(
    list(name = "A", dgp = "NULL", model = "M0", obj = results$null_M0),
    list(name = "B", dgp = "NULL", model = "M1", obj = results$null_M1),
    list(name = "C", dgp = "INT",  model = "M0", obj = results$int_M0),
    list(name = "D", dgp = "INT",  model = "M1", obj = results$int_M1)
  )

  rows <- lapply(cells, function(cell) {
    rm <- cell$obj$rmse_main
    row <- data.frame(
      p             = p,
      n             = settings$n,
      M             = settings$M,
      n_iter        = settings$n_iter,
      n_burn        = settings$n_burn,
      c_int_true    = settings$c_int,
      cell          = cell$name,
      DGP           = cell$dgp,
      Model         = cell$model,
      RMSE_f1       = rm[1],
      RMSE_f2       = if (length(rm) >= 2) rm[2] else NA,
      RMSE_f3       = if (length(rm) >= 3) rm[3] else NA,
      RMSE_f12      = cell$obj$rmse_f12,
      cov_f12       = cell$obj$cov_f12,
      RMSE_y        = cell$obj$rmse_y,
      sigma2_mean   = cell$obj$sigma2_mean,
      tau2_mean     = cell$obj$tau2_mean,
      rho_mean      = cell$obj$rho_mean,
      WAIC          = cell$obj$waic,
      p_WAIC        = cell$obj$p_waic,
      accept_sigma2 = cell$obj$accept[1],
      accept_tau2   = cell$obj$accept[2],
      accept_rho    = cell$obj$accept[3],
      elapsed_sec   = cell$obj$elapsed,
      stringsAsFactors = FALSE
    )
    rownames(row) <- NULL
    row
  })
  do.call(rbind, rows)
}


# ============================================================
# plot_marginals_2x2(): full-page marginal recovery plots with bands
#
# For each main effect f_j (j = 1..p), produces a single page with a 2x2
# layout of cells: (NULL+M0, NULL+M1, INT+M0, INT+M1).
# Each panel: truth (solid black) + posterior mean (dashed red) + 95% band (shaded).
# ============================================================
plot_marginals_2x2 <- function(res_obj, p, outfile) {
  pdf(outfile, width = 11, height = 8.5)

  results <- res_obj$results
  truth_null <- res_obj$truth_grid_null
  truth_int  <- res_obj$truth_grid_int

  cells <- list(
    list(name = "NULL + M0 (Cell A)", obj = results$null_M0, truth = truth_null),
    list(name = "NULL + M1 (Cell B)", obj = results$null_M1, truth = truth_null),
    list(name = "INT  + M0 (Cell C)", obj = results$int_M0,  truth = truth_int),
    list(name = "INT  + M1 (Cell D)", obj = results$int_M1,  truth = truth_int)
  )

  for (j in seq_len(p)) {
    fname <- paste0("f", j)

    # Compute consistent y-range across all 4 cells for this main effect
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
      if (is.null(m)) {
        plot.new(); next
      }
      x <- m$x
      true_y <- cell$truth[[fname]]

      plot(x, true_y, type = "n", ylim = ylim,
           xlab = paste0("X", j), ylab = paste0("f", j, "(X", j, ")"),
           main = sprintf("%s   RMSE = %.3f", cell$name,
                          cell$obj$rmse_main[j]))
      # Shaded 95% band
      polygon(c(x, rev(x)), c(m$post_lo, rev(m$post_hi)),
              col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
      # Posterior mean
      lines(x, m$post_mean, col = "red", lwd = 2.5, lty = 2)
      # Truth
      lines(x, true_y, col = "black", lwd = 2.5)

      legend("topright",
             c("truth", "posterior mean", "95% CI"),
             col = c("black", "red", adjustcolor("steelblue", alpha.f = 0.5)),
             lwd = c(2.5, 2.5, 8), lty = c(1, 2, 1),
             bty = "n", cex = 0.95)
      grid()
    }
    mtext(sprintf("Marginal recovery: f%d  (p = %d)", j, p),
          outer = TRUE, cex = 1.5, font = 2)
  }

  dev.off()
  cat(sprintf("  Marginal plots written to: %s\n", outfile))
}


# ============================================================
# plot_interaction_surfaces(): full-page f12 surface plots
#
# Produces a multi-page PDF:
#   Page 1: TRUTH surface for INT DGP
#   Page 2: Posterior mean f12 from INT + M1 (Cell D)
#   Page 3: Pointwise band width (uncertainty heatmap) for Cell D
#   Page 4: Posterior mean f12 from NULL + M1 (Cell B) -- should be near 0
#   Page 5: Difference (Cell D fitted - truth)
# ============================================================
plot_interaction_surfaces <- function(res_obj, p, outfile) {
  pdf(outfile, width = 9, height = 8)

  truth_int  <- res_obj$truth_grid_int$f12
  results    <- res_obj$results

  surfD <- results$int_M1$surface12     # Cell D
  surfB <- results$null_M1$surface12    # Cell B (NULL+M1, should be ~0)

  if (is.null(surfD)) {
    cat("  Warning: no surface12 stored for Cell D; skipping interaction plots.\n")
    dev.off()
    return(invisible(NULL))
  }

  x1 <- surfD$x1; x2 <- surfD$x2

  # Common color scale across truth + fitted + null
  z_all <- c(as.vector(truth_int), as.vector(surfD$post_mean),
             if (!is.null(surfB)) as.vector(surfB$post_mean) else numeric(0))
  z_lim <- max(abs(z_all))
  z_breaks <- seq(-z_lim, z_lim, length.out = 51)
  pal <- hcl.colors(50, "Blue-Red 2")

  par(mar = c(5, 5, 4, 6), cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)

  # Page 1: TRUTH
  image(x1, x2, truth_int, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("TRUE f_12(X1,X2)  (INT DGP, p = %d)\nrange [%.2f, %.2f]",
                       p, min(truth_int), max(truth_int)))
  contour(x1, x2, truth_int, add = TRUE, lwd = 1.2, labcex = 0.9)

  # Page 2: Cell D posterior mean
  image(x1, x2, surfD$post_mean, col = pal, breaks = z_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("POSTERIOR MEAN f_12  (Cell D: INT + M1)\nRMSE = %.3f, coverage of truth = %.2f",
                       results$int_M1$rmse_f12, results$int_M1$cov_f12))
  contour(x1, x2, surfD$post_mean, add = TRUE, lwd = 1.2, labcex = 0.9)

  # Page 3: uncertainty (band width)
  pal_unc <- hcl.colors(50, "YlOrRd", rev = TRUE)
  image(x1, x2, surfD$band_width, col = pal_unc,
        xlab = "X1", ylab = "X2",
        main = sprintf("UNCERTAINTY: 95%% band width  (Cell D)\nmean width = %.3f, max = %.3f",
                       mean(surfD$band_width), max(surfD$band_width)))
  contour(x1, x2, surfD$band_width, add = TRUE, lwd = 1.2, labcex = 0.9)

  # Page 4: Cell B posterior mean (NULL DGP shrinkage check)
  if (!is.null(surfB)) {
    image(x1, x2, surfB$post_mean, col = pal, breaks = z_breaks,
          xlab = "X1", ylab = "X2",
          main = sprintf("POSTERIOR MEAN f_12  (Cell B: NULL + M1, truth = 0)\nRMSE = %.4f, coverage of 0 = %.2f",
                         results$null_M1$rmse_f12, results$null_M1$cov_f12))
    contour(x1, x2, surfB$post_mean, add = TRUE, lwd = 1.2, labcex = 0.9)
  }

  # Page 5: Cell D residual surface
  resid_surf <- surfD$post_mean - truth_int
  r_lim <- max(abs(resid_surf))
  r_breaks <- seq(-r_lim, r_lim, length.out = 51)
  image(x1, x2, resid_surf, col = pal, breaks = r_breaks,
        xlab = "X1", ylab = "X2",
        main = sprintf("RESIDUAL: posterior mean - truth  (Cell D)\nrange [%.2f, %.2f]",
                       min(resid_surf), max(resid_surf)))
  contour(x1, x2, resid_surf, add = TRUE, lwd = 1.2, labcex = 0.9)

  dev.off()
  cat(sprintf("  Interaction surface plots written to: %s\n", outfile))
}


# ============================================================
# plot_2x2_diagnostics(): recovery and trace plots, all to one PDF
# ============================================================
plot_2x2_diagnostics <- function(res_obj, p, outfile) {
  pdf(outfile, width = 11, height = 8.5)

  results <- res_obj$results
  sim_null <- res_obj$sim_null
  sim_int  <- res_obj$sim_int

  # ---- Page 1: variance component traces, INT DGP, M1 cell (Cell D) ----
  gs_D <- attr(results$int_M1, "gs")  # may be NULL if not stored; fall back
  # (We didn't store gs in summarise_cell; recover from saved obj if present)

  # Instead use the stored summaries: just draw bar chart of WAIC
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  waic_vec <- c(A = results$null_M0$waic, B = results$null_M1$waic,
                C = results$int_M0$waic,  D = results$int_M1$waic)
  barplot(waic_vec, main = sprintf("WAIC by cell (p=%d)", p),
          ylab = "WAIC (lower better)",
          col = c("steelblue", "steelblue", "tomato", "tomato"),
          names.arg = c("A:N+M0", "B:N+M1", "C:I+M0", "D:I+M1"))

  # RMSE_y by cell
  rmse_y_vec <- c(A = results$null_M0$rmse_y, B = results$null_M1$rmse_y,
                  C = results$int_M0$rmse_y,  D = results$int_M1$rmse_y)
  barplot(rmse_y_vec, main = sprintf("In-sample RMSE_y (p=%d)", p),
          ylab = "RMSE",
          col = c("steelblue", "steelblue", "tomato", "tomato"),
          names.arg = c("A:N+M0", "B:N+M1", "C:I+M0", "D:I+M1"))

  # sigma2 estimates vs truth
  sig_vec <- c(A = results$null_M0$sigma2_mean, B = results$null_M1$sigma2_mean,
               C = results$int_M0$sigma2_mean,  D = results$int_M1$sigma2_mean)
  barplot(sig_vec, main = sprintf("Posterior mean sigma2 (truth=%.2f)",
                                   sim_null$true_sigma2),
          ylab = "sigma2",
          col = c("steelblue", "steelblue", "tomato", "tomato"),
          names.arg = c("A:N+M0", "B:N+M1", "C:I+M0", "D:I+M1"))
  abline(h = sim_null$true_sigma2, lty = 2, col = "darkgreen", lwd = 2)

  # Elapsed time
  t_vec <- c(A = results$null_M0$elapsed, B = results$null_M1$elapsed,
             C = results$int_M0$elapsed,  D = results$int_M1$elapsed)
  barplot(t_vec, main = sprintf("Elapsed time (sec), p=%d", p),
          ylab = "seconds",
          col = c("steelblue", "steelblue", "tomato", "tomato"),
          names.arg = c("A:N+M0", "B:N+M1", "C:I+M0", "D:I+M1"))

  # ---- Page 2: f_12 recovery -- INT DGP both models ----
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # Cell D: M1 fitted f12 vs true f12 (scatter)
  if (!is.null(results$int_M1$f12_hat)) {
    plot(sim_int$true_f12, results$int_M1$f12_hat,
         pch = 20, col = "tomato",
         main = sprintf("INT DGP: M1 f12 recovery\nRMSE=%.4f",
                        results$int_M1$rmse_f12),
         xlab = "true f12", ylab = "fitted f12 (post mean)")
    abline(0, 1, lty = 2, col = "darkgreen")
  }

  # Cell B: M1 fitted f12 vs zero (NULL DGP)
  if (!is.null(results$null_M1$f12_hat)) {
    plot(sim_null$true_f12, results$null_M1$f12_hat,
         pch = 20, col = "steelblue",
         main = sprintf("NULL DGP: M1 f12 (truth=0)\nRMSE=%.4f cov=%.2f",
                        results$null_M1$rmse_f12, results$null_M1$cov_f12),
         xlab = "true f12 (=0)", ylab = "fitted f12 (post mean)",
         ylim = range(c(results$null_M1$f12_lo, results$null_M1$f12_hi)))
    abline(h = 0, lty = 2, col = "darkgreen")
  }

  # Cell D: f12 pointwise CI bands, sorted
  if (!is.null(results$int_M1$f12_hat)) {
    ord <- order(sim_int$true_f12)
    plot(seq_along(ord), sim_int$true_f12[ord], type = "l", lwd = 2,
         ylim = range(c(results$int_M1$f12_lo, results$int_M1$f12_hi,
                        sim_int$true_f12)),
         main = "INT DGP: M1 f12 with 95% CI",
         xlab = "obs (sorted by truth)", ylab = "f12")
    lines(seq_along(ord), results$int_M1$f12_hat[ord], col = "tomato", lwd = 2)
    lines(seq_along(ord), results$int_M1$f12_lo[ord],  col = "tomato", lty = 2)
    lines(seq_along(ord), results$int_M1$f12_hi[ord],  col = "tomato", lty = 2)
    legend("topleft", c("true", "fitted", "95% CI"),
           col = c("black", "tomato", "tomato"),
           lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n")
  }

  # Cell B: f12 pointwise CI bands (should bracket 0)
  if (!is.null(results$null_M1$f12_hat)) {
    ord <- order(results$null_M1$f12_hat)
    plot(seq_along(ord), results$null_M1$f12_hat[ord], type = "l",
         col = "steelblue", lwd = 2,
         ylim = range(c(results$null_M1$f12_lo, results$null_M1$f12_hi)),
         main = sprintf("NULL DGP: M1 f12 bands (cov of 0 = %.2f)",
                        results$null_M1$cov_f12),
         xlab = "obs (sorted by post mean)", ylab = "f12")
    lines(seq_along(ord), results$null_M1$f12_lo[ord], lty = 2, col = "steelblue")
    lines(seq_along(ord), results$null_M1$f12_hi[ord], lty = 2, col = "steelblue")
    abline(h = 0, col = "darkgreen", lwd = 2)
  }

  dev.off()
  cat(sprintf("  Plots written to: %s\n", outfile))
}


# ============================================================
# print_2x2_table(): pretty 2x2 summary
# ============================================================
print_2x2_table <- function(results, p) {
  cells <- list(
    A = results$null_M0,    # NULL DGP, M0
    B = results$null_M1,    # NULL DGP, M1
    C = results$int_M0,     # INT DGP, M0
    D = results$int_M1      # INT DGP, M1
  )

  cat("\n========================================================\n")
  cat(sprintf(" 2x2 RESULTS: p = %d\n", p))
  cat("========================================================\n\n")

  cat(sprintf("%-12s | %-22s | %-22s\n", "", "M0 (no interaction)", "M1 (with interaction)"))
  cat(strrep("-", 64), "\n", sep = "")
  for (dgp in c("NULL", "INT")) {
    cell_M0 <- if (dgp == "NULL") cells$A else cells$C
    cell_M1 <- if (dgp == "NULL") cells$B else cells$D
    rmse_main_str_M0 <- paste(sprintf("%.3f", cell_M0$rmse_main), collapse = ",")
    rmse_main_str_M1 <- paste(sprintf("%.3f", cell_M1$rmse_main), collapse = ",")
    cat(sprintf("%-12s | RMSE_main=%-12s | RMSE_main=%-12s\n",
                paste0("DGP-", dgp), rmse_main_str_M0, rmse_main_str_M1))
    cat(sprintf("%-12s | RMSE_y   =%-12.4f | RMSE_y   =%-12.4f\n",
                "", cell_M0$rmse_y, cell_M1$rmse_y))
    cat(sprintf("%-12s | sigma2   =%-12.3f | sigma2   =%-12.3f\n",
                "", cell_M0$sigma2_mean, cell_M1$sigma2_mean))
    cat(sprintf("%-12s | tau2     =%-12.3f | tau2     =%-12.3f\n",
                "", cell_M0$tau2_mean, cell_M1$tau2_mean))
    cat(sprintf("%-12s | rho      =%-12.3f | rho      =%-12.3f\n",
                "", cell_M0$rho_mean, cell_M1$rho_mean))
    cat(sprintf("%-12s | WAIC     =%-12.1f | WAIC     =%-12.1f\n",
                "", cell_M0$waic, cell_M1$waic))
    if (!is.na(cell_M1$rmse_f12)) {
      cat(sprintf("%-12s | f12      = (n/a)        | RMSE_f12 =%-12.4f\n",
                  "", cell_M1$rmse_f12))
      cat(sprintf("%-12s |                         | cov_f12  =%-12.3f\n",
                  "", cell_M1$cov_f12))
    }
    cat(sprintf("%-12s | time(s)  =%-12.1f | time(s)  =%-12.1f\n",
                "", cell_M0$elapsed, cell_M1$elapsed))
    cat(strrep("-", 64), "\n", sep = "")
  }

  # Headline comparisons
  cat("\nKey comparisons:\n")
  cat(sprintf("  WAIC delta (NULL: M1-M0) = %+.1f  [>0 means M0 preferred -> good]\n",
              cells$B$waic - cells$A$waic))
  cat(sprintf("  WAIC delta (INT:  M1-M0) = %+.1f  [<0 means M1 preferred -> good]\n",
              cells$D$waic - cells$C$waic))
  cat(sprintf("  Time multiplier M1/M0 (INT DGP) = %.2fx\n",
              cells$D$elapsed / cells$C$elapsed))
}


# ============================================================
# run_2x2_p2(): full 2x2 experiment for p=2
# ============================================================
run_2x2_p2 <- function(n = 300, M = 8, n_iter = 3000, n_burn = 1000,
                        c_int = 1.5,
                        data_seed_null = 42, data_seed_int = 43,
                        mcmc_seed = 100, verbose = FALSE) {
  cat("\n##########  p = 2 EXPERIMENT  ##########\n")
  cat(sprintf("n=%d  M=%d  n_iter=%d  n_burn=%d  c_int=%.2f\n",
              n, M, n_iter, n_burn, c_int))

  # Two DGPs
  sim_null <- simulate_p2(n, c_int = 0,     data_seed = data_seed_null)
  sim_int  <- simulate_p2(n, c_int = c_int, data_seed = data_seed_int)

  # Designs (both DGPs share design structure; only y differs)
  des_M0 <- build_design_p2(sim_null$X1, sim_null$X2, M, include_int = FALSE)
  des_M1 <- build_design_p2(sim_null$X1, sim_null$X2, M, include_int = TRUE)
  des_M0_int <- build_design_p2(sim_int$X1, sim_int$X2, M, include_int = FALSE)
  des_M1_int <- build_design_p2(sim_int$X1, sim_int$X2, M, include_int = TRUE)

  cat("\nFitting Cell A (NULL + M0)...\n")
  fit_A <- fit_one_cell(sim_null$y, des_M0$H, sim_null$D,
                         des_M0$col_map_main, des_M0$K_main_list,
                         des_M0$col_map_int, des_M0$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_A$elapsed))

  cat("Fitting Cell B (NULL + M1)...\n")
  fit_B <- fit_one_cell(sim_null$y, des_M1$H, sim_null$D,
                         des_M1$col_map_main, des_M1$K_main_list,
                         des_M1$col_map_int, des_M1$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_B$elapsed))

  cat("Fitting Cell C (INT + M0)...\n")
  fit_C <- fit_one_cell(sim_int$y, des_M0_int$H, sim_int$D,
                         des_M0_int$col_map_main, des_M0_int$K_main_list,
                         des_M0_int$col_map_int, des_M0_int$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_C$elapsed))

  cat("Fitting Cell D (INT + M1)...\n")
  fit_D <- fit_one_cell(sim_int$y, des_M1_int$H, sim_int$D,
                         des_M1_int$col_map_main, des_M1_int$K_main_list,
                         des_M1_int$col_map_int, des_M1_int$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_D$elapsed))

  # ---- Grids for marginal/surface plotting ----
  g1d <- seq(0, 1, length.out = 200)
  grid_main <- list(g1d, g1d)               # f1, f2 grids
  g2d <- seq(0, 1, length.out = 30)         # 30x30 surface grid
  gg2d <- expand.grid(x1 = g2d, x2 = g2d)
  grid_int <- list(x1 = g2d, x2 = g2d,
                   X1_long = gg2d$x1, X2_long = gg2d$x2)

  # Truth on grids (for plotting)
  truth_grid <- list(
    f1  = sin(2 * pi * g1d),
    f2  = (2 * g1d - 1)^2,
    f12 = matrix(c_int * sin(pi * gg2d$x1) * cos(pi * gg2d$x2),
                 length(g2d), length(g2d))
  )
  # NULL DGP truths (for clarity in plots)
  truth_grid_null <- list(
    f1  = sin(2 * pi * g1d),
    f2  = (2 * g1d - 1)^2,
    f12 = matrix(0, length(g2d), length(g2d))
  )

  results <- list(
    null_M0 = summarise_cell(fit_A, des_M0,     sim_null, p = 2,
                              grid_main = grid_main, grid_int = NULL),
    null_M1 = summarise_cell(fit_B, des_M1,     sim_null, p = 2,
                              grid_main = grid_main, grid_int = grid_int),
    int_M0  = summarise_cell(fit_C, des_M0_int, sim_int,  p = 2,
                              grid_main = grid_main, grid_int = NULL),
    int_M1  = summarise_cell(fit_D, des_M1_int, sim_int,  p = 2,
                              grid_main = grid_main, grid_int = grid_int)
  )

  print_2x2_table(results, p = 2)

  invisible(list(results = results,
                 sim_null = sim_null, sim_int = sim_int,
                 settings = list(n = n, M = M, n_iter = n_iter, n_burn = n_burn,
                                 c_int = c_int),
                 truth_grid_null = truth_grid_null,
                 truth_grid_int  = truth_grid))
}


# ============================================================
# run_2x2_p3(): full 2x2 experiment for p=3
# ============================================================
run_2x2_p3 <- function(n = 300, M = 8, n_iter = 3000, n_burn = 1000,
                        c_int = 1.5,
                        data_seed_null = 142, data_seed_int = 143,
                        mcmc_seed = 200, verbose = FALSE) {
  cat("\n##########  p = 3 EXPERIMENT  ##########\n")
  cat(sprintf("n=%d  M=%d  n_iter=%d  n_burn=%d  c_int=%.2f\n",
              n, M, n_iter, n_burn, c_int))
  cat("Note: M1 has 3 main + 3 pairwise interactions; only X1xX2 is true-nonzero in INT DGP\n")

  sim_null <- simulate_p3(n, c_int = 0,     data_seed = data_seed_null)
  sim_int  <- simulate_p3(n, c_int = c_int, data_seed = data_seed_int)

  des_M0     <- build_design_p3(sim_null$X1, sim_null$X2, sim_null$X3, M, include_int = FALSE)
  des_M1     <- build_design_p3(sim_null$X1, sim_null$X2, sim_null$X3, M, include_int = TRUE)
  des_M0_int <- build_design_p3(sim_int$X1,  sim_int$X2,  sim_int$X3,  M, include_int = FALSE)
  des_M1_int <- build_design_p3(sim_int$X1,  sim_int$X2,  sim_int$X3,  M, include_int = TRUE)

  cat(sprintf("\n  H dims: M0 = %d x %d  | M1 = %d x %d\n",
              nrow(des_M0$H), ncol(des_M0$H),
              nrow(des_M1$H), ncol(des_M1$H)))

  cat("\nFitting Cell A (NULL + M0)...\n")
  fit_A <- fit_one_cell(sim_null$y, des_M0$H, sim_null$D,
                         des_M0$col_map_main, des_M0$K_main_list,
                         des_M0$col_map_int, des_M0$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_A$elapsed))

  cat("Fitting Cell B (NULL + M1)...\n")
  fit_B <- fit_one_cell(sim_null$y, des_M1$H, sim_null$D,
                         des_M1$col_map_main, des_M1$K_main_list,
                         des_M1$col_map_int, des_M1$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_B$elapsed))

  cat("Fitting Cell C (INT + M0)...\n")
  fit_C <- fit_one_cell(sim_int$y, des_M0_int$H, sim_int$D,
                         des_M0_int$col_map_main, des_M0_int$K_main_list,
                         des_M0_int$col_map_int, des_M0_int$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_C$elapsed))

  cat("Fitting Cell D (INT + M1)...\n")
  fit_D <- fit_one_cell(sim_int$y, des_M1_int$H, sim_int$D,
                         des_M1_int$col_map_main, des_M1_int$K_main_list,
                         des_M1_int$col_map_int, des_M1_int$K_int_list,
                         n_iter, n_burn, mcmc_seed, verbose = verbose)
  cat(sprintf("  done in %.1f sec\n", fit_D$elapsed))

  # ---- Grids for marginal/surface plotting ----
  g1d <- seq(0, 1, length.out = 200)
  grid_main <- list(g1d, g1d, g1d)         # f1, f2, f3 grids
  g2d <- seq(0, 1, length.out = 30)
  gg2d <- expand.grid(x1 = g2d, x2 = g2d)
  grid_int <- list(x1 = g2d, x2 = g2d,
                   X1_long = gg2d$x1, X2_long = gg2d$x2)

  truth_grid <- list(
    f1  = sin(2 * pi * g1d),
    f2  = (2 * g1d - 1)^2,
    f3  = 0.8 * (g1d - 0.5),
    f12 = matrix(c_int * sin(pi * gg2d$x1) * cos(pi * gg2d$x2),
                 length(g2d), length(g2d))
  )
  truth_grid_null <- list(
    f1  = sin(2 * pi * g1d),
    f2  = (2 * g1d - 1)^2,
    f3  = 0.8 * (g1d - 0.5),
    f12 = matrix(0, length(g2d), length(g2d))
  )

  results <- list(
    null_M0 = summarise_cell(fit_A, des_M0,     sim_null, p = 3,
                              grid_main = grid_main, grid_int = NULL),
    null_M1 = summarise_cell(fit_B, des_M1,     sim_null, p = 3,
                              grid_main = grid_main, grid_int = grid_int),
    int_M0  = summarise_cell(fit_C, des_M0_int, sim_int,  p = 3,
                              grid_main = grid_main, grid_int = NULL),
    int_M1  = summarise_cell(fit_D, des_M1_int, sim_int,  p = 3,
                              grid_main = grid_main, grid_int = grid_int)
  )

  print_2x2_table(results, p = 3)

  invisible(list(results = results,
                 sim_null = sim_null, sim_int = sim_int,
                 settings = list(n = n, M = M, n_iter = n_iter, n_burn = n_burn,
                                 c_int = c_int),
                 truth_grid_null = truth_grid_null,
                 truth_grid_int  = truth_grid))
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
  cat(" 2x2 INTERACTION COMPARISON\n")
  cat(" {NULL DGP, INT DGP} x {Model M0 (no int), Model M1 (with int)}\n")
  cat("==============================================================\n")

  t_total <- proc.time()

  res_p2 <- run_2x2_p2(n = 300, M = 8, n_iter = 3000, n_burn = 1000,
                        c_int = 1.5, mcmc_seed = 100)

  res_p3 <- run_2x2_p3(n = 300, M = 8, n_iter = 3000, n_burn = 1000,
                        c_int = 1.5, mcmc_seed = 200)

  total_elapsed <- as.numeric((proc.time() - t_total)["elapsed"])
  cat(sprintf("\n##########  TOTAL ELAPSED: %.1f sec (%.1f min)  ##########\n",
              total_elapsed, total_elapsed / 60))

  # Save: RDS (full), CSV (tidy summary), PDFs (diagnostics)
  saveRDS(list(p2 = res_p2, p3 = res_p3, total_sec = total_elapsed),
          file = "interaction_2x2_results.rds")

  df_p2 <- results_to_dataframe(res_p2$results, p = 2, settings = res_p2$settings)
  df_p3 <- results_to_dataframe(res_p3$results, p = 3, settings = res_p3$settings)
  df_all <- rbind(df_p2, df_p3)
  write.csv(df_all, file = "interaction_2x2_summary.csv", row.names = FALSE)
  cat("CSV summary saved to interaction_2x2_summary.csv\n")

  plot_2x2_diagnostics(res_p2, p = 2, outfile = "interaction_2x2_p2.pdf")
  plot_2x2_diagnostics(res_p3, p = 3, outfile = "interaction_2x2_p3.pdf")

  plot_marginals_2x2(res_p2, p = 2, outfile = "interaction_2x2_p2_marginals.pdf")
  plot_marginals_2x2(res_p3, p = 3, outfile = "interaction_2x2_p3_marginals.pdf")

  plot_interaction_surfaces(res_p2, p = 2, outfile = "interaction_2x2_p2_surfaces.pdf")
  plot_interaction_surfaces(res_p3, p = 3, outfile = "interaction_2x2_p3_surfaces.pdf")

  cat("Results saved to interaction_2x2_results.rds\n")
}
