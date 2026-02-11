# ============================================================
# marginal_utils.R
# For LS-basis spatial additive model fitted by fit_ls_spatial()
# Outputs:
#   (1) importance table on observed data points
#   (2) curve error table on x-grid in [0,1]
#   (3) marginal curve plots: f_hat_j(x) vs f_true_j(x)
# ============================================================

rmse <- function(a, b) sqrt(mean((a - b)^2))

# ------------------------------------------------------------
# Extract per-covariate coefficient blocks beta_j from fitted beta
# obj: output from fit_ls_spatial()
# Returns list of beta blocks (length p), each numeric vector
# ------------------------------------------------------------
extract_beta_blocks <- function(obj) {
  beta_all <- as.numeric(obj$fit$beta)
  col_map <- obj$des$col_map
  p <- length(col_map)
  beta_blocks <- vector("list", p)
  for (j in 1:p) {
    # +1 because beta_all[1] is intercept (mu)
    beta_blocks[[j]] <- beta_all[1 + col_map[[j]]]
  }
  beta_blocks
}

# ------------------------------------------------------------
# Compute f_hat_j(X_j(s_r)) at observed data points (r=1..n)
# Returns list of length p, each is n-vector of contributions
# ------------------------------------------------------------
fhat_at_data <- function(obj) {
  W <- obj$des$W
  col_map <- obj$des$col_map
  beta_blocks <- extract_beta_blocks(obj)
  p <- length(col_map)
  
  out <- vector("list", p)
  for (j in 1:p) {
    Wj <- W[, col_map[[j]], drop = FALSE]
    out[[j]] <- as.numeric(Wj %*% beta_blocks[[j]])
  }
  out
}

# ------------------------------------------------------------
# Importance table (sd and mean abs) based on f_hat at data points
# ------------------------------------------------------------
importance_table <- function(obj, var_names = NULL) {
  fh_list <- fhat_at_data(obj)
  p <- length(fh_list)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  data.frame(
    var = var_names,
    sd_fhat = sapply(fh_list, sd),
    mean_abs_fhat = sapply(fh_list, function(v) mean(abs(v))),
    range_fhat = sapply(fh_list, function(v) max(v) - min(v))
  )
}

# ------------------------------------------------------------
# Build marginal curves on a grid x in [0,1]
# truth_f_list: optional list of true functions (length p), each f(x) vectorized
# Returns a list with:
#   grid: x_grid
#   fhat_grid: list length p, each length(grid)
#   ftrue_grid: list length p or NULL
# ------------------------------------------------------------
marginal_curves <- function(obj, x_grid = seq(0, 1, length.out = 101),
                            truth_f_list = NULL, clip = TRUE,
                            center_truth = TRUE) {
  des <- obj$des
  beta_blocks <- extract_beta_blocks(obj)
  p <- length(des$objs)
  
  if (!is.null(truth_f_list)) stopifnot(length(truth_f_list) == p)
  
  fhat_grid <- vector("list", p)
  ftrue_grid <- if (is.null(truth_f_list)) NULL else vector("list", p)
  
  stopifnot(!center_truth || !is.null(obj$X_raw_for_marginal))
  if (center_truth) {
    Xraw <- as.matrix(obj$X_raw_for_marginal)
    stopifnot(ncol(Xraw) == p)
  }
  
  for (j in 1:p) {
    Wj_grid <- des$objs[[j]]$design_new(x_grid, type = "W", clip = clip)
    fhat_grid[[j]] <- as.numeric(Wj_grid %*% beta_blocks[[j]])
    
    if (!is.null(truth_f_list)) {
      ytrue <- as.numeric(truth_f_list[[j]](x_grid))
      
      if (center_truth) {
        shift <- mean(truth_f_list[[j]](Xraw[, j]))
        ytrue <- ytrue - shift
      }
      ftrue_grid[[j]] <- ytrue
    }
  }
  
  list(grid = x_grid, fhat_grid = fhat_grid, ftrue_grid = ftrue_grid)
}



# ------------------------------------------------------------
# Curve error table on grid (simulation-only when truth exists)
# Computes RMSE and correlation between fhat_grid and ftrue_grid
# ------------------------------------------------------------
curve_error_table <- function(curves, var_names = NULL) {
  if (is.null(curves$ftrue_grid)) {
    stop("No truth_f_list provided, cannot compute curve errors.")
  }
  p <- length(curves$fhat_grid)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  rmse_vec <- numeric(p)
  cor_vec  <- rep(NA_real_, p)
  
  for (j in 1:p) {
    yhat  <- curves$fhat_grid[[j]]
    ytrue <- curves$ftrue_grid[[j]]
    rmse_vec[j] <- rmse(yhat, ytrue)
    
    # only compute correlation when both have variability
    if (sd(yhat) > 0 && sd(ytrue) > 0) {
      cor_vec[j] <- cor(yhat, ytrue)
    } else {
      cor_vec[j] <- NA_real_
    }
  }
  
  data.frame(
    var = var_names,
    rmse_curve = rmse_vec,
    cor_curve  = cor_vec
  )
}


# ------------------------------------------------------------
# Plot marginal curves (one plot per covariate)
# If truth exists: overlay true and estimated
# base R plot (no ggplot)
# ------------------------------------------------------------
plot_marginals <- function(curves, var_names = NULL, show_rmse = TRUE,
                           ylim_shared = TRUE, ylim = NULL,
                           symmetric = TRUE, pad = 0.05) {
  x <- curves$grid
  p <- length(curves$fhat_grid)
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  
  # decide common ylim
  if (is.null(ylim) && ylim_shared) {
    all_y <- unlist(curves$fhat_grid)
    if (!is.null(curves$ftrue_grid)) all_y <- c(all_y, unlist(curves$ftrue_grid))
    rng <- range(all_y, finite = TRUE)
    if (symmetric) {
      a <- max(abs(rng))
      rng <- c(-a, a)
    }
    # small padding
    span <- diff(rng); if (span == 0) span <- 1
    rng <- rng + c(-1, 1) * pad * span
    ylim <- rng
  }
  
  for (j in 1:p) {
    yhat <- curves$fhat_grid[[j]]
    
    plot(x, yhat, type = "l",
         main = paste0("Marginal: ", var_names[j]),
         xlab = "x", ylab = "f_j(x)",
         ylim = ylim)
    
    abline(h = 0, lty = 3)
    
    if (!is.null(curves$ftrue_grid)) {
      ytrue <- curves$ftrue_grid[[j]]
      lines(x, ytrue, lty = 2)
      
      if (show_rmse) {
        r <- rmse(yhat, ytrue)
        mtext(sprintf("curve RMSE = %.4f", r), side = 3, line = 0.2, cex = 0.85)
      }
      legend("topleft", legend = c("estimated", "true"),
             lty = c(1, 2), bty = "n")
    }
  }
}


# ------------------------------------------------------------
# Default truth functions for your current Sim1/2/3 setup (p=4)
# f1(x)=2 sin(pi x), f2(x)=1.5 exp(x-0.5), f3(x)=0.7 x^2, f4(x)=0.5 sin(2 pi x)
# Returns a list of 4 functions.
# ------------------------------------------------------------
default_truth_f_list <- function() {
  list(
    function(x) 2 * sin(pi * x),
    function(x) 1.5 * exp(x - 0.5),
    function(x) 0.7 * (x^2),
    function(x) 0.5 * sin(2 * pi * x)
  )
}
