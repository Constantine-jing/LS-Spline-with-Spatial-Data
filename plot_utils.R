# ============================================================
# plot_utils.R
#   Plot helpers for marginal curves
#   - plot_marginals_fixed_ylim(): same y-scale across panels
#   - plot_marginals(): per-panel auto y-scale (fallback)
#
#   Robust to curves missing x_grid:
#     - you can pass x_grid explicitly, OR
#     - it will try curves$x_grid / curves$x, OR
#     - it will infer from length(fhat) and use seq(0,1, ...)
# ============================================================

.get_curves_field <- function(curves, name_candidates) {
  for (nm in name_candidates) {
    if (!is.null(curves[[nm]])) return(curves[[nm]])
  }
  NULL
}

# ---- robust x-grid getter ----
get_x_grid <- function(curves, x_grid = NULL, default_domain = c(0,1)) {
  if (!is.null(x_grid)) return(x_grid)
  
  xg <- .get_curves_field(curves, c("x_grid", "x", "grid", "xg"))
  if (!is.null(xg)) return(xg)
  
  # infer from estimated curve length
  est <- .get_curves_field(curves, c("estimated", "fhat_grid", "fhat"))
  if (is.null(est)) stop("curves object has no 'estimated' (or fhat_grid) field.")
  L <- length(est[[1]])
  seq(default_domain[1], default_domain[2], length.out = L)
}

compute_global_ylim <- function(curves, include_true = TRUE) {
  est <- .get_curves_field(curves, c("estimated", "fhat_grid", "fhat"))
  tru <- if (include_true) .get_curves_field(curves, c("true", "ftrue_grid", "ftrue")) else NULL
  
  if (is.null(est)) stop("curves object has no 'estimated' (or fhat_grid) field.")
  p <- length(est)
  
  ys <- c()
  for (j in 1:p) {
    ys <- c(ys, est[[j]])
    if (!is.null(tru)) ys <- c(ys, tru[[j]])
  }
  range(ys, na.rm = TRUE)
}

# ------------------------------------------------------------
# Fixed y-scale across panels
# ------------------------------------------------------------
plot_marginals_fixed_ylim <- function(curves,
                                      var_names = NULL,
                                      show_rmse = TRUE,
                                      ylim = NULL,
                                      include_true = TRUE,
                                      x_grid = NULL,
                                      default_domain = c(0,1),
                                      mfrow = NULL,
                                      add_zero_line = TRUE) {
  
  est <- .get_curves_field(curves, c("estimated", "fhat_grid", "fhat"))
  tru <- if (include_true) .get_curves_field(curves, c("true", "ftrue_grid", "ftrue")) else NULL
  rmse_vec <- .get_curves_field(curves, c("rmse", "rmse_curve"))
  
  if (is.null(est)) stop("curves object has no 'estimated' (or fhat_grid) field.")
  p <- length(est)
  
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  stopifnot(length(var_names) == p)
  
  xg <- get_x_grid(curves, x_grid = x_grid, default_domain = default_domain)
  
  if (is.null(ylim)) ylim <- compute_global_ylim(curves, include_true = include_true)
  
  if (is.null(mfrow)) {
    if (p <= 4) mfrow <- c(2, 2) else mfrow <- c(2, 3)
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = mfrow, mar = mar)
  
  for (j in 1:p) {
    fhat <- est[[j]]
    ftru <- if (is.null(tru)) NULL else tru[[j]]
    
    main_txt <- paste0("Marginal: ", var_names[j])
    if (show_rmse && !is.null(rmse_vec) && length(rmse_vec) >= j && !is.na(rmse_vec[j])) {
      main_txt <- paste0(main_txt, "\ncurve RMSE = ", sprintf("%.4f", rmse_vec[j]))
    }
    
    plot(xg, fhat, type = "l", lwd = 2,
         xlab = "x", ylab = "f_j(x)",
         main = main_txt,
         ylim = ylim)
    
    if (add_zero_line) abline(h = 0, lty = 3)
    
    if (!is.null(ftru)) {
      lines(xg, ftru, lty = 2, lwd = 2)
      legend("topleft", legend = c("estimated", "true"),
             lty = c(1, 2), bty = "n")
    } else {
      legend("topleft", legend = c("estimated"),
             lty = 1, bty = "n")
    }
  }
  
  invisible(list(ylim = ylim, mfrow = mfrow))
}

# ------------------------------------------------------------
# Per-panel auto y-scale (fallback)
# ------------------------------------------------------------
plot_marginals <- function(curves,
                           var_names = NULL,
                           show_rmse = TRUE,
                           include_true = TRUE,
                           x_grid = NULL,
                           default_domain = c(0,1),
                           mfrow = NULL,
                           add_zero_line = TRUE) {
  
  est <- .get_curves_field(curves, c("estimated", "fhat_grid", "fhat"))
  tru <- if (include_true) .get_curves_field(curves, c("true", "ftrue_grid", "ftrue")) else NULL
  rmse_vec <- .get_curves_field(curves, c("rmse", "rmse_curve"))
  
  if (is.null(est)) stop("curves object has no 'estimated' (or fhat_grid) field.")
  p <- length(est)
  
  if (is.null(var_names)) var_names <- paste0("X", 1:p)
  stopifnot(length(var_names) == p)
  
  xg <- get_x_grid(curves, x_grid = x_grid, default_domain = default_domain)
  
  if (is.null(mfrow)) {
    if (p <= 4) mfrow <- c(2, 2) else mfrow <- c(2, 3)
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = mfrow, mar = mar)
  
  for (j in 1:p) {
    fhat <- est[[j]]
    ftru <- if (is.null(tru)) NULL else tru[[j]]
    
    ys <- c(fhat, if (!is.null(ftru)) ftru else NULL)
    ylim <- range(ys, na.rm = TRUE)
    
    main_txt <- paste0("Marginal: ", var_names[j])
    if (show_rmse && !is.null(rmse_vec) && length(rmse_vec) >= j && !is.na(rmse_vec[j])) {
      main_txt <- paste0(main_txt, "\ncurve RMSE = ", sprintf("%.4f", rmse_vec[j]))
    }
    
    plot(xg, fhat, type = "l", lwd = 2,
         xlab = "x", ylab = "f_j(x)",
         main = main_txt,
         ylim = ylim)
    
    if (add_zero_line) abline(h = 0, lty = 3)
    
    if (!is.null(ftru)) {
      lines(xg, ftru, lty = 2, lwd = 2)
      legend("topleft", legend = c("estimated", "true"),
             lty = c(1, 2), bty = "n")
    } else {
      legend("topleft", legend = c("estimated"),
             lty = 1, bty = "n")
    }
  }
  
  invisible(list(mfrow = mfrow))
}
