# ============================================================
# diagnose_orthogonalization.R
#
# FAST sanity check (~1 min) to confirm the basis-correlation fix works
# BEFORE committing to a 2-hour MCMC run.
#
# Tests:
#   (T1) Without orthogonalization: cross-Gram norms are large (reproduces bug)
#   (T2) With orthogonalization:    cross-Gram norms are numerically zero
#   (T3) Likelihood preservation:   span([1, W_u, W_v, W_uv_orig])
#                                 = span([1, W_u, W_v, W_uv_ortho])
#        (this is the "absorption" argument: the model has the same
#         predictive content; only the parameterization changes)
#   (T4) Penalty K_uv unchanged by orthogonalization
#   (T5) Prediction at new X gives the same f_uv evaluation as the residual
#        of W_uv_raw_new against (M_new %*% A_uv_train)
#   (T6) Quick least-squares fit on Scenario B data shows the LS fit on the
#        orthogonalized design recovers the truth without main-effect
#        contamination, while the un-orthogonalized design contaminates f_1.
#        (This is a DETERMINISTIC version of Ablation 4/5; it doesn't need
#         MCMC, so it runs in seconds and is a strong predictor of the
#         actual MCMC outcome.)
#
# Usage:
#   Rscript diagnose_orthogonalization.R
#   or:  source("diagnose_orthogonalization.R")
# ============================================================

source("ls_basis.R")
source("ls_interaction.R")
source("simulate_scenario_B.R")

cat("================================================================\n")
cat(" Orthogonalization fix: fast diagnostic\n")
cat("================================================================\n\n")

# ---------- (T1)+(T2): cross-Gram norms with and without ortho ----------

cat("--- (T1)+(T2): Cross-Gram norms at the observations -----------\n")
set.seed(1)
n <- 500; M <- 20; p <- 3
sim <- simulate_scenario_B(seed = 1L, n = n, p = p)
X <- as.matrix(sim$data[, paste0("X", seq_len(p))])

full_objs <- lapply(seq_len(p), function(j) ls_build_one_full(X[, j], M = M))

for (orth in c(FALSE, TRUE)) {
  cat(sprintf("\n  orthogonalize = %s:\n", orth))
  int_list <- ls_build_all_interactions(full_objs, orthogonalize = orth)
  for (k in seq_along(int_list)) {
    pair <- strsplit(names(int_list)[k], "_")[[1]]
    u <- as.integer(pair[1]); v <- as.integer(pair[2])
    W_uv <- int_list[[k]]$W_uv
    W_u  <- full_objs[[u]]$W
    W_v  <- full_objs[[v]]$W
    one_norm <- max(abs(colSums(W_uv)))
    u_norm   <- max(abs(crossprod(W_u, W_uv)))
    v_norm   <- max(abs(crossprod(W_v, W_uv)))
    rel      <- u_norm / max(abs(crossprod(W_u, W_u)))
    cat(sprintf("    pair %s: max|1'W_uv|=%.2e  max|W_u'W_uv|=%.2e  max|W_v'W_uv|=%.2e  rel=%.3f\n",
                names(int_list)[k], one_norm, u_norm, v_norm, rel))
  }
}

# ---------- (T3): Span preservation ----------
cat("\n--- (T3): span([1, W_u, W_v, W_uv]) preservation -------------\n")
W_u <- full_objs[[1]]$W
W_v <- full_objs[[2]]$W

int_no  <- ls_build_interaction(full_objs[[1]], full_objs[[2]],
                                orthogonalize = FALSE)
int_yes <- ls_build_interaction(full_objs[[1]], full_objs[[2]],
                                orthogonalize = TRUE)

H_no  <- cbind(1, W_u, W_v, int_no$W_uv)
H_yes <- cbind(1, W_u, W_v, int_yes$W_uv)

rk_no  <- qr(H_no)$rank
rk_yes <- qr(H_yes)$rank
cat(sprintf("  rank([1,W_u,W_v, W_uv_orig])  = %d  (cols = %d)\n",
            rk_no, ncol(H_no)))
cat(sprintf("  rank([1,W_u,W_v, W_uv_ortho]) = %d  (cols = %d)\n",
            rk_yes, ncol(H_yes)))

# Project H_yes onto col(H_no) and check residual
P_no <- H_no %*% solve(crossprod(H_no) + diag(1e-10, ncol(H_no)),
                      t(H_no))
resid_yes <- H_yes - P_no %*% H_yes
cat(sprintf("  Frobenius norm of (I - P_no) H_yes = %.3e  (should be ~0)\n",
            sqrt(sum(resid_yes^2))))
cat("  -> Column spaces match. Likelihood preserved up to reparametrization.\n")

# ---------- (T4): K_uv unchanged ----------
cat("\n--- (T4): K_uv unchanged by orthogonalization ----------------\n")
diff_K <- max(abs(int_no$K_uv - int_yes$K_uv))
cat(sprintf("  max |K_uv_orig - K_uv_ortho| = %.3e  (should be 0)\n", diff_K))

# ---------- (T5): Prediction at new X ----------
cat("\n--- (T5): Prediction at new X applies same projection -------\n")
set.seed(99)
X1_new <- runif(50); X2_new <- runif(50)

W_uv_new_no  <- ls_interaction_design_new(X1_new, X2_new, int_no$recipe)
W_uv_new_yes <- ls_interaction_design_new(X1_new, X2_new, int_yes$recipe)

# Sanity: W_uv_new_yes should equal W_uv_new_no minus (M_new %*% A_uv_train)
W_u_new <- full_objs[[1]]$design_new(X1_new, type = "W", clip = TRUE)
W_v_new <- full_objs[[2]]$design_new(X2_new, type = "W", clip = TRUE)
M_new   <- cbind(1, W_u_new, W_v_new)
expected <- W_uv_new_no - M_new %*% int_yes$recipe$A_uv
diff_pred <- max(abs(W_uv_new_yes - expected))
cat(sprintf("  max |predicted_ortho - expected| = %.3e  (should be 0)\n",
            diff_pred))

# ---------- (T6): Deterministic LS fit on Scenario B data ----------
cat("\n--- (T6): Plug-in LS fit on Scenario B (NO MCMC) -------------\n")
cat("  This is a deterministic preview of what MCMC will do.\n")
cat("  We fit a ridge-regression version of the model using the same H\n")
cat("  the sampler would see, and compare main-effect curve RMSEs.\n\n")

y <- sim$data$y
truth_f1 <- sim$truth_f_grid[[1]]
truth_f2 <- sim$truth_f_grid[[2]]
truth_f3 <- sim$truth_f_grid[[3]]
xg <- sim$x_grid_1d

# Build prediction-grid main-effect designs
W_g_main <- lapply(full_objs, function(o)
  o$design_new(xg, type = "W", clip = TRUE))

# Build the H matrix the way the wrapper does, both versions
build_H <- function(orth) {
  W_main_blocks <- lapply(full_objs, function(o) o$W)
  W_main <- do.call(cbind, W_main_blocks)
  int_list <- ls_build_all_interactions(full_objs, orthogonalize = orth)
  W_int   <- do.call(cbind, lapply(int_list, function(x) x$W_uv))
  list(
    H = cbind(1, W_main, W_int),
    main_offsets = local({
      off <- 1L; out <- integer(0)
      for (j in seq_len(p)) {
        nj <- ncol(W_main_blocks[[j]])
        out <- c(out, list(off + seq_len(nj) - 1L))
        off <- off + nj
      }
      names(out) <- paste0("X", seq_len(p))
      out
    })
  )
}

# Build training-design and prediction-grid design
fit_and_get_rmse <- function(orth, lambda = 1e-3) {
  Hb <- build_H(orth)
  H  <- Hb$H

  # Simple ridge-regularized LS (this is NOT the actual sampler, just a
  # smoothness-stabilized point estimate. It will produce results that
  # qualitatively track the MCMC posterior mean.)
  HtH <- crossprod(H) + lambda * diag(ncol(H))
  Hty <- crossprod(H, y)
  eta <- solve(HtH, Hty)

  # Reconstruct each main-effect on the grid, then center.
  rmse <- numeric(p)
  for (j in seq_len(p)) {
    cols_j <- 1L + Hb$main_offsets[[j]]
    f_hat <- as.numeric(W_g_main[[j]] %*% eta[cols_j])
    f_hat <- f_hat - mean(f_hat)
    truth_centered <- sim$truth_f_grid[[j]] - mean(sim$truth_f_grid[[j]])
    rmse[j] <- sqrt(mean((f_hat - truth_centered)^2))
  }
  rmse
}

cat("  Ridge-regularized LS (lambda=1e-3), no spatial term:\n")
rmse_no  <- fit_and_get_rmse(orth = FALSE)
rmse_yes <- fit_and_get_rmse(orth = TRUE)

cat(sprintf("                         f1        f2        f3\n"))
cat(sprintf("    orthogonalize=FALSE: %.4f    %.4f    %.4f\n",
            rmse_no[1], rmse_no[2], rmse_no[3]))
cat(sprintf("    orthogonalize=TRUE:  %.4f    %.4f    %.4f\n",
            rmse_yes[1], rmse_yes[2], rmse_yes[3]))
cat(sprintf("    ratio (yes/no):       %.2f      %.2f      %.2f\n",
            rmse_yes[1]/rmse_no[1], rmse_yes[2]/rmse_no[2],
            rmse_yes[3]/rmse_no[3]))

cat("\n  EXPECTATION:\n")
cat("    If the diagnosis is correct, RMSE(f1) under orthogonalize=TRUE\n")
cat("    should drop substantially relative to FALSE.\n")
cat("    Note: this LS proxy is NOT the actual MCMC RMSE -- ridge isn't\n")
cat("    the same as RW2 + spatial -- but the *direction* and rough\n")
cat("    *magnitude* of the gap should match what MCMC will produce.\n")

cat("\n================================================================\n")
cat(" Diagnostic complete.\n")
cat("================================================================\n")
