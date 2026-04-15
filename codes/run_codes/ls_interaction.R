# ============================================================
# ls_interaction.R
# Tensor product LS interaction basis for bivariate surface f_{u,v}(x_u, x_v)
#
# Implements Part II–III of ls_tensor_product_construction.pdf
#
# Depends on: ls_basis.R   (ls_build_one_full, ls_contrast_T, etc.)
#
# Main functions:
#   ls_build_interaction()   -- build W_uv and penalty K_uv for one pair (u,v)
#   ls_build_all_interactions() -- build all p*(p-1)/2 pairs
#   ls_interaction_design_new() -- evaluate W_uv at new x values
#   build_rw2_2d_penalty()   -- 2D RW2 Kronecker sum penalty (pure R fallback)
#
# Usage (p=2 proof-of-concept):
#   source("ls_basis.R")
#   source("ls_interaction.R")
#   bu <- ls_build_one_full(X[,1], M = M)
#   bv <- ls_build_one_full(X[,2], M = M)
#   int12 <- ls_build_interaction(bu, bv)
#   # int12$W_uv  : n x (M-1)^2 identified interaction design
#   # int12$K_uv  : (M-1)^2 x (M-1)^2 identified 2D RW2 penalty
# ============================================================


# ============================================================
# build_rw2_penalty_1d(d)
# Standard 1D RW2 penalty, size d x d. Same as in gibbs_stage_c_full.R.
# Duplicated here so ls_interaction.R is self-contained when sourced alone.
# ============================================================
build_rw2_penalty_1d <- function(d) {
  if (d < 3) return(matrix(0, d, d))
  D2 <- matrix(0, d - 2, d)
  for (i in seq_len(d - 2)) {
    D2[i, i]     <-  1
    D2[i, i + 1] <- -2
    D2[i, i + 2] <-  1
  }
  crossprod(D2)   # (d-2 x d) x (d x d-2) = d x d, rank d-2
}


# ============================================================
# build_rw2_2d_penalty(K_u_raw, K_v_raw)
#
# Constructs the 2D Kronecker-sum RW2 penalty on the (M_u * M_v) raw
# coefficient grid BEFORE the T_u x T_v reparametrisation:
#
#   K_uv_raw = K_u_raw ⊗ I_{M_v}  +  I_{M_u} ⊗ K_v_raw
#
# Then transforms to the identified parametrisation (M_u-1)(M_v-1):
#
#   K_uv_beta = (T_u ⊗ T_v)^T  K_uv_raw  (T_u ⊗ T_v)
#
# Inputs:
#   K_u_raw : M_u x M_u  RW2 penalty on raw theta_u  (rank M_u-2)
#   K_v_raw : M_v x M_v  RW2 penalty on raw theta_v  (rank M_v-2)
#   T_u     : M_u x (M_u-1) contrast matrix
#   T_v     : M_v x (M_v-1) contrast matrix
#
# Returns:
#   K_uv_raw  : M_u*M_v x M_u*M_v Kronecker-sum penalty (unidentified)
#   K_uv_beta : (M_u-1)*(M_v-1) x (M_u-1)*(M_v-1) identified penalty
# ============================================================
build_rw2_2d_penalty <- function(K_u_raw, K_v_raw, T_u, T_v) {
  M_u <- nrow(K_u_raw)
  M_v <- nrow(K_v_raw)

  # Kronecker sum:  K_u ⊗ I_v  +  I_u ⊗ K_v
  K_uv_raw <- kronecker(K_u_raw, diag(M_v)) + kronecker(diag(M_u), K_v_raw)

  # Reparametrise: theta_uv = (T_u ⊗ T_v) beta_uv
  T_uv <- kronecker(T_u, T_v)   # M_u*M_v x (M_u-1)*(M_v-1)
  K_uv_beta <- t(T_uv) %*% K_uv_raw %*% T_uv

  list(K_uv_raw = K_uv_raw, K_uv_beta = K_uv_beta, T_uv = T_uv)
}


# ============================================================
# khatri_rao_rowwise(W_u, W_v)
#
# Computes the row-wise Kronecker product (Khatri-Rao product):
#   result[i, ] = W_u[i, ] ⊗ W_v[i, ]    (length (M_u-1)*(M_v-1))
#
# Output: n x (M_u-1)*(M_v-1) matrix = W_uv  (the identified interaction design)
#
# Pure R version — superseded by the Rcpp version in ls_interaction_core.cpp
# when available, but kept as fallback.
# ============================================================
khatri_rao_rowwise_R <- function(W_u, W_v) {
  n   <- nrow(W_u)
  d_u <- ncol(W_u)
  d_v <- ncol(W_v)
  stopifnot(nrow(W_v) == n)

  # Efficient row-wise Kronecker via rep/rep_each tricks
  # W_u[i,] has d_u elements; W_v[i,] has d_v elements
  # result[i, (a-1)*d_v + b] = W_u[i,a] * W_v[i,b]
  W_u_rep <- W_u[, rep(seq_len(d_u), each = d_v), drop = FALSE]   # n x d_u*d_v
  W_v_rep <- W_v[, rep(seq_len(d_v), times = d_u), drop = FALSE]  # n x d_u*d_v
  W_u_rep * W_v_rep   # element-wise product
}


# ============================================================
# ls_build_interaction(obj_u, obj_v, X_u, X_v, use_cpp = TRUE)
#
# Builds everything needed for the interaction surface f_{u,v}(x_u, x_v):
#
#   W_uv  : n x (M_u-1)*(M_v-1)  identified design matrix (Khatri-Rao of W_u, W_v)
#   K_uv  : (M_u-1)*(M_v-1) x (M_u-1)*(M_v-1)  identified 2D RW2 penalty
#   T_uv  : M_u*M_v x (M_u-1)*(M_v-1)  combined contrast (for reference)
#   recipe: stores tau, AinvC, T for prediction on new data
#
# Inputs:
#   obj_u : output of ls_build_one_full(X[,u], M=M_u)  (or ls_build_one_train)
#   obj_v : output of ls_build_one_full(X[,v], M=M_v)
#   X_u, X_v : optional, only needed if obj_u/obj_v don't store W directly.
#              If obj_u$W exists, X_u is ignored.
#   use_cpp  : use Rcpp Khatri-Rao if ls_interaction_core.so is loaded
#
# Note: obj_u and obj_v must have been built from the same X rows (same n).
# ============================================================
ls_build_interaction <- function(obj_u, obj_v, use_cpp = FALSE) {
  # --- Extract W matrices ---
  W_u <- obj_u$W   # n x (M_u - 1)
  W_v <- obj_v$W   # n x (M_v - 1)
  stopifnot(!is.null(W_u), !is.null(W_v))
  stopifnot(nrow(W_u) == nrow(W_v))

  M_u <- obj_u$M
  M_v <- obj_v$M
  T_u <- obj_u$T    # M_u x (M_u-1)
  T_v <- obj_v$T    # M_v x (M_v-1)

  # --- 1. Row-wise Kronecker (Khatri-Rao) product: W_uv ---
  if (use_cpp && exists("khatri_rao_cpp")) {
    W_uv <- khatri_rao_cpp(W_u, W_v)
  } else {
    W_uv <- khatri_rao_rowwise_R(W_u, W_v)
  }

  # --- 2. Raw 1D RW2 penalties on theta_u and theta_v ---
  K_u_raw <- build_rw2_penalty_1d(M_u)   # M_u x M_u
  K_v_raw <- build_rw2_penalty_1d(M_v)   # M_v x M_v

  # --- 3. 2D Kronecker-sum penalty ---
  pen <- build_rw2_2d_penalty(K_u_raw, K_v_raw, T_u, T_v)
  K_uv  <- pen$K_uv_beta   # (M_u-1)*(M_v-1) x (M_u-1)*(M_v-1)
  T_uv  <- pen$T_uv

  # --- 4. Recipe for prediction ---
  recipe <- list(
    M_u = M_u, M_v = M_v,
    T_u = T_u, T_v = T_v,
    AinvC_u = obj_u$AinvC,
    AinvC_v = obj_v$AinvC,
    tau_u = obj_u$tau,
    tau_v = obj_v$tau,
    design_new_u = obj_u$design_new,
    design_new_v = obj_v$design_new
  )

  list(
    W_uv  = W_uv,           # n x (M_u-1)*(M_v-1)   design matrix
    K_uv  = K_uv,           # (M_u-1)*(M_v-1)^2     identified penalty
    T_uv  = T_uv,           # for reference
    K_u_raw = K_u_raw,
    K_v_raw = K_v_raw,
    M_u = M_u, M_v = M_v,
    recipe = recipe
  )
}


# ============================================================
# ls_interaction_design_new(X_u_new, X_v_new, recipe, clip = TRUE)
#
# Build the identified interaction design matrix at NEW data points,
# using the stored recipe from ls_build_interaction().
#
# Returns: n_new x (M_u-1)*(M_v-1) matrix W_uv_new
# ============================================================
ls_interaction_design_new <- function(X_u_new, X_v_new, recipe, clip = TRUE) {
  W_u_new <- recipe$design_new_u(X_u_new, type = "W", clip = clip)
  W_v_new <- recipe$design_new_v(X_v_new, type = "W", clip = clip)
  khatri_rao_rowwise_R(W_u_new, W_v_new)
}


# ============================================================
# ls_build_all_interactions(X, M_vec, objs)
#
# Build interaction designs for all p*(p-1)/2 covariate pairs.
# Useful after ls_additive_build() has already built objs.
#
# Inputs:
#   X      : n x p covariate matrix
#   objs   : list of length p, from ls_additive_build() (each has $W, $M, $T, etc.)
#            Note: ls_additive_build returns compact objs WITHOUT $W stored.
#            Pass the full ls_build_one_full() objects instead for p=2 POC.
#   M_vec  : integer vector of length p (same M for all if scalar)
#
# Returns:
#   A list indexed by pair "(u,v)":
#     $W_uv  : n x (M_u-1)*(M_v-1)
#     $K_uv  : identified 2D RW2 penalty
#     $recipe: for prediction
# ============================================================
ls_build_all_interactions <- function(full_objs) {
  p <- length(full_objs)
  if (p < 2) stop("Need at least 2 covariates for interactions.")
  pairs <- combn(p, 2, simplify = FALSE)
  int_list <- vector("list", length(pairs))
  names_list <- character(length(pairs))

  for (k in seq_along(pairs)) {
    u <- pairs[[k]][1]
    v <- pairs[[k]][2]
    int_list[[k]] <- ls_build_interaction(full_objs[[u]], full_objs[[v]])
    names_list[k] <- paste0(u, "_", v)
  }
  names(int_list) <- names_list
  int_list
}


# ============================================================
# Assemble full design matrix H for the model with interactions
#
# Model: y = [1 | W_main | W_interactions] eta + b + eps
#
# Inputs:
#   intercept_col : n x 1 column of ones
#   W_main        : n x sum(M_j-1)   main effect design (from ls_additive_build)
#   int_list      : list from ls_build_all_interactions()
#
# Returns:
#   H        : full n x (1 + sum(M_j-1) + sum((M_u-1)*(M_v-1))) design matrix
#   col_map  : column index map (for Gibbs sampler prior block structure)
#              List with: $intercept, $main[[j]], $interaction[[k]]
# ============================================================
ls_assemble_full_design <- function(W_main, int_list, n) {
  intercept <- matrix(1, n, 1)
  int_blocks <- lapply(int_list, function(x) x$W_uv)

  H <- cbind(intercept, W_main, do.call(cbind, int_blocks))

  # Build col_map (0-indexed from col 1 = intercept at position 0)
  p_main <- ncol(W_main)
  p_int_each <- sapply(int_blocks, ncol)

  col_map <- list()
  col_map$intercept <- 0L  # column index 1 in 1-indexed = 0 in offset

  # main effect blocks: cols 1..(p_main) in H (0-indexed: 1..p_main)
  col_map$main <- list()
  start <- 1L
  for (j in seq_len(ncol(W_main) > 0)) {
    # This should mirror the structure of how ls_additive_build sets col_map
    col_map$main[[j]] <- start  # placeholder; caller should merge with additive col_map
  }

  # interaction blocks
  col_map$interaction <- vector("list", length(int_list))
  int_start <- 1L + p_main  # after intercept and main effects (0-indexed offset)
  for (k in seq_along(int_list)) {
    col_map$interaction[[k]] <- int_start + seq_len(p_int_each[k]) - 1L
    int_start <- int_start + p_int_each[k]
  }

  list(H = H, col_map = col_map, p_main = p_main, p_int_each = p_int_each)
}


# ============================================================
# Sanity checks
# ============================================================
ls_interaction_tests <- function(n = 200, M = 8, seed = 42) {
  set.seed(seed)
  cat("=== ls_interaction sanity checks ===\n")

  X1 <- runif(n, 0, 1)
  X2 <- runif(n, 0, 1)

  bu <- ls_build_one_full(X1, M = M)
  bv <- ls_build_one_full(X2, M = M)

  int <- ls_build_interaction(bu, bv)

  # Dimension checks
  expected_cols <- (M - 1)^2
  stopifnot(nrow(int$W_uv) == n)
  stopifnot(ncol(int$W_uv) == expected_cols)
  cat(sprintf("  W_uv: %d x %d  [expected %d x %d] OK\n",
              nrow(int$W_uv), ncol(int$W_uv), n, expected_cols))

  # Penalty dimension
  stopifnot(nrow(int$K_uv) == expected_cols)
  stopifnot(ncol(int$K_uv) == expected_cols)
  cat(sprintf("  K_uv: %d x %d  OK\n", nrow(int$K_uv), ncol(int$K_uv)))

  # K_uv should be symmetric PSD
  ev <- eigen(int$K_uv, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(ev > -1e-8))
  cat(sprintf("  K_uv eigenvalues: min=%.2e, max=%.2e  (PSD OK)\n",
              min(ev), max(ev)))

  # Rank check: expected (M-1)^2 - 1 = (M-1)^2 - 1
  r <- qr(int$K_uv)$rank
  cat(sprintf("  K_uv rank: %d  [expected %d]\n", r, expected_cols - 1))

  # Identifiability: column sums of W_uv beta should match
  # W_uv = W_u ⊙ W_v; each row is a Kronecker of identified (sum-to-zero) rows
  # Check: W_uv %*% 1 should NOT be zero generally (it's data-dependent)
  # But: contrast enforced via T_uv — verify T_u sum-to-zero property
  T_u <- bu$T
  ones_Mu <- rep(1, M)
  stopifnot(max(abs(t(ones_Mu) %*% T_u)) < 1e-10)
  cat("  T_u sum-to-zero constraint: OK\n")

  # Prediction on new data
  X1_new <- runif(50, 0, 1)
  X2_new <- runif(50, 0, 1)
  W_uv_new <- ls_interaction_design_new(X1_new, X2_new, int$recipe)
  stopifnot(nrow(W_uv_new) == 50, ncol(W_uv_new) == expected_cols)
  cat(sprintf("  Prediction design: %d x %d  OK\n", nrow(W_uv_new), ncol(W_uv_new)))

  cat("=== All checks passed ===\n")
  invisible(int)
}
