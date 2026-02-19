# ============================================================
# test_wald_selection.R
#
# Test the block Wald variable selection on Sim5-type data.
#
# Experiments:
#   (A) Single run at n=1000, 10 real + 10 garbage:
#       -> Print Wald table, check real vs garbage separation
#
#   (B) Garbage sweep: fix n=1000, 10 real covariates,
#       vary p_garbage in {10, 20, 50, 100}
#       -> Track: power (fraction of real covariates selected)
#                 FPR   (fraction of garbage covariates selected)
#
#   (C) Sample-size sweep: fix 10 real + 50 garbage,
#       vary n in {200, 400, 1000, 2000}
#       -> Same metrics
#
# All use BH-corrected p-values at alpha = 0.05.
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("sim5_sweep.R")
source("variable_selection.R")

library(mvtnorm)

# ============================================================
# (A) Single diagnostic run
# ============================================================
cat("\n========== (A) Single run: n=1000, 10 real + 10 garbage ==========\n\n")

n <- 1000; M <- 6; seed <- 42
p_real <- 10; p_garbage <- 10; p <- p_real + p_garbage

dat <- simulate_sim5(n = n, p_real = p_real, p_garbage = p_garbage,
                     sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15, seed = seed)

coords <- as.matrix(dat[, c("x", "y_coord")])
X_raw  <- as.matrix(dat[, paste0("X", 1:p)])
y      <- dat$Y

cat("Fitting model with", p, "covariates (", p_real, "real +", p_garbage, "garbage ) ...\n")
obj <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                      M_vec = rep(M, p), nu = 1.5,
                      rho_init = 0.2, lambda_init = 0.15/0.8, verbose = TRUE)

var_names <- paste0("X", 1:p)
sel <- select_covariates(obj, var_names = var_names, alpha = 0.05, method = "BH")

# Add a column indicating truth
sel$truth <- c(rep("real", p_real), rep("garbage", p_garbage))

cat("\n--- Full Wald table ---\n")
print(sel[, c("var", "truth", "F_stat", "p_value", "p_adjusted", "selected")], digits = 4)

cat("\nPower (real selected):  ",
    mean(sel$selected[sel$truth == "real"]), "\n")
cat("FPR (garbage selected):",
    mean(sel$selected[sel$truth == "garbage"]), "\n")

# Also print importance table for comparison
cat("\n--- Importance table (sd of fitted marginal) ---\n")
obj$X_raw_for_marginal <- X_raw
imp <- importance_table(obj, var_names = var_names)
imp$truth <- c(rep("real", p_real), rep("garbage", p_garbage))
print(imp, digits = 4)


# ============================================================
# (B) Garbage sweep: fix n=1000, vary p_garbage
# ============================================================
cat("\n\n========== (B) Garbage sweep: n=1000, p_real=10 ==========\n\n")

garbage_vec <- c(10, 20, 50, 100)
sweep_B <- data.frame()

for (pg in garbage_vec) {
  p_tot <- p_real + pg
  cat(sprintf("  p_garbage = %d (total p = %d) ... ", pg, p_tot))

  dat_b <- simulate_sim5(n = 1000, p_real = 10, p_garbage = pg,
                         sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15, seed = 42)

  coords_b <- as.matrix(dat_b[, c("x", "y_coord")])
  X_raw_b  <- as.matrix(dat_b[, paste0("X", 1:p_tot)])
  y_b      <- dat_b$Y

  # Check feasibility: need n > total fixed-effect columns = 1 + p_tot*(M-1)
  k_total <- 1 + p_tot * (M - 1)
  if (k_total >= 1000) {
    cat("SKIP (k =", k_total, ">= n = 1000)\n")
    next
  }

  obj_b <- fit_ls_spatial(y = y_b, X_raw = X_raw_b, coords = coords_b,
                          M_vec = rep(M, p_tot), nu = 1.5,
                          rho_init = 0.2, lambda_init = 0.15/0.8, verbose = FALSE)

  sel_b <- select_covariates(obj_b, var_names = paste0("X", 1:p_tot),
                              alpha = 0.05, method = "BH")
  sel_b$truth <- c(rep("real", 10), rep("garbage", pg))

  power_b <- mean(sel_b$selected[sel_b$truth == "real"])
  fpr_b   <- mean(sel_b$selected[sel_b$truth == "garbage"])

  # Which real covariates were missed?
  missed <- sel_b$var[sel_b$truth == "real" & !sel_b$selected]

  cat(sprintf("power = %.3f, FPR = %.3f", power_b, fpr_b))
  if (length(missed) > 0) cat("  [missed:", paste(missed, collapse = ", "), "]")
  cat("\n")

  sweep_B <- rbind(sweep_B, data.frame(
    n = 1000, p_real = 10, p_garbage = pg,
    k_total = k_total,
    power = power_b, FPR = fpr_b,
    n_missed = length(missed),
    missed = paste(missed, collapse = ",")
  ))
}

cat("\n--- Garbage sweep summary ---\n")
print(sweep_B, row.names = FALSE)


# ============================================================
# (C) Sample-size sweep: fix 10 real + 50 garbage, vary n
# ============================================================
cat("\n\n========== (C) Sample-size sweep: p_real=10, p_garbage=50 ==========\n\n")

n_vec_c <- c(400, 1000, 2000, 5000)
pg_c <- 50
p_tot_c <- 10 + pg_c
k_total_c <- 1 + p_tot_c * (M - 1)
cat("Fixed-effect columns k =", k_total_c, "\n\n")

sweep_C <- data.frame()

for (nn in n_vec_c) {
  if (nn <= k_total_c) {
    cat(sprintf("  n = %d: SKIP (n <= k = %d)\n", nn, k_total_c))
    next
  }

  cat(sprintf("  n = %d ... ", nn))

  dat_c <- simulate_sim5(n = nn, p_real = 10, p_garbage = pg_c,
                         sigma2 = 0.8, rho = 0.2, nu = 1.5, tau2 = 0.15, seed = 42)

  coords_c <- as.matrix(dat_c[, c("x", "y_coord")])
  X_raw_c  <- as.matrix(dat_c[, paste0("X", 1:p_tot_c)])
  y_c      <- dat_c$Y

  obj_c <- fit_ls_spatial(y = y_c, X_raw = X_raw_c, coords = coords_c,
                          M_vec = rep(M, p_tot_c), nu = 1.5,
                          rho_init = 0.2, lambda_init = 0.15/0.8, verbose = FALSE)

  sel_c <- select_covariates(obj_c, var_names = paste0("X", 1:p_tot_c),
                              alpha = 0.05, method = "BH")
  sel_c$truth <- c(rep("real", 10), rep("garbage", pg_c))

  power_c <- mean(sel_c$selected[sel_c$truth == "real"])
  fpr_c   <- mean(sel_c$selected[sel_c$truth == "garbage"])
  missed  <- sel_c$var[sel_c$truth == "real" & !sel_c$selected]

  cat(sprintf("power = %.3f, FPR = %.3f", power_c, fpr_c))
  if (length(missed) > 0) cat("  [missed:", paste(missed, collapse = ", "), "]")
  cat("\n")

  sweep_C <- rbind(sweep_C, data.frame(
    n = nn, p_real = 10, p_garbage = pg_c,
    k_total = k_total_c,
    power = power_c, FPR = fpr_c,
    n_missed = length(missed),
    missed = paste(missed, collapse = ",")
  ))
}

cat("\n--- Sample-size sweep summary ---\n")
print(sweep_C, row.names = FALSE)


cat("\n\nDone.\n")
