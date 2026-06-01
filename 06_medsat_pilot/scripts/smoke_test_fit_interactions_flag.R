## ============================================================
## scripts/smoke_test_fit_interactions_flag.R
##
## Smoke test for the new `fit_interactions` switch on
## fit_ours_interaction. Tiny MCMC budget (100 iters total) in both
## modes. Verifies:
##   - additive mode (fit_interactions = FALSE) completes; f_int /
##     tau2_int posteriors are empty lists
##   - all-pairs mode (fit_interactions = TRUE)  completes as before
##   - additive mode is faster per iteration than all-pairs
## ============================================================

SIM_RDS <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"
stopifnot(file.exists(SIM_RDS))
sim <- readRDS(SIM_RDS)

# Compile C++ kernels for fair runtime comparison.
if (requireNamespace("Rcpp", quietly = TRUE)) {
  tryCatch(Rcpp::sourceCpp("R/ls_interaction_core.cpp", env = globalenv()),
           error = function(e) message("C++ compile failed: ",
                                       conditionMessage(e)))
}
stopifnot(exists("khatri_rao_cpp"))

env <- new.env(parent = globalenv())
sys.source("R/fit_ours_interaction_wrapper.R", chdir = TRUE, envir = env)

tiny <- list(
  M       = 10L,
  nu      = 1.0,
  n_iter  = 100L,
  n_burn  = 50L,
  n_thin  = 1L,
  n_draws = 50L,
  orthogonalize = TRUE,
  verbose = TRUE
)

# ---- inner-source() workaround for fit_spatial_reml.R ----
run_with_chdir <- function(fn, ...) {
  old <- getwd(); setwd("R"); on.exit(setwd(old))
  fn(...)
}

cat("\n=========== Mode A: fit_interactions = FALSE ===========\n")
t0 <- Sys.time()
fit_add <- run_with_chdir(env$fit_ours_interaction, sim,
  settings = c(tiny, list(fit_interactions = FALSE)))
t_add <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Mode A wall time: %.2fs (%.3fs / iter)\n",
            t_add, t_add / tiny$n_iter))

cat("\n=========== Mode B: fit_interactions = TRUE  ===========\n")
t0 <- Sys.time()
fit_full <- run_with_chdir(env$fit_ours_interaction, sim,
  settings = c(tiny, list(fit_interactions = TRUE)))
t_full <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Mode B wall time: %.2fs (%.3fs / iter)\n",
            t_full, t_full / tiny$n_iter))

# ----------------------------------------------------------------
# Assertions
# ----------------------------------------------------------------
cat("\n=========== Assertions ==================================\n")
ok <- TRUE
checkthat <- function(cond, msg) {
  if (!isTRUE(cond)) {
    cat("FAIL: ", msg, "\n", sep = ""); ok <<- FALSE
  } else {
    cat("OK  : ", msg, "\n", sep = "")
  }
}

# A: additive mode
checkthat(is.list(fit_add$f_int)             && length(fit_add$f_int)            == 0L,
          "Mode A: fit_add$f_int is empty list")
checkthat(is.list(fit_add$var_comp$tau2_int) && length(fit_add$var_comp$tau2_int) == 0L,
          "Mode A: fit_add$var_comp$tau2_int is empty list")
checkthat(isFALSE(fit_add$settings$fit_interactions),
          "Mode A: fit_add$settings$fit_interactions == FALSE")
checkthat(identical(fit_add$settings$int_keys, character(0)),
          "Mode A: fit_add$settings$int_keys == character(0)")
checkthat(t_add < 120,
          sprintf("Mode A: wall time %.2fs < 120s", t_add))

# B: all-pairs mode (backward-compat check)
expected_keys <- sim$int_keys
checkthat(is.list(fit_full$f_int) &&
            setequal(names(fit_full$f_int), expected_keys) &&
            length(fit_full$f_int) == length(expected_keys),
          sprintf("Mode B: fit_full$f_int has keys {%s}",
                  paste(expected_keys, collapse = ",")))
checkthat(is.list(fit_full$var_comp$tau2_int) &&
            setequal(names(fit_full$var_comp$tau2_int), expected_keys),
          "Mode B: tau2_int keys match sim$int_keys")
checkthat(isTRUE(fit_full$settings$fit_interactions),
          "Mode B: fit_full$settings$fit_interactions == TRUE")
checkthat(t_full < 120,
          sprintf("Mode B: wall time %.2fs < 120s", t_full))

# Per-iter relative speed
ratio <- (t_full / tiny$n_iter) / (t_add / tiny$n_iter)
checkthat(t_add < t_full,
          sprintf("Mode A per-iter (%.3fs) < Mode B per-iter (%.3fs);  speedup = %.2fx",
                  t_add / tiny$n_iter, t_full / tiny$n_iter, ratio))

cat("\n", if (ok) "ALL ASSERTIONS PASSED" else "FAILED",
    " (additive ", round(t_add, 1), "s, full ", round(t_full, 1), "s; ",
    "speedup ", round(ratio, 2), "x)\n", sep = "")

if (!ok) quit(status = 1L)
