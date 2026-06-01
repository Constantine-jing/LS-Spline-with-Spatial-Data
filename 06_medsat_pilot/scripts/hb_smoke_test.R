## ============================================================
## scripts/hb_smoke_test.R
##
## Hellbender environment smoke test. Run this FIRST after transferring
## the project to Hellbender to verify that:
##   1. R can find and load all required packages.
##   2. The C++ kernels (Rcpp/RcppEigen) compile and link.
##   3. The sampler completes a tiny run on the n=500 v2 sim.
##
## Uses M=10, n_iter=100, n_burn=50, n_thin=1, n_draws=50 with
## fit_interactions=TRUE. Expected wall time: < 3 minutes on any
## modern node. Run with:
##
##   Rscript scripts/hb_smoke_test.R
##
## A clean exit (status 0) means the environment is ready for the
## production run (scripts/run_n1500_allpairs.R).
## ============================================================

cat("============================================================\n")
cat(sprintf(" Hellbender smoke test  |  %s\n", Sys.time()))
cat(sprintf(" Hostname: %s\n", Sys.info()[["nodename"]]))
cat(sprintf(" R version: %s\n", R.version.string))
cat("============================================================\n\n")

SIM_RDS <- "output/sim_medsat_london_asthma_n500_seed1_v2.rds"

# ---- 1. Package availability ----
cat("--- 1. Package availability ---\n")
required_pkgs <- c("Rcpp", "RcppEigen", "coda")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs,
                                        requireNamespace,
                                        logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(sprintf("Missing packages: %s\n  Install with: install.packages(c(%s))",
               paste(missing_pkgs, collapse = ", "),
               paste(shQuote(missing_pkgs), collapse = ", ")))
}
cat("  OK: Rcpp, RcppEigen, coda found.\n")

# ---- 2. Load sim object ----
cat("\n--- 2. Load sim (n=500 v2) ---\n")
if (!file.exists(SIM_RDS)) {
  stop(sprintf("Sim RDS not found: %s\n  Transfer it to Hellbender first.", SIM_RDS))
}
sim <- readRDS(SIM_RDS)
stopifnot(sim$n == 500L, sim$p == 3L, identical(sim$truth_params$nu, 1.0))
cat(sprintf("  OK: sim loaded (n=%d  p=%d  nu=%g)\n", sim$n, sim$p, sim$truth_params$nu))

# ---- 3. Compile C++ kernels ----
cat("\n--- 3. Compile C++ kernels ---\n")
t_cpp <- system.time(
  tryCatch(
    Rcpp::sourceCpp("R/ls_interaction_core.cpp", env = globalenv()),
    error = function(e) stop("Rcpp::sourceCpp failed: ", conditionMessage(e))
  )
)
if (!exists("khatri_rao_cpp")) {
  stop("khatri_rao_cpp not found after sourceCpp -- compilation likely failed.")
}
cat(sprintf("  OK: C++ compiled in %.1f s\n", t_cpp[3]))

# ---- 4. Source wrapper ----
cat("\n--- 4. Source wrapper ---\n")
env <- new.env(parent = globalenv())
sys.source("R/fit_ours_interaction_wrapper.R", chdir = TRUE, envir = env)
stopifnot(exists("fit_ours_interaction", envir = env))
cat("  OK: fit_ours_interaction sourced\n")

# ---- 5. Tiny MCMC run ----
cat("\n--- 5. Tiny MCMC (M=10, 100 iters, fit_interactions=TRUE) ---\n")
tiny_settings <- list(
  M                = 10L,
  nu               = 1.0,
  n_iter           = 100L,
  n_burn           = 50L,
  n_thin           = 1L,
  n_draws          = 50L,
  orthogonalize    = TRUE,
  fit_interactions = TRUE,
  log_rho_mu       = log(0.2),
  log_rho_sd       = 1.0,
  verbose          = TRUE
)

run_with_chdir <- function(fn, ...) {
  old <- getwd(); setwd("R"); on.exit(setwd(old), add = TRUE)
  fn(...)
}

t0 <- proc.time()
fit <- run_with_chdir(env$fit_ours_interaction, sim, settings = tiny_settings)
elapsed <- (proc.time() - t0)[3]

# ---- 6. Assertions ----
cat("\n--- 6. Assertions ---\n")
ok <- TRUE
chk <- function(cond, msg) {
  if (!isTRUE(cond)) { cat("FAIL: ", msg, "\n", sep = ""); ok <<- FALSE }
  else                { cat("OK  : ", msg, "\n", sep = "") }
}

chk(is.list(fit),                               "fit is a list")
chk(length(fit$f_main)   == 3L,                "f_main has 3 elements")
chk(length(fit$f_int)    == 3L,                "f_int has 3 elements (all-pairs)")
chk(length(fit$var_comp$rho) == 50L,            "50 rho draws")
chk(all(fit$var_comp$rho > 0),                 "rho draws positive")
chk(all(fit$var_comp$sigma2 > 0),              "sigma2 draws positive")
chk(isTRUE(fit$settings$fit_interactions),     "fit_interactions = TRUE")
chk(isTRUE(fit$settings$orthogonalize),        "orthogonalize = TRUE")
chk(elapsed < 300,
    sprintf("wall time %.1f s < 300 s", elapsed))

# ---- 7. Summary ----
cat("\n============================================================\n")
if (ok) {
  cat(" SMOKE TEST PASSED\n")
} else {
  cat(" SMOKE TEST FAILED -- fix errors above before launching production run\n")
}
cat(sprintf(" Wall time:   %.1f s\n", elapsed))
cat(sprintf(" Per iter:    %.3f s  (expect << 1 s at M=10)\n",
            elapsed / tiny_settings$n_iter))
cat(sprintf(" rho mean:    %.4f\n", mean(fit$var_comp$rho)))
cat(sprintf(" sigma2 mean: %.4f\n", mean(fit$var_comp$sigma2)))
cat(sprintf(" MH accepts:  sigma2=%.2f  tau2=%.2f  rho=%.2f\n",
            fit$convergence$mh_accept[["sigma2"]],
            fit$convergence$mh_accept[["tau2"]],
            fit$convergence$mh_accept[["rho"]]))
cat(sprintf(" Hostname:    %s\n", Sys.info()[["nodename"]]))
cat("============================================================\n")

if (!ok) quit(status = 1L)
