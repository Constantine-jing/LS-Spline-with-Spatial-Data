# ============================================================
# diagnose_centering_bug.R
# 
# Hypothesis: per-draw recentering of f_main and f_int in
# fit_ours_interaction_wrapper.R is degrading recovery, especially
# for f1 (the largest-amplitude main effect under interaction
# collinearity).
# 
# We CAN'T un-do the centering on existing fits (the centered
# values are already saved). So we redo the per-smooth posterior
# extraction WITHOUT centering, using the saved gs object.
# 
# Wait -- the wrapper saves the centered f_main but NOT the raw
# eta_samples. We need to refit. UNLESS...
# 
# Actually we can! The centering subtracted colMeans which is
# a numeric subtraction. We just need the per-draw column means
# we removed, but we don't have them either.
# 
# Cleanest path: rerun the wrapper on seed=1 with centering
# DISABLED, and compare. That requires a small wrapper change.
# This file just patches the wrapper inline and reruns.
# ============================================================

source("simulate_scenario_B.R")
source("fit_ours_interaction_wrapper.R")

# Patch: replace the per-draw centering with a pass-through.
# We do this by overwriting fit_ours_interaction with a no-centering
# version. Cleanest is to copy the function and edit.

# Read the source of fit_ours_interaction:
src_fn <- deparse(fit_ours_interaction)

# Find and remove the two per-draw centering blocks. They look like:
#   col_means <- colMeans(f_main[[j]])
#   f_main[[j]] <- sweep(f_main[[j]], 2, col_means, FUN = "-")
# and similarly for f_int.

src <- paste(src_fn, collapse = "\n")
src <- gsub(
  pattern = "col_means <- colMeans\\(f_main\\[\\[j\\]\\]\\)\\s*\\n\\s*f_main\\[\\[j\\]\\] <- sweep\\(f_main\\[\\[j\\]\\], 2, col_means, FUN = \"-\"\\)",
  replacement = "# centering removed for diagnostic",
  x = src
)
src <- gsub(
  pattern = "col_means_2d <- colMeans\\(f_int\\[\\[k\\]\\]\\)\\s*\\n\\s*f_int\\[\\[k\\]\\]\\s*<- sweep\\(f_int\\[\\[k\\]\\], 2, col_means_2d, FUN = \"-\"\\)",
  replacement = "# centering removed for diagnostic",
  x = src
)

# Eval the patched function in a new env that knows about source helpers
patched_env <- new.env(parent = globalenv())
fit_ours_interaction_nocenter <- eval(parse(text = src), envir = patched_env)
environment(fit_ours_interaction_nocenter) <- globalenv()

# Sanity: check the patched body actually no longer contains "sweep"
body_str <- paste(deparse(fit_ours_interaction_nocenter), collapse = "\n")
if (grepl("sweep\\(f_main", body_str)) {
  warning("patched function STILL contains f_main centering -- regex failed!")
}
if (grepl("sweep\\(f_int", body_str)) {
  warning("patched function STILL contains f_int centering -- regex failed!")
}

# ---- run on seed=1 with the same settings as the failed REF case ----
sim <- simulate_scenario_B(seed = 1L, n = 500L)
cat(sprintf("[diag]  scenario B, seed=%d, n=%d, M=20 (matches REF), CENTERING DISABLED\n",
            sim$seed, sim$n))

cat("\n[fit]   running interaction Gibbs without per-draw centering...\n")
fit <- fit_ours_interaction_nocenter(sim)

cat("\n--- timing ---\n")
cat(sprintf("  total_sec = %.1f  (%.1f min)\n",
            fit$timing$total_sec, fit$timing$total_sec / 60))

cat("\n--- main-effect recovery WITHOUT per-draw centering ---\n")
for (j in seq_len(sim$p)) {
  truth <- sim$truth_f_grid[[j]]
  draws <- fit$f_main[[j]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  cat(sprintf("  f%d:  RMSE = %.3f   cov95 = %.2f\n", j, rmse, cov95))
}

cat("\n--- interaction recovery WITHOUT per-draw centering ---\n")
for (k in sim$int_keys) {
  truth <- sim$truth_f_int[[k]]
  draws <- fit$f_int[[k]]
  pmean <- rowMeans(draws)
  qq    <- apply(draws, 1, quantile, probs = c(0.025, 0.975))
  rmse  <- sqrt(mean((pmean - truth)^2))
  cov95 <- mean(truth >= qq[1, ] & truth <= qq[2, ])
  width <- mean(qq[2, ] - qq[1, ])
  is_null <- isTRUE(all(truth == 0))
  tag <- if (is_null) "[NULL]" else "[TRUE]"
  cat(sprintf("  f_%-3s %s:  RMSE = %.3f   cov95 = %.2f   width = %.3f\n",
              k, tag, rmse, cov95, width))
}

# ---- save ----
out_path <- file.path("..", "output", "fits",
                      sprintf("fit_ours_scenarioB_seed1_NOCENTER.rds"))
saveRDS(fit, out_path)

cat(sprintf("\n[save]  %s\n", normalizePath(out_path, mustWork = FALSE)))
cat("\n--- COMPARISON ---\n")
cat("REF (centered)   : f1 RMSE = 0.337, cov = 0.56\n")
cat("Anchor validated : f1 RMSE ~ 0.17 (but at lower noise + nu=1.5)\n")
cat("If diag <= 0.20  : centering bug confirmed.\n")
cat("If diag > 0.30   : centering not the (only) issue.\n")
