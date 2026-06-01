## ============================================================
## scripts/pilot_p3_inspect.R
##
## Quick non-destructive read of the n=1500 Hellbender fit and its
## matching sim object. Confirms structure, key dims, draw counts,
## ESS, and basic var-component summaries. Prints; does not save.
## ============================================================

FIT_RDS <- "output/fit_medsat_n1500_seed1_allpairs_hb.rds"
SIM_RDS <- "output/sim_medsat_london_asthma_n1500_seed1_v2.rds"
stopifnot(file.exists(FIT_RDS), file.exists(SIM_RDS))

fit <- readRDS(FIT_RDS)
sim <- readRDS(SIM_RDS)

cat("\n[sim]\n")
cat(sprintf("  n=%d  p=%d  scenario=%s\n", sim$n, sim$p, sim$scenario))
cat(sprintf("  int_keys: {%s}\n", paste(sim$int_keys, collapse = ", ")))
cat(sprintf("  X_scale: %s\n",
            paste(vapply(sim$X_scale, function(s)
              sprintf("%s[%.3f,%.3f]", s$name, s$lo, s$hi), character(1)),
              collapse = ", ")))
cat(sprintf("  coords_scale: center=(%.0f,%.0f) scale=%.0f\n",
            sim$coords_scale$center[1], sim$coords_scale$center[2],
            sim$coords_scale$scale))

cat("\n[fit] top-level names:\n")
print(names(fit))

cat("\n[fit$settings$fit_interactions]: ", fit$settings$fit_interactions, "\n")
cat("[fit$settings$M]: ",                  fit$settings$M, "\n")
cat("[fit$settings$n_draws]: ",            fit$settings$n_draws, "\n")
cat("[fit$settings$orthogonalize]: ",      fit$settings$orthogonalize, "\n")

cat("\n[fit] f_main: list of length", length(fit$f_main),
    "; per-element dims:\n")
for (j in seq_along(fit$f_main)) {
  cat(sprintf("  f_main[[%d]]: %s\n", j,
              paste(dim(fit$f_main[[j]]), collapse=" x ")))
}

cat("\n[fit] f_int: keys=", paste(names(fit$f_int), collapse=", "), "\n")
for (k in names(fit$f_int)) {
  cat(sprintf("  f_int[['%s']]: %s\n", k,
              paste(dim(fit$f_int[[k]]), collapse=" x ")))
}

cat("\n[fit] s_obs: ",
    if (is.null(fit$s_obs)) "NULL" else paste(dim(fit$s_obs), collapse=" x "),
    "\n")

cat("\n[fit$var_comp] names:\n")
print(names(fit$var_comp))

cat("\n[fit$convergence$mh_accept]:\n")
print(fit$convergence$mh_accept)

cat("\n[fit$timing]:\n")
print(fit$timing)

cat("\n[fit$reproducibility]:\n")
print(fit$reproducibility)

# ESS summary
if (requireNamespace("coda", quietly = TRUE)) {
  cat("\n[ESS] (out of", fit$settings$n_draws, "draws)\n")
  vc <- fit$var_comp
  ess_one <- function(x) as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
  rows <- list()
  rows[["sigma2"]]      <- ess_one(vc$sigma2)
  rows[["tau2_s_spat"]] <- ess_one(vc$tau2_s)
  rows[["rho"]]         <- ess_one(vc$rho)
  for (j in seq_along(vc$tau2))     rows[[paste0("tau2_s_main[",j,"]")]] <- ess_one(vc$tau2[[j]])
  for (k in names(vc$tau2_int))     rows[[paste0("tau2_int[",k,"]")]]    <- ess_one(vc$tau2_int[[k]])
  print(round(unlist(rows), 1))
}

cat("\n[posterior quick] sigma2 mean=", round(mean(fit$var_comp$sigma2),4),
    "  tau2_s mean=", round(mean(fit$var_comp$tau2_s),4),
    "  rho mean=",    round(mean(fit$var_comp$rho),5), "\n")

cat("\nDone.\n")
