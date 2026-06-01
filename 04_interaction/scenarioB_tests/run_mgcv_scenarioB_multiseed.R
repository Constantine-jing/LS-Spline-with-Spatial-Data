# ============================================================
# run_mgcv_scenarioB_multiseed.R
#
# Fit mgcv on all 10 Scenario B seeds so the comparison table
# against the post-fix multi-seed results is symmetric.
#
# mgcv fits in ~1 minute per seed -> ~10-15 min total.
#
# Output:
#   comparison/output/scenarioB_postfix_summary/mgcv_per_seed.csv
# ============================================================

rm(list = ls()); gc()

source("simulate_scenario_B.R")
library(mgcv)

out_dir <- file.path("comparison", "output", "scenarioB_postfix_summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seeds <- 1:10

rmse_curve <- function(pm, truth) {
  pm <- pm - mean(pm)
  tt <- truth - mean(truth)
  sqrt(mean((pm - tt)^2))
}

rows <- list()

for (sd in seeds) {
  cat(sprintf("Seed %d ... ", sd))
  sim <- simulate_scenario_B(seed = sd)
  df  <- sim$data
  xg  <- sim$x_grid_1d

  t0 <- Sys.time()
  fit_m <- gam(
    y ~ s(X1, bs = "ps", k = 20) +
        s(X2, bs = "ps", k = 20) +
        s(X3, bs = "ps", k = 20) +
        ti(X1, X2, bs = "ps", k = c(10, 10)) +
        ti(X1, X3, bs = "ps", k = c(10, 10)) +
        ti(X2, X3, bs = "ps", k = c(10, 10)) +
        s(lon, lat, bs = "gp"),
    data = df, method = "REML"
  )
  secs <- as.numeric(Sys.time() - t0, units = "secs")

  # main-effect curves on the 1d grid
  row <- list(seed = sd)
  for (j in 1:3) {
    vary <- xg
    dfp <- data.frame(
      X1  = if (j == 1) vary else mean(df$X1),
      X2  = if (j == 2) vary else mean(df$X2),
      X3  = if (j == 3) vary else mean(df$X3),
      lon = mean(df$lon),
      lat = mean(df$lat)
    )
    term <- sprintf("s(X%d)", j)
    pr <- predict(fit_m, newdata = dfp, type = "terms", terms = term)
    fj <- as.numeric(pr)
    row[[sprintf("RMSE_f%d", j)]] <- rmse_curve(fj, sim$truth_f_grid[[j]])
  }

  # interaction surfaces on the 2d grid
  flat2d <- sim$flat_grid_2d
  for (key in names(sim$truth_f_int)) {
    pair <- as.integer(strsplit(key, "_")[[1]])
    u <- pair[1]; v <- pair[2]
    dfp2 <- data.frame(
      X1 = mean(df$X1), X2 = mean(df$X2), X3 = mean(df$X3),
      lon = mean(df$lon), lat = mean(df$lat)
    )
    dfp2 <- dfp2[rep(1, nrow(flat2d)), ]
    dfp2[[sprintf("X%d", u)]] <- flat2d[, "u"]
    dfp2[[sprintf("X%d", v)]] <- flat2d[, "v"]
    term2 <- sprintf("ti(X%d,X%d)", u, v)
    pr2 <- tryCatch(
      as.numeric(predict(fit_m, newdata = dfp2, type = "terms", terms = term2)),
      error = function(e) rep(NA_real_, nrow(flat2d))
    )
    row[[sprintf("RMSE_f%s", key)]] <-
      if (all(is.na(pr2))) NA_real_
      else rmse_curve(pr2, sim$truth_f_int[[key]])
  }

  row$time_sec <- secs
  rows[[length(rows) + 1]] <- row
  cat(sprintf("done (%.0f s)  RMSE_f1=%.4f\n", secs, row$RMSE_f1))
}

df_mgcv <- do.call(rbind, lapply(rows, as.data.frame))
write.csv(df_mgcv, file.path(out_dir, "mgcv_per_seed.csv"), row.names = FALSE)

cat("\n=== mgcv per-seed ===\n")
print(round(df_mgcv, 4))

cat("\n=== mgcv summary (mean / se) ===\n")
mc <- setdiff(names(df_mgcv), "seed")
for (k in mc) {
  v <- df_mgcv[[k]]
  cat(sprintf("  %-12s mean=%.4f  se=%.4f\n",
              k, mean(v, na.rm = TRUE),
              sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))))
}

cat(sprintf("\nSaved: %s\n", file.path(out_dir, "mgcv_per_seed.csv")))
