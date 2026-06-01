# ============================================================
# prepare_medsat_london.R
#
# Build a sim-shaped object from the MEDSAT (Scepanovic et al., 2023)
# LSOA-level public-health x environment dataset, restricted to Greater
# London. The returned object plugs directly into
# fit_ours_interaction(sim, settings = ...) WITHOUT any wrapper changes
# beyond the small Option-A patch that tolerates NA truth_params.
#
# Real-data fields that have no analogue in simulation are set to NULL
# (truth_f_grid, truth_s_obs, truth_f_int, truth_f12_obs) or NA
# (truth_params$rho, $sigma2, $tau2_s). This is intentional: anyone
# reading the saved RDS can immediately see this is not a simulation,
# and any sim-only metric script (compare_seed1_metrics.R etc.) will
# fail loudly rather than silently treating estimates as truth.
#
# Usage:
#   source("prepare_medsat_london.R")
#   sim <- prepare_medsat_london(
#     medsat_path = "data/medsat/2019_spatial_raw_master.csv",
#     imd_path    = "data/medsat/imd_2019_lsoa21.csv",
#     outcome     = "asthma",
#     covariates  = c("NO2", "NDVI", "IMD"),
#     n_subsample = 500L,
#     seed        = 1L
#   )
#   saveRDS(sim, "output/sim_medsat_london_asthma_n500_seed1.rds")
#
# Then later, with no wrapper code changes beyond the Option-A patch:
#   fit <- fit_ours_interaction(sim, settings = list(n_iter = 3000,
#                                                    n_burn = 1000,
#                                                    n_thin = 1,
#                                                    n_draws = 2000,
#                                                    nu = 1.0))
#
# ------------------------------------------------------------
# IMPORTANT — input file expectations
# ------------------------------------------------------------
# Default column names match the real MEDSAT 2019 master CSV
# (data/medsat/2019_spatial_raw_master.csv):
#
#     "geography code"               -- LSOA21CD (note the space)
#     "LSOA21NM"                     -- LSOA name; borough is derived from
#                                       this when no 'borough' column exists,
#                                       via regex "\\s\\d{3}[A-Z]$"
#                                       (e.g. "Camden 001A" -> "Camden")
#     "o_asthma_quantity_per_capita" -- outcome
#     "e_NO2", "e_ndvi"              -- environmental covariates
#     "centroid_x", "centroid_y"     -- LSOA centroid, OSGB (metres)
#
# IMD 2019 is NOT in the MEDSAT master file. Supply it via `imd_path`
# pointing at the LSOA21-resolved IMD lookup produced by
# scripts/build_imd_lsoa21.R:
#
#     "LSOA21CD", "imd_score_2019", "imd_source_chgind"
#
# When "IMD" appears in `covariates` and the resolved IMD column is not
# already present in the master file, `imd_path` is required.
#
# If your local download has different column names, override via the
# `col_map` argument. Example:
#   prepare_medsat_london(medsat_path = "foo.csv",
#                         imd_path    = "imd.csv",
#                         col_map = list(outcome   = "asthma_2019",
#                                        NO2       = "no2_mean_2019",
#                                        NDVI      = "ndvi_mean_2019",
#                                        IMD       = "imd_score",
#                                        lsoa      = "LSOA11CD",
#                                        easting   = "X_centroid",
#                                        northing  = "Y_centroid",
#                                        borough   = "LAD21NM"))
#
# ============================================================

source("spatial_utils.R")   # pairdist() used downstream by fit_ours_interaction;
                            # sourcing here so test scripts have it loaded if they
                            # call the fitter on the returned object. Remove this
                            # line if you prefer to source it separately.

# ------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------

# Min-max scale a numeric vector to [0,1], returning the vector
# and the (min, max) so we can back-transform for plotting.
.scale_unit <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(".scale_unit: input contains non-finite values; impute or drop first.")
  }
  lo <- min(x); hi <- max(x)
  if (hi - lo < .Machine$double.eps) {
    stop(".scale_unit: input has zero range, cannot scale to [0,1].")
  }
  list(x = (x - lo) / (hi - lo), lo = lo, hi = hi)
}

# Scale 2D coordinates to the unit square jointly so aspect ratio is
# preserved (Matern depends on Euclidean distance and ρ-prior is
# calibrated to a unit-square domain).
.scale_coords_unit <- function(coords) {
  stopifnot(is.matrix(coords), ncol(coords) == 2)
  rng <- apply(coords, 2, range)
  span <- max(rng[2, ] - rng[1, ])    # joint scale: largest axis -> 1
  if (span < .Machine$double.eps) {
    stop(".scale_coords_unit: zero coordinate range.")
  }
  ctr <- rng[1, ]
  scaled <- sweep(coords, 2, ctr, FUN = "-") / span
  list(coords = scaled, center = ctr, scale = span)
}

# Read MEDSAT file in whatever format the user has on disk.
# check.names = FALSE so non-syntactic names like "geography code" survive.
.read_medsat_table <- function(path) {
  ext <- tools::file_ext(path)
  ext <- tolower(ext)
  if (ext == "csv") {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "tab")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "rds") {
    readRDS(path)
  } else {
    stop(sprintf("Unsupported file extension '%s'. Use .csv, .tsv, or .rds.", ext))
  }
}

# Parse borough (= LAD22NM equivalent) from LSOA21NM.
# LSOA21NM convention: "<LAD name> NNNL", e.g. "Barking and Dagenham 016A".
# Strips trailing space + 3 digits + 1 capital letter.
.parse_borough_from_lsoa21nm <- function(x) {
  sub("\\s\\d{3}[A-Z]$", "", as.character(x))
}

# OSGB easting/northing detector. Heuristic: OSGB easting is roughly
# 0..700000 metres, northing is 0..1300000. Lat/lon would be tiny by
# comparison. We use this only to warn, not to fail.
.detect_coord_system <- function(coords) {
  rng <- range(coords, na.rm = TRUE)
  if (rng[2] > 1000) "projected" else "geographic"
}

# Default column-name map for the real MEDSAT 2019 master CSV.
# Canonical name -> actual column in the master file. Overridable via
# `col_map`. The "borough" column doesn't exist in MEDSAT; it is derived
# from LSOA21NM inside prepare_medsat_london() and then accessed via
# cm$borough = "borough".
.default_col_map <- function(outcome, covariates) {
  outcome_default_for <- list(
    asthma       = "o_asthma_quantity_per_capita",
    diabetes     = "o_diabetes_quantity_per_capita",
    hypertension = "o_hypertension_quantity_per_capita",
    depression   = "o_depression_quantity_per_capita",
    anxiety      = "o_anxiety_quantity_per_capita",
    opioids      = "o_opioids_quantity_per_capita",
    total        = "o_total_quantity_per_capita"
  )
  covariate_default_for <- list(
    NO2  = "e_NO2",
    NDVI = "e_ndvi",
    IMD  = "imd_score_2019"
  )
  cm <- list(
    outcome  = outcome_default_for[[outcome]] %||% outcome,
    lsoa     = "geography code",
    easting  = "centroid_x",
    northing = "centroid_y",
    borough  = "borough"
  )
  for (cv in covariates) cm[[cv]] <- covariate_default_for[[cv]] %||% cv
  cm
}

# null-coalesce helper (R has no `%||%` until 4.4; define our own).
`%||%` <- function(a, b) if (is.null(a)) b else a


# ------------------------------------------------------------
# Stratified subsample by borough (proportional to borough sizes).
# Falls back to simple random subsampling if no borough column.
# ------------------------------------------------------------
.stratified_subsample <- function(df, n_sub, stratify_col = "borough",
                                  seed = 1L) {
  set.seed(seed)
  if (is.null(stratify_col) || !(stratify_col %in% names(df))) {
    idx <- sort(sample.int(nrow(df), n_sub))
    return(list(idx = idx, strata = NULL))
  }
  strata <- df[[stratify_col]]
  tab    <- table(strata)
  # Proportional allocation (largest-remainder rounding so they sum to n_sub)
  raw_alloc <- as.numeric(tab) * n_sub / sum(tab)
  base     <- floor(raw_alloc)
  remainder <- n_sub - sum(base)
  if (remainder > 0L) {
    frac_order <- order(raw_alloc - base, decreasing = TRUE)
    base[frac_order[seq_len(remainder)]] <- base[frac_order[seq_len(remainder)]] + 1L
  }
  names(base) <- names(tab)

  # Cap allocation at the actual stratum size (small boroughs).
  base <- pmin(base, as.integer(tab))
  # Top up if capping reduced the total.
  while (sum(base) < n_sub) {
    deficit <- n_sub - sum(base)
    room    <- as.integer(tab) - base
    add     <- pmin(room, deficit)
    if (all(add == 0L)) break
    # Add greedily to the strata with the most headroom
    add_order <- order(room, decreasing = TRUE)
    take <- min(deficit, sum(room))
    for (i in add_order) {
      if (take <= 0L) break
      step <- min(room[i], take)
      base[i] <- base[i] + step
      room[i] <- room[i] - step
      take <- take - step
    }
  }

  picked <- integer(0)
  for (b in names(base)) {
    n_b   <- base[[b]]
    if (n_b <= 0L) next
    rows  <- which(strata == b)
    chosen <- if (length(rows) <= n_b) rows else sample(rows, n_b)
    picked <- c(picked, chosen)
  }
  picked <- sort(picked)
  list(idx = picked, strata = strata[picked])
}


# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------
prepare_medsat_london <- function(
    medsat_path,
    imd_path            = NULL,
    outcome             = "asthma",
    covariates          = c("NO2", "NDVI", "IMD"),
    int_pairs           = list(c(1L, 2L)),   # which interactions to FIT
    n_subsample         = 500L,
    seed                = 1L,
    nu                  = 1.0,
    n_grid_1d           = 101L,
    n_grid_2d_per       = 30L,
    response_transform  = c("identity", "log", "sqrt"),
    london_filter       = c("borough_list", "lad_code_prefix", "none"),
    drop_ndvi_ceiling   = TRUE,
    # ^ MEDSAT's `e_ndvi` saturates at exactly 1.000 for a substantial
    # fraction of LSOAs (~9% of London) — a Sentinel-2 aggregation
    # artifact, not a real continuous greenness value. When TRUE, rows
    # with the NDVI column equal to 1.000 exact are dropped AFTER the
    # London filter and the NA drop, BEFORE subsampling. The count is
    # always recorded in sim$meta$ndvi_ceiling_dropped (0 when FALSE
    # or when no NDVI ceiling rows are present).
    col_map             = NULL,
    verbose             = TRUE
) {
  response_transform <- match.arg(response_transform)
  london_filter      <- match.arg(london_filter)
  p <- length(covariates)
  stopifnot(p >= 2L, n_subsample >= 50L)

  # -- merge user col_map with defaults --
  cm <- .default_col_map(outcome, covariates)
  if (!is.null(col_map)) {
    for (nm in names(col_map)) cm[[nm]] <- col_map[[nm]]
  }

  if (verbose) cat(sprintf("[prepare_medsat_london] reading %s\n", medsat_path))
  df <- .read_medsat_table(medsat_path)

  # ---- derive borough from LSOA21NM if no borough column ----
  if (!(cm$borough %in% names(df))) {
    if ("LSOA21NM" %in% names(df)) {
      df[["borough"]] <- .parse_borough_from_lsoa21nm(df[["LSOA21NM"]])
      cm$borough <- "borough"
      if (verbose) cat("  derived 'borough' from LSOA21NM\n")
    } else if (london_filter == "borough_list") {
      stop(sprintf(
        "london_filter='borough_list' requires either a '%s' column or an ",
        cm$borough),
        "'LSOA21NM' column to derive borough from. ",
        "Pass london_filter='lad_code_prefix' if you have LAD21CD, ",
        "or 'none' if you have pre-filtered upstream.")
    } else {
      if (verbose) cat("  no borough column and no LSOA21NM; subsampling will be unstratified.\n")
    }
  }

  # ---- IMD join (if IMD is a covariate and the score column is not in df) ----
  if ("IMD" %in% covariates && !(cm[["IMD"]] %in% names(df))) {
    if (is.null(imd_path)) {
      stop(sprintf(paste0(
        "IMD column '%s' is not in %s, and `imd_path` was not supplied.\n",
        "  The MEDSAT 2019 master file does not carry IMD 2019; build the\n",
        "  LSOA21-resolved lookup with scripts/build_imd_lsoa21.R, then pass\n",
        "  its CSV path via imd_path=. Alternatively, drop \"IMD\" from\n",
        "  `covariates`, or override col_map$IMD to a column that IS in the\n",
        "  master file."),
        cm[["IMD"]], medsat_path))
    }
    if (!(cm$lsoa %in% names(df))) {
      stop(sprintf(
        "Cannot join IMD: LSOA key column '%s' not found in %s.",
        cm$lsoa, medsat_path))
    }
    if (verbose) cat(sprintf("  joining IMD from %s\n", imd_path))
    imd_tbl <- utils::read.csv(imd_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
    if (!all(c("LSOA21CD", "imd_score_2019") %in% names(imd_tbl))) {
      stop(sprintf(
        "imd_path file '%s' must have columns 'LSOA21CD' and 'imd_score_2019' ",
        imd_path),
        "(typically also 'imd_source_chgind'). Got: ",
        paste(names(imd_tbl), collapse = ", "))
    }
    df_before <- nrow(df)
    df <- merge(df, imd_tbl, by.x = cm$lsoa, by.y = "LSOA21CD",
                all.x = FALSE, sort = FALSE)
    if (verbose) {
      cat(sprintf("  IMD inner-join: %d -> %d rows (%d dropped: no IMD match)\n",
                  df_before, nrow(df), df_before - nrow(df)))
    }
  }

  # ---- column-existence check ----
  needed <- c(cm$outcome, vapply(covariates, function(cv) cm[[cv]], character(1)))
  has_proj <- all(c(cm$easting,  cm$northing) %in% names(df))
  has_geo  <- all(c("longitude", "latitude")  %in% names(df))
  if (!has_proj && !has_geo) {
    stop(sprintf("Coordinates missing: need either (%s, %s) or (longitude, latitude).",
                 cm$easting, cm$northing))
  }
  needed <- c(needed,
              if (has_proj) c(cm$easting, cm$northing) else c("longitude", "latitude"))
  miss <- setdiff(needed, names(df))
  if (length(miss) > 0L) {
    stop(sprintf("Missing required columns in %s: %s",
                 medsat_path, paste(miss, collapse = ", ")))
  }

  # ---- LSOA-key column (informational only; not used for joining downstream) ----
  if (!(cm$lsoa %in% names(df))) {
    if (verbose) cat(sprintf("  note: LSOA key column '%s' not found; continuing.\n",
                             cm$lsoa))
  }

  # ---- restrict to London ----
  if (london_filter == "borough_list") {
    london_boroughs <- .london_boroughs()
    if (!(cm$borough %in% names(df))) {
      stop(sprintf(
        "london_filter='borough_list' but column '%s' not found. ",
        cm$borough),
        "Pass london_filter='lad_code_prefix' if you have LAD21CD, ",
        "or 'none' if you have already filtered upstream.")
    }
    keep <- df[[cm$borough]] %in% london_boroughs
    df <- df[keep, , drop = FALSE]
  } else if (london_filter == "lad_code_prefix") {
    if (!("LAD21CD" %in% names(df))) {
      stop("london_filter='lad_code_prefix' requires a LAD21CD column.")
    }
    # London LAD codes all start with E090000XX where XX in 01..33
    london_lads <- sprintf("E090000%02d", 1:33)
    df <- df[df$LAD21CD %in% london_lads, , drop = FALSE]
  }
  if (nrow(df) < n_subsample) {
    stop(sprintf(
      "After London filter, only %d LSOAs remain but n_subsample=%d requested.",
      nrow(df), n_subsample))
  }
  if (verbose) {
    cat(sprintf("  after London filter: %d LSOAs\n", nrow(df)))
  }

  # ---- NA diagnostics computed BEFORE listwise deletion ----
  use_cols <- c(cm$outcome, vapply(covariates, function(cv) cm[[cv]], character(1)),
                if (has_proj) c(cm$easting, cm$northing) else c("longitude", "latitude"),
                if (cm$borough %in% names(df)) cm$borough else NULL)
  complete  <- stats::complete.cases(df[, use_cols, drop = FALSE])
  n_dropped <- sum(!complete)

  na_by_col <- vapply(use_cols, function(c) sum(is.na(df[[c]])), integer(1))
  names(na_by_col) <- use_cols

  if (cm$borough %in% names(df)) {
    br        <- df[[cm$borough]]
    totals    <- as.data.frame(table(borough = br),
                               responseName = "n_total",
                               stringsAsFactors = FALSE)
    if (n_dropped > 0L) {
      dropped_b <- as.data.frame(table(borough = br[!complete]),
                                 responseName = "n_dropped",
                                 stringsAsFactors = FALSE)
      na_by_borough <- merge(totals, dropped_b, by = "borough", all.x = TRUE)
    } else {
      na_by_borough <- totals
      na_by_borough$n_dropped <- 0L
    }
    na_by_borough$n_dropped[is.na(na_by_borough$n_dropped)] <- 0L
    na_by_borough$pct_dropped <- round(
      100 * na_by_borough$n_dropped / na_by_borough$n_total, 2)
    na_by_borough <- na_by_borough[order(-na_by_borough$pct_dropped), ,
                                   drop = FALSE]
    rownames(na_by_borough) <- NULL
  } else {
    na_by_borough <- data.frame(
      borough     = character(0),
      n_total     = integer(0),
      n_dropped   = integer(0),
      pct_dropped = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # ---- drop rows with missing values in the columns we need ----
  if (n_dropped > 0L && verbose) {
    cat(sprintf("  dropping %d LSOAs with NA in required columns\n", n_dropped))
  }
  df <- df[complete, , drop = FALSE]
  if (nrow(df) < n_subsample) {
    stop(sprintf("After NA drop, only %d LSOAs remain but n_subsample=%d requested.",
                 nrow(df), n_subsample))
  }

  # ---- drop NDVI saturation ceiling (Sentinel-2 aggregation artifact) ----
  # Always count; only drop if drop_ndvi_ceiling is TRUE.
  ndvi_col <- if ("NDVI" %in% covariates) cm[["NDVI"]] else NA_character_
  n_ndvi_ceiling_dropped <- 0L
  if (!is.na(ndvi_col) && (ndvi_col %in% names(df))) {
    is_ceiling <- !is.na(df[[ndvi_col]]) & df[[ndvi_col]] == 1.0
    n_at_ceiling <- sum(is_ceiling)
    if (drop_ndvi_ceiling && n_at_ceiling > 0L) {
      df <- df[!is_ceiling, , drop = FALSE]
      n_ndvi_ceiling_dropped <- n_at_ceiling
      if (verbose) {
        cat(sprintf("  dropping %d LSOAs with %s == 1.000 (Sentinel-2 saturation artifact)\n",
                    n_at_ceiling, ndvi_col))
      }
    } else if (n_at_ceiling > 0L && verbose) {
      cat(sprintf("  note: %d LSOAs have %s == 1.000 (drop_ndvi_ceiling=FALSE, kept)\n",
                  n_at_ceiling, ndvi_col))
    }
    if (nrow(df) < n_subsample) {
      stop(sprintf("After NDVI-ceiling drop, only %d LSOAs remain but n_subsample=%d requested.",
                   nrow(df), n_subsample))
    }
  }

  # ---- subsample ----
  strat_col <- if (cm$borough %in% names(df)) cm$borough else NULL
  ss <- .stratified_subsample(df, n_sub = n_subsample,
                              stratify_col = strat_col, seed = seed)
  df_sub <- df[ss$idx, , drop = FALSE]
  if (verbose) {
    cat(sprintf("  subsampled to n=%d", nrow(df_sub)))
    if (!is.null(ss$strata)) {
      cat(sprintf(" across %d strata\n", length(unique(ss$strata))))
    } else {
      cat(" (no stratification)\n")
    }
  }

  # ---- extract & transform response ----
  y_raw <- df_sub[[cm$outcome]]
  y_transformed <- switch(
    response_transform,
    identity = y_raw,
    log      = {
      if (any(y_raw <= 0)) stop("log transform: response has non-positive values.")
      log(y_raw)
    },
    sqrt     = {
      if (any(y_raw < 0))  stop("sqrt transform: response has negative values.")
      sqrt(y_raw)
    }
  )

  # ---- extract & scale covariates to [0,1] ----
  X_raw <- vapply(covariates,
                  function(cv) as.numeric(df_sub[[cm[[cv]]]]),
                  numeric(nrow(df_sub)))
  X_raw <- matrix(X_raw, nrow = nrow(df_sub), ncol = p,
                  dimnames = list(NULL, covariates))
  X_scaled <- matrix(NA_real_, nrow = nrow(df_sub), ncol = p)
  X_scale_meta <- vector("list", p)
  for (j in seq_len(p)) {
    s <- .scale_unit(X_raw[, j])
    X_scaled[, j] <- s$x
    X_scale_meta[[j]] <- list(name = covariates[j], lo = s$lo, hi = s$hi)
  }
  colnames(X_scaled) <- paste0("X", seq_len(p))

  # ---- coordinates: scale jointly to unit square ----
  if (has_proj) {
    coords_raw <- as.matrix(df_sub[, c(cm$easting, cm$northing), drop = FALSE])
    coord_system <- "OSGB"
  } else {
    coords_raw <- as.matrix(df_sub[, c("longitude", "latitude"), drop = FALSE])
    coord_system <- "WGS84"
  }
  colnames(coords_raw) <- c("east", "north")
  cs_det <- .detect_coord_system(coords_raw)
  if (cs_det == "geographic" && coord_system == "OSGB") {
    warning("Easting/northing columns look like degrees, not metres. ",
            "Check your input file.")
  }
  cscaled <- .scale_coords_unit(coords_raw)
  # In the sim convention, the coordinate columns of $data are named lon, lat
  # (this is what the wrapper expects).
  coords_unit <- cscaled$coords
  colnames(coords_unit) <- c("lon", "lat")

  # ---- assemble sim$data ----
  sim_data <- data.frame(X_scaled,
                         lon = coords_unit[, "lon"],
                         lat = coords_unit[, "lat"],
                         y   = as.numeric(y_transformed))

  # ---- evaluation grids ----
  x_grid_1d <- seq(0, 1, length.out = n_grid_1d)
  u_grid    <- seq(0, 1, length.out = n_grid_2d_per)
  v_grid    <- seq(0, 1, length.out = n_grid_2d_per)
  flat2d    <- as.matrix(expand.grid(u = u_grid, v = v_grid,
                                     KEEP.OUT.ATTRS = FALSE))

  # ---- int_keys from the requested pairs ----
  int_keys <- vapply(int_pairs, function(pr) {
    pr <- sort(as.integer(pr))
    stopifnot(length(pr) == 2L, pr[1] >= 1L, pr[2] <= p, pr[1] < pr[2])
    paste(pr, collapse = "_")
  }, character(1))
  # ls_build_all_interactions builds ALL p*(p-1)/2 pairs, so int_keys must
  # match what the wrapper will produce. For now we enforce all-pairs; the
  # int_pairs argument is recorded but ignored at fitting time. (Filtering
  # to a subset of pairs would require a small change in the wrapper.)
  all_keys <- character(0)
  for (a in 1:(p - 1)) for (b in (a + 1):p) {
    all_keys <- c(all_keys, paste(a, b, sep = "_"))
  }
  if (!setequal(int_keys, all_keys) && verbose) {
    cat(sprintf("  note: int_pairs requested {%s} but wrapper fits all pairs {%s};\n",
                paste(int_keys, collapse = ","), paste(all_keys, collapse = ",")))
    cat("        sim$int_keys set to ALL pairs to match.\n")
  }
  int_keys <- all_keys

  # ---- record low-y rows (y_raw < 1.0) for honest reporting ----
  low_idx <- which(y_raw < 1.0)
  get_raw <- function(name) {
    col <- cm[[name]]
    if (!is.null(col) && col %in% names(df_sub)) {
      as.numeric(df_sub[[col]])[low_idx]
    } else {
      rep(NA_real_, length(low_idx))
    }
  }
  low_y_rows <- data.frame(
    borough          = if (!is.null(ss$strata)) {
                         as.character(ss$strata)[low_idx]
                       } else {
                         rep(NA_character_, length(low_idx))
                       },
    LSOA21NM         = if ("LSOA21NM" %in% names(df_sub)) {
                         df_sub[["LSOA21NM"]][low_idx]
                       } else {
                         rep(NA_character_, length(low_idx))
                       },
    y_raw            = y_raw[low_idx],
    NO2_raw          = get_raw("NO2"),
    NDVI_raw         = get_raw("NDVI"),
    IMD_raw          = get_raw("IMD"),
    stringsAsFactors = FALSE
  )

  # ---- meta for back-transforming plots + honest-reporting diagnostics ----
  meta <- list(
    outcome_name        = outcome,
    outcome_raw_summary = summary(y_raw),
    response_transform  = response_transform,
    covariates          = covariates,
    coord_system        = coord_system,
    n_pre_subsample     = nrow(df),
    london_filter       = london_filter,
    na_drops            = list(
      total      = as.integer(n_dropped),
      by_column  = na_by_col,
      by_borough = na_by_borough
    ),
    ndvi_ceiling_dropped = as.integer(n_ndvi_ceiling_dropped),
    low_y_rows           = low_y_rows
  )

  # ---- assemble sim object ----
  sim <- list(
    scenario      = "medsat_london",
    seed          = as.integer(seed),
    n             = nrow(df_sub),
    p             = p,

    data          = sim_data,

    # ---- REAL DATA: no ground truth ----
    truth_f_grid  = NULL,
    truth_s_obs   = NULL,
    truth_f_int   = NULL,
    truth_f12_obs = NULL,

    truth_params  = list(
      sigma2 = NA_real_,
      tau2   = NA_real_,
      tau2_s = NA_real_,
      rho    = NA_real_,
      nu     = as.numeric(nu),
      c_int  = NA_real_
    ),

    x_grid_1d     = x_grid_1d,
    x_grid_2d     = list(u = u_grid, v = v_grid),
    flat_grid_2d  = flat2d,

    int_keys      = int_keys,

    # ---- real-data-only metadata, kept on the object for plotting/back-transform ----
    coords_raw    = coords_raw,
    coords_scale  = list(center = cscaled$center, scale = cscaled$scale),
    X_scale       = X_scale_meta,
    borough       = if (!is.null(ss$strata)) as.character(ss$strata) else NULL,
    y_raw         = y_raw,
    meta          = meta,

    settings      = list(
      medsat_path        = medsat_path,
      imd_path           = imd_path,
      outcome            = outcome,
      covariates         = covariates,
      int_pairs          = int_pairs,
      n_subsample        = n_subsample,
      seed               = as.integer(seed),
      nu                 = as.numeric(nu),
      n_grid_1d          = n_grid_1d,
      n_grid_2d_per      = n_grid_2d_per,
      response_transform = response_transform,
      london_filter      = london_filter,
      drop_ndvi_ceiling  = drop_ndvi_ceiling,
      col_map            = cm
    )
  )

  if (verbose) {
    cat(sprintf("[prepare_medsat_london] done: n=%d  p=%d  int_keys={%s}\n",
                sim$n, sim$p, paste(sim$int_keys, collapse = ", ")))
    cat(sprintf("  y (%s, %s-transformed): mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
                outcome, response_transform,
                mean(sim$data$y), sd(sim$data$y),
                min(sim$data$y), max(sim$data$y)))
    for (j in seq_len(p)) {
      raw_rng <- c(X_scale_meta[[j]]$lo, X_scale_meta[[j]]$hi)
      cat(sprintf("  %s (X%d): original range=[%.3f, %.3f]  -> scaled [0,1]\n",
                  covariates[j], j, raw_rng[1], raw_rng[2]))
    }
    cat(sprintf("  coords: %s scaled to unit square (joint scale=%.0f)\n",
                coord_system, cscaled$scale))
  }

  sim
}


# ------------------------------------------------------------
# Lookup: 33 London borough names (canonical ONS spellings).
# Used by london_filter = 'borough_list'.
# ------------------------------------------------------------
.london_boroughs <- function() {
  c(
    "Barking and Dagenham", "Barnet", "Bexley", "Brent", "Bromley",
    "Camden", "City of London", "Croydon", "Ealing", "Enfield",
    "Greenwich", "Hackney", "Hammersmith and Fulham", "Haringey",
    "Harrow", "Havering", "Hillingdon", "Hounslow", "Islington",
    "Kensington and Chelsea", "Kingston upon Thames", "Lambeth",
    "Lewisham", "Merton", "Newham", "Redbridge",
    "Richmond upon Thames", "Southwark", "Sutton", "Tower Hamlets",
    "Waltham Forest", "Wandsworth", "Westminster"
  )
}


# ------------------------------------------------------------
# Quick smoke test: construct a tiny synthetic master CSV that mirrors
# the real MEDSAT 2019 column conventions (incl. space-named LSOA code,
# LSOA21NM suffix pattern, no borough column, no IMD column), plus a
# synthetic LSOA21 IMD CSV, and round-trip them through the prep
# function. Run with:
#   source("prepare_medsat_london.R"); prepare_medsat_london_test()
# ------------------------------------------------------------
prepare_medsat_london_test <- function() {
  set.seed(123)
  boroughs <- .london_boroughs()
  # 30 LSOAs per borough so every borough is well-represented and
  # subsampling stratifies cleanly.
  per_borough <- 30L
  fake <- do.call(rbind, lapply(boroughs, function(b) {
    suffixes <- sprintf("%03d%s", seq_len(per_borough),
                        sample(LETTERS, per_borough, replace = TRUE))
    data.frame(
      `geography code` = sprintf("E01%06d",
                                 sample.int(999999L, per_borough)),
      LSOA21NM         = paste(b, suffixes),
      o_asthma_quantity_per_capita =
        pmax(0.1, rnorm(per_borough, mean = 1.0, sd = 0.4)),
      e_NO2            = runif(per_borough, 1e-4, 2e-4),
      e_ndvi           = runif(per_borough, 0.05, 0.85),
      centroid_x       = runif(per_borough, 500000, 560000),
      centroid_y       = runif(per_borough, 160000, 200000),
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  }))
  # Force unique LSOA21CDs after sampling
  fake[["geography code"]] <- sprintf("E01%06d", seq_len(nrow(fake)))

  # Inject NDVI = 1.000 ceiling rows mimicking the Sentinel-2 saturation
  # artifact. ~10% of rows, fixed seed so the count is reproducible.
  set.seed(456)
  n_ceiling_inject <- as.integer(round(0.10 * nrow(fake)))   # 99 rows out of 990
  ceiling_idx <- sample.int(nrow(fake), n_ceiling_inject)
  fake$e_ndvi[ceiling_idx] <- 1.000

  tmp_master <- tempfile(fileext = ".csv")
  utils::write.csv(fake, tmp_master, row.names = FALSE)

  # Synthetic IMD lookup: one score per LSOA21CD in the master.
  imd_fake <- data.frame(
    LSOA21CD          = fake[["geography code"]],
    imd_score_2019    = runif(nrow(fake), 0.5, 90),
    imd_source_chgind = "U",
    stringsAsFactors  = FALSE
  )
  tmp_imd <- tempfile(fileext = ".csv")
  utils::write.csv(imd_fake, tmp_imd, row.names = FALSE)

  sim <- prepare_medsat_london(
    medsat_path       = tmp_master,
    imd_path          = tmp_imd,
    outcome           = "asthma",
    covariates        = c("NO2", "NDVI", "IMD"),
    n_subsample       = 500L,
    seed              = 1L,
    drop_ndvi_ceiling = TRUE
  )

  stopifnot(sim$n == 500L, sim$p == 3L)
  stopifnot(all(c("X1", "X2", "X3", "lon", "lat", "y") %in% names(sim$data)))
  stopifnot(min(sim$data$X1) >= 0, max(sim$data$X1) <= 1)
  stopifnot(min(sim$data$lon) >= 0, max(sim$data$lon) <= 1)
  stopifnot(is.null(sim$truth_f_grid), is.null(sim$truth_s_obs))
  stopifnot(is.na(sim$truth_params$rho), !is.na(sim$truth_params$nu))
  stopifnot(setequal(sim$int_keys, c("1_2", "1_3", "2_3")))

  # Borough parsing: every subsample row should carry a recognised borough.
  stopifnot(!is.null(sim$borough), length(sim$borough) == sim$n)
  parsed <- sort(unique(sim$borough))
  stopifnot(identical(parsed, sort(boroughs)),
            length(parsed) == 33L)

  # Settings should record both paths and the flag.
  stopifnot(identical(sim$settings$medsat_path, tmp_master),
            identical(sim$settings$imd_path,    tmp_imd),
            isTRUE(sim$settings$drop_ndvi_ceiling))

  # NDVI ceiling: dropped count matches injection, no surviving row hits 1.0,
  # and the raw NDVI max in the kept pool is strictly < 1.0.
  stopifnot(identical(sim$meta$ndvi_ceiling_dropped, n_ceiling_inject))
  stopifnot(sim$X_scale[[2]]$hi < 1.0)

  # na_drops shape: total integer, by_column named integer vector, by_borough
  # is a data.frame with the right four columns.
  stopifnot(is.list(sim$meta$na_drops),
            all(c("total", "by_column", "by_borough") %in% names(sim$meta$na_drops)),
            is.integer(sim$meta$na_drops$total),
            is.integer(sim$meta$na_drops$by_column),
            !is.null(names(sim$meta$na_drops$by_column)),
            is.data.frame(sim$meta$na_drops$by_borough),
            all(c("borough", "n_total", "n_dropped", "pct_dropped") %in%
                names(sim$meta$na_drops$by_borough)))
  # Synthetic data has no NAs.
  stopifnot(sim$meta$na_drops$total == 0L,
            all(sim$meta$na_drops$by_column == 0L))

  # low_y_rows shape: data.frame with exactly these columns.
  stopifnot(is.data.frame(sim$meta$low_y_rows),
            identical(names(sim$meta$low_y_rows),
                      c("borough","LSOA21NM","y_raw",
                        "NO2_raw","NDVI_raw","IMD_raw")))
  # Every low_y row's y_raw must be < 1.0 (vacuously true if empty).
  if (nrow(sim$meta$low_y_rows) > 0L) {
    stopifnot(all(sim$meta$low_y_rows$y_raw < 1.0))
  }

  cat("\n[prepare_medsat_london_test] all assertions passed.\n")
  invisible(sim)
}
