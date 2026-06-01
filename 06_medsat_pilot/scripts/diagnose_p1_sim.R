## ============================================================
## scripts/diagnose_p1_sim.R
##
## Read-only diagnostic pass on
##   output/sim_medsat_london_asthma_n500_seed1.rds
## plus a quick re-read of the master + IMD to attribute the 335
## NA-drops to specific columns and boroughs.
##
## Does NOT modify the sim object.
## ============================================================

SIM_RDS     <- "output/sim_medsat_london_asthma_n500_seed1.rds"
MASTER_PATH <- "data/medsat/2019_spatial_raw_master.csv"
IMD_PATH    <- "data/medsat/imd_2019_lsoa21.csv"

suppressPackageStartupMessages({
  stopifnot(requireNamespace("data.table", quietly = TRUE))
})

sim <- readRDS(SIM_RDS)
stopifnot(sim$n == 500L, sim$p == 3L)

cat("==============================================================\n")
cat(" Diagnostic pass on sim_medsat_london_asthma_n500_seed1.rds\n")
cat("==============================================================\n")

# Convenience: build a per-row table for the n=500 subsample with raw values
# back out of the scaled ones. X_scale carries lo/hi per covariate.
X1_lo <- sim$X_scale[[1]]$lo; X1_hi <- sim$X_scale[[1]]$hi
X2_lo <- sim$X_scale[[2]]$lo; X2_hi <- sim$X_scale[[2]]$hi
X3_lo <- sim$X_scale[[3]]$lo; X3_hi <- sim$X_scale[[3]]$hi

sub <- data.table::data.table(
  borough     = sim$borough,
  y_raw       = sim$y_raw,
  NO2_raw     = sim$data$X1 * (X1_hi - X1_lo) + X1_lo,
  NDVI_raw    = sim$data$X2 * (X2_hi - X2_lo) + X2_lo,
  IMD_raw     = sim$data$X3 * (X3_hi - X3_lo) + X3_lo,
  centroid_x  = sim$coords_raw[, 1],
  centroid_y  = sim$coords_raw[, 2]
)

# To attach LSOA21NM we re-read the master file's name+coord columns and
# join on (centroid_x, centroid_y), which uniquely identify LSOAs.
master_keys <- data.table::fread(
  MASTER_PATH,
  select = c("geography code", "LSOA21NM", "centroid_x", "centroid_y",
             "c_total population"),
  check.names = FALSE
)
data.table::setnames(master_keys, "geography code", "LSOA21CD")
data.table::setnames(master_keys, "c_total population", "population")
data.table::setkey(master_keys, centroid_x, centroid_y)

sub <- merge(sub, master_keys, by = c("centroid_x", "centroid_y"),
             all.x = TRUE, sort = FALSE)
stopifnot(nrow(sub) == sim$n)

# ------------------------------------------------------------
# 1) Rows with y_raw < 1.0
# ------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat(" 1) Subsample rows with y_raw < 1.0\n")
cat("--------------------------------------------------------------\n")
low <- sub[y_raw < 1.0][order(y_raw)]
cat(sprintf("Found %d row(s) with y_raw < 1.0 in the n=500 subsample.\n",
            nrow(low)))
old_w <- getOption("width"); options(width = 200)
print(low[, .(LSOA21CD, LSOA21NM, borough, y_raw, NO2_raw, NDVI_raw,
              IMD_raw, centroid_x, centroid_y, population)])
options(width = old_w)
cat(sprintf("\nMinimum y_raw in subsample: %.4f\n", min(sub$y_raw)))
cat(sprintf("median y_raw in subsample : %.4f\n", median(sub$y_raw)))

# ------------------------------------------------------------
# 2) NDVI = 1.000 / NDVI > 0.95
# ------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat(" 2) High-NDVI rows in subsample and in pre-subsample London\n")
cat("--------------------------------------------------------------\n")
hi_ndvi_sub <- sub[NDVI_raw >= 0.95][order(-NDVI_raw)]
cat(sprintf("Subsample rows with e_ndvi >= 0.95: %d\n", nrow(hi_ndvi_sub)))
cat(sprintf("Subsample rows with e_ndvi == 1.0 (exact): %d\n",
            sum(sub$NDVI_raw == 1.0)))
if (nrow(hi_ndvi_sub) > 0L) {
  print(hi_ndvi_sub[, .(LSOA21CD, LSOA21NM, borough,
                        NDVI_raw, NO2_raw, IMD_raw, y_raw)])
}

# Replicate prepare_medsat_london's London-pre-subsample stage for comparison.
master <- data.table::fread(MASTER_PATH, check.names = FALSE)
data.table::setnames(master, "geography code", "LSOA21CD")
master[, borough_parsed := sub("\\s\\d{3}[A-Z]$", "", LSOA21NM)]
LONDON_33 <- c(
  "Barking and Dagenham","Barnet","Bexley","Brent","Bromley",
  "Camden","City of London","Croydon","Ealing","Enfield",
  "Greenwich","Hackney","Hammersmith and Fulham","Haringey",
  "Harrow","Havering","Hillingdon","Hounslow","Islington",
  "Kensington and Chelsea","Kingston upon Thames","Lambeth",
  "Lewisham","Merton","Newham","Redbridge",
  "Richmond upon Thames","Southwark","Sutton","Tower Hamlets",
  "Waltham Forest","Wandsworth","Westminster"
)
master_l <- master[borough_parsed %in% LONDON_33]
imd <- data.table::fread(IMD_PATH)
master_l_imd <- merge(master_l, imd, by = "LSOA21CD", all.x = FALSE, sort = FALSE)
cat(sprintf("\nLondon LSOAs after IMD inner-join: %d\n", nrow(master_l_imd)))

# After-NA pre-subsample pool used by the prep:
required <- c("o_asthma_quantity_per_capita","e_NO2","e_ndvi",
              "imd_score_2019","centroid_x","centroid_y")
keep <- stats::complete.cases(master_l_imd[, ..required])
pool <- master_l_imd[keep]
cat(sprintf("London pre-subsample pool (after NA-drop): %d\n", nrow(pool)))

hi_ndvi_pool <- pool[e_ndvi >= 0.95]
cat(sprintf("Pre-subsample London-pool rows with e_ndvi >= 0.95: %d\n",
            nrow(hi_ndvi_pool)))
cat(sprintf("Pre-subsample London-pool rows with e_ndvi == 1.0 (exact): %d\n",
            sum(pool$e_ndvi == 1.0)))

# ------------------------------------------------------------
# 3) Attribute the 335 NA drops
# ------------------------------------------------------------
cat("\n--------------------------------------------------------------\n")
cat(" 3) NA-drops from London-after-IMD (n=4994) -> pool (n=4659)\n")
cat("--------------------------------------------------------------\n")
cat(sprintf("Total London-after-IMD rows: %d ; pool after NA-drop: %d ; dropped: %d\n",
            nrow(master_l_imd), nrow(pool),
            nrow(master_l_imd) - nrow(pool)))

cat("\nNA counts by column (over London-after-IMD, n=4994):\n")
na_tbl <- data.table::data.table(
  column = required,
  n_NA   = vapply(required, function(c) sum(is.na(master_l_imd[[c]])),
                  integer(1))
)
print(na_tbl[order(-n_NA)])

# Combinations of NA columns among the dropped rows
dropped <- master_l_imd[!keep]
dropped[, pattern := apply(
  is.na(.SD), 1L,
  function(z) paste(required[z], collapse = "+")
), .SDcols = required]
cat("\nCombined NA patterns (only across the 335 dropped rows):\n")
print(dropped[, .N, by = pattern][order(-N)])

# Boroughs: NA-drops per borough vs total London-after-IMD rows per borough
totals  <- master_l_imd[, .(n_total   = .N), by = borough_parsed]
drops_b <- dropped     [, .(n_dropped = .N), by = borough_parsed]
brog <- merge(totals, drops_b, by = "borough_parsed", all.x = TRUE)
brog[is.na(n_dropped), n_dropped := 0L]
brog[, pct_dropped := round(100 * n_dropped / n_total, 1)]
cat("\nNA-drops by borough (sorted by drop rate):\n")
print(brog[order(-pct_dropped)])

cat("\nDone (no modifications).\n")
