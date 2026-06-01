## ============================================================
## scripts/build_imd_lsoa21.R
##
## Build an LSOA21-level IMD 2019 table by reconciling:
##   - IMD_2019_File7.csv          (LSOA11-level IMD 2019, England)
##   - LSOA11_to_LSOA21_lookup.csv (ONS 2011 -> 2021 LSOA mapping)
##
## CHGIND handling:
##   U  Unchanged       : take the single LSOA11 IMD score directly.
##   S  Split           : each LSOA21 child inherits the parent LSOA11 score
##                        (1 row per child in the lookup, so this is just a
##                        copy after the merge).
##   M  Merged          : multiple LSOA11 parents fold into one LSOA21 child;
##                        we take a simple MEAN over parents.
##   X  Irregular       : dropped, with row counts reported.
##
## Output:
##   data/medsat/imd_2019_lsoa21.csv   (LSOA21CD, imd_score_2019,
##                                      imd_source_chgind)
## ============================================================

suppressPackageStartupMessages({
  stopifnot(requireNamespace("data.table", quietly = TRUE))
})
DT <- data.table::data.table   # silence R CMD style

DATA_DIR    <- "data/medsat"
imd_path    <- file.path(DATA_DIR, "IMD_2019_File7.csv")
lookup_path <- file.path(DATA_DIR, "LSOA11_to_LSOA21_lookup.csv")
out_path    <- file.path(DATA_DIR, "imd_2019_lsoa21.csv")

stopifnot(file.exists(imd_path), file.exists(lookup_path))

# ---- 33 London boroughs (ONS canonical, copied from
# R/prepare_medsat_london.R .london_boroughs() to keep this script
# standalone-runnable). ---------------------------------------------------
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
stopifnot(length(LONDON_33) == 33L)

# ---- 1. Read IMD File 7 ----------------------------------------------
imd_raw <- data.table::fread(imd_path, check.names = FALSE)
imd <- imd_raw[, .(
  LSOA11CD  = `LSOA code (2011)`,
  IMD_Score = `Index of Multiple Deprivation (IMD) Score`,
  IMD_Rank  = `Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)`
)]
cat(sprintf("[1] IMD File 7: %d LSOA11 rows.\n", nrow(imd)))

# ---- 2. Read lookup --------------------------------------------------
lu <- data.table::fread(lookup_path, check.names = FALSE)
cat(sprintf("[2] Lookup: %d rows.\n", nrow(lu)))
cat("    Row-level CHGIND counts in lookup:\n")
print(table(lu$CHGIND))

# ---- 3. Left-join IMD onto lookup ------------------------------------
mer <- merge(lu, imd, by = "LSOA11CD", all.x = TRUE, sort = FALSE)
n_na  <- sum(is.na(mer$IMD_Score))
n_welsh_l11 <- sum(startsWith(mer$LSOA11CD, "W"))
cat(sprintf("[3] After left-join: %d lookup rows missing IMD; of those, %d have a Welsh (W01) LSOA11 code.\n",
            n_na, sum(is.na(mer$IMD_Score) & startsWith(mer$LSOA11CD, "W"))))

# ---- 4. Drop X (irregular) and Welsh / NA-IMD rows -------------------
n_X    <- sum(mer$CHGIND == "X")
n_na_e <- sum(is.na(mer$IMD_Score) & !startsWith(mer$LSOA11CD, "W"))
mer <- mer[CHGIND != "X" & !is.na(IMD_Score)]
cat(sprintf("[4] Dropped %d X-rows and %d Welsh-only rows missing IMD; %d non-X English rows still missing IMD.\n",
            n_X, n_welsh_l11, n_na_e))
cat(sprintf("    %d rows remain for aggregation.\n", nrow(mer)))

# ---- 5. Aggregate to LSOA21 level ------------------------------------
imd_l21 <- mer[, .(
  imd_score_2019    = mean(IMD_Score),
  imd_source_chgind = {
    u <- unique(CHGIND)
    if (length(u) == 1L) u else paste(sort(u), collapse = ";")
  },
  n_parents = .N
), by = LSOA21CD]
cat(sprintf("[5] After aggregation: %d unique LSOA21s.\n", nrow(imd_l21)))

# ---- 6. Dedup verification ------------------------------------------
dup <- imd_l21[, .N, by = LSOA21CD][N > 1L]
cat(sprintf("[6] LSOA21CD multiplicity after aggregation: max %d (any >1 is a bug).\n",
            max(imd_l21[, .N, by = LSOA21CD]$N)))
stopifnot(nrow(dup) == 0L)

# ---- 7. CHGIND counts at LSOA21 level: national + London-only -------
# Attach LAD22NM (unique per LSOA21 - verified below) for London filter.
lad_lookup <- unique(lu[, .(LSOA21CD, LAD22NM)])
n_lad_dup <- lad_lookup[, .N, by = LSOA21CD][N > 1L]
stopifnot(nrow(n_lad_dup) == 0L)
imd_l21_lad <- merge(imd_l21, lad_lookup, by = "LSOA21CD", all.x = TRUE,
                     sort = FALSE)

cat("\n[7a] LSOA21-level CHGIND counts (national):\n")
print(imd_l21_lad[, .N, by = imd_source_chgind][order(imd_source_chgind)])

london_idx <- imd_l21_lad$LAD22NM %in% LONDON_33
cat(sprintf("\n[7b] London LSOA21s with IMD: %d (of 33 boroughs).\n",
            sum(london_idx)))
cat("     LSOA21-level CHGIND counts (London-only):\n")
print(imd_l21_lad[london_idx, .N, by = imd_source_chgind][order(imd_source_chgind)])

# ---- 8. Score summary -----------------------------------------------
cat("\n[8] imd_score_2019 summary (national):\n")
print(summary(imd_l21$imd_score_2019))
cat(sprintf("    n missing: %d ; sd: %.3f\n",
            sum(is.na(imd_l21$imd_score_2019)),
            sd(imd_l21$imd_score_2019, na.rm = TRUE)))

# ---- 9. Save --------------------------------------------------------
out <- imd_l21[, .(LSOA21CD, imd_score_2019, imd_source_chgind)]
data.table::setkey(out, LSOA21CD)
data.table::fwrite(out, out_path)
cat(sprintf("\n[9] Wrote %s (%d rows, %d cols).\n",
            out_path, nrow(out), ncol(out)))
