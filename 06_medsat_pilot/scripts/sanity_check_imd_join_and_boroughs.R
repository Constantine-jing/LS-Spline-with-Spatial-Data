## ============================================================
## scripts/sanity_check_imd_join_and_boroughs.R
##
## Step 4 -- inner-join data/medsat/imd_2019_lsoa21.csv onto the MEDSAT
##           2019 master file (`geography code` == LSOA21CD). Report
##           coverage nationally and London-only. Show first 5 rows of
##           the joined London subset. The merged file is NOT saved.
##
## Step 5 -- Borough parsing on the 4,994 London LSOAs:
##              regex:  sub("\\s\\d{3}[A-Z]$", "", LSOA21NM)
##           Verify: 0 parse failures, exactly 33 unique boroughs, all
##           match .london_boroughs() in R/prepare_medsat_london.R.
## ============================================================

suppressPackageStartupMessages({
  stopifnot(requireNamespace("data.table", quietly = TRUE))
})

master_path <- "data/medsat/2019_spatial_raw_master.csv"
imd_path    <- "data/medsat/imd_2019_lsoa21.csv"
helper_path <- "R/prepare_medsat_london.R"

stopifnot(file.exists(master_path),
          file.exists(imd_path),
          file.exists(helper_path))

# ---- Load only the columns we need from MEDSAT master ---------------
need <- c("geography code","LSOA21NM",
          "o_asthma_quantity_per_capita","e_NO2","e_ndvi",
          "centroid_x","centroid_y")
master <- data.table::fread(master_path, select = need, check.names = FALSE)
data.table::setnames(master, "geography code", "LSOA21CD")
cat(sprintf("MEDSAT master 2019: %d LSOAs.\n", nrow(master)))

imd <- data.table::fread(imd_path)
cat(sprintf("imd_2019_lsoa21:    %d LSOAs.\n\n", nrow(imd)))

# ---- Load the canonical 33-borough list from the helper.
# The helper file sources spatial_utils.R via a relative path, so
# sys.source with chdir=TRUE so that relative source() inside resolves.
helper_env <- new.env()
sys.source(helper_path, chdir = TRUE, envir = helper_env)
london_boroughs <- helper_env$.london_boroughs()
stopifnot(length(london_boroughs) == 33L)

# ---- Step 4: inner-join ---------------------------------------------
joined <- merge(master, imd, by = "LSOA21CD", all.x = FALSE, all.y = FALSE,
                sort = FALSE)
cat(sprintf("[Step 4] Inner-join MEDSAT x IMD: %d LSOAs.\n", nrow(joined)))

missing <- master[!LSOA21CD %in% imd$LSOA21CD]
cat(sprintf("MEDSAT LSOAs WITHOUT matching IMD: %d\n", nrow(missing)))
if (nrow(missing) > 0L) {
  cat("  All such codes (LSOA21CD, LSOA21NM):\n")
  print(missing[, .(LSOA21CD, LSOA21NM)])
  welsh <- missing[startsWith(LSOA21CD, "W")]
  cat(sprintf("  ...of which %d are Welsh (W01 prefix).\n", nrow(welsh)))
}

# London-only views
master[, borough := sub("\\s\\d{3}[A-Z]$", "", LSOA21NM)]
master_l <- master[borough %in% london_boroughs]
joined[,  borough := sub("\\s\\d{3}[A-Z]$", "", LSOA21NM)]
joined_l <- joined[borough %in% london_boroughs]
miss_l   <- master_l[!LSOA21CD %in% imd$LSOA21CD]
cat(sprintf("\nLondon-only: MEDSAT has %d LSOAs, joined %d, missing %d.\n",
            nrow(master_l), nrow(joined_l), nrow(miss_l)))
if (nrow(miss_l) > 0L) {
  cat("  London LSOAs missing IMD:\n")
  print(miss_l[, .(LSOA21CD, LSOA21NM, borough)])
}

# First 5 rows of the joined London subset
show_cols <- c("LSOA21CD","LSOA21NM","o_asthma_quantity_per_capita",
               "e_NO2","e_ndvi","imd_score_2019","imd_source_chgind",
               "centroid_x","centroid_y")
cat("\n[Step 4] First 5 rows of the joined London subset:\n")
# print wide so columns aren't truncated
old_width <- getOption("width"); options(width = 200)
print(head(joined_l[, ..show_cols], 5))
options(width = old_width)

# ---- Step 5: borough parsing on London LSOAs ------------------------
PARSE_REGEX <- "\\s\\d{3}[A-Z]$"
cat(sprintf("\n[Step 5] Regex used to parse borough from LSOA21NM:\n"))
cat(sprintf('   sub("%s", "", LSOA21NM)\n', PARSE_REGEX))
cat("   strips a trailing space + 3 digits + 1 capital letter,\n")
cat("   e.g. 'Barking and Dagenham 016A' -> 'Barking and Dagenham'.\n\n")

parsed <- sub(PARSE_REGEX, "", joined_l$LSOA21NM)
n_fail <- sum(parsed == joined_l$LSOA21NM)
ub <- sort(unique(parsed))
extras  <- setdiff(ub, london_boroughs)
missing_b <- setdiff(london_boroughs, ub)

cat(sprintf("Applied to the %d London LSOAs:\n", nrow(joined_l)))
cat(sprintf("  parse failures : %d  (expect 0)\n", n_fail))
cat(sprintf("  unique boroughs: %d  (expect 33)\n", length(ub)))
cat(sprintf("  identical to .london_boroughs()? %s\n",
            identical(ub, sort(london_boroughs))))
if (length(extras)    > 0L) cat("  EXTRA names:   ", paste(extras,    collapse=", "), "\n")
if (length(missing_b) > 0L) cat("  MISSING names: ", paste(missing_b, collapse=", "), "\n")

cat("\nDone. The merged file was NOT saved.\n")
