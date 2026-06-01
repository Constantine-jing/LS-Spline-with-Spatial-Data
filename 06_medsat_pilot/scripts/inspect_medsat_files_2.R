## Follow-up inspection: dump ALL columns, look harder for LSOA code + IMD,
## verify borough parseability from LSOA21NM.

suppressPackageStartupMessages({
  stopifnot(requireNamespace("data.table", quietly = TRUE),
            requireNamespace("readxl",     quietly = TRUE))
})

master_path   <- "data/medsat/2019_spatial_raw_master.csv"
codebook_path <- "data/medsat/MedSat Variables.xlsx"

# --- Full column dump (header only, very fast) -----------------------------
hdr <- data.table::fread(master_path, nrows = 0)
all_cols <- names(hdr)
cat("Total columns:", length(all_cols), "\n\n")

# Bucket by prefix
prefix <- ifelse(grepl("^[a-z]_", all_cols), substr(all_cols, 1, 2), "other")
cat("Column prefix counts:\n")
print(table(prefix))
cat("\n")

# Dump non-prefixed columns (likely IDs/coords)
cat("Non-prefixed columns (likely IDs / geometry / outcome):\n")
print(all_cols[!grepl("^[a-z]_", all_cols)])
cat("\n")

# Sociodemographic c_ columns with 'depriv', 'imd', 'income', 'rank', 'decile'
cat("c_* columns that *might* be IMD / deprivation related:\n")
print(grep("depriv|imd|income|rank|decile|class|edu",
           all_cols, ignore.case = TRUE, value = TRUE))
cat("\n")

# Anything resembling an LSOA code (E01xxxxxxx)
cat("Read first row, scan for E01-prefixed values:\n")
row1 <- data.table::fread(master_path, nrows = 1L)
e01_cols <- names(row1)[sapply(row1, function(v) {
  is.character(v) && length(v) >= 1 && grepl("^E0[12]\\d{7}$", v[1])
})]
print(e01_cols)
cat("\n")

# --- Codebook auxiliary sheet, full content --------------------------------
aux <- readxl::read_excel(codebook_path, sheet = "auxiliary", .name_repair = "minimal")
cat("Auxiliary sheet (full):\n")
print(aux)
cat("\n")

# --- Codebook outcomes sheet, full content ---------------------------------
out <- readxl::read_excel(codebook_path, sheet = "outcomes", .name_repair = "minimal")
cat("Outcomes sheet (full):\n")
print(out)
cat("\n")

# --- Borough parseability check from LSOA21NM ------------------------------
nm_sample <- data.table::fread(master_path, select = "LSOA21NM", nrows = 50000L)
nms <- nm_sample$LSOA21NM
# LSOA21NM convention: "<Local Authority District Name> NNNL" e.g. "Camden 001A"
borough <- sub("\\s\\d{3}[A-Z]$", "", nms)
n_parse_fail <- sum(borough == nms)
cat(sprintf("LSOA21NM parsed borough; %d / %d names failed to strip the trailing 'NNNL' suffix.\n",
            n_parse_fail, length(nms)))
cat("Unique borough names (first 50):\n")
ub <- sort(unique(borough))
print(head(ub, 50))
cat(sprintf("Total unique boroughs/LADs in 50k-row sample: %d\n\n", length(ub)))

# Check whether the 33 ONS London boroughs are present
london33 <- c("City of London","Barking and Dagenham","Barnet","Bexley","Brent",
              "Bromley","Camden","Croydon","Ealing","Enfield","Greenwich","Hackney",
              "Hammersmith and Fulham","Haringey","Harrow","Havering","Hillingdon",
              "Hounslow","Islington","Kensington and Chelsea","Kingston upon Thames",
              "Lambeth","Lewisham","Merton","Newham","Redbridge","Richmond upon Thames",
              "Southwark","Sutton","Tower Hamlets","Waltham Forest","Wandsworth",
              "Westminster")
cat("London boroughs found in sample:\n")
found <- london33[london33 %in% ub]
missing <- london33[!london33 %in% ub]
cat(sprintf("  found: %d / 33\n", length(found)))
if (length(missing) > 0) cat("  missing:", paste(missing, collapse=", "), "\n")

# Within London-only, count LSOAs
cat("\nLSOA row count by 'is in London-33' (50k sample):\n")
in_london <- borough %in% london33
print(table(in_london))
cat("\n")
cat("Likely full-file London LSOA count (~33 * mean per borough):\n")
# Get exact count by reading the full LSOA21NM column
full_nm <- data.table::fread(master_path, select = "LSOA21NM")
full_borough <- sub("\\s\\d{3}[A-Z]$", "", full_nm$LSOA21NM)
cat(sprintf("  Full-file rows: %d\n", nrow(full_nm)))
cat(sprintf("  Rows whose parsed borough is in the 33 London LADs: %d\n",
            sum(full_borough %in% london33)))

# --- Quick sanity peek at the IMD-percent-deprivation columns --------------
imd_like <- grep("depriv", all_cols, ignore.case = TRUE, value = TRUE)
cat("\nColumns matching 'depriv':\n")
print(imd_like)
if (length(imd_like) > 0) {
  peek <- data.table::fread(master_path, select = imd_like, nrows = 5L)
  cat("First 5 rows:\n")
  print(peek)
}
