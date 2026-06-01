## Read-only inspection of MEDSAT files (2019 only).
## Run from project root:  Rscript scripts/inspect_medsat_files.R

suppressPackageStartupMessages({
  ok_readxl <- requireNamespace("readxl", quietly = TRUE)
  ok_dt     <- requireNamespace("data.table", quietly = TRUE)
})
stopifnot(ok_readxl, ok_dt)

DATA_DIR <- "data/medsat"
codebook_path <- file.path(DATA_DIR, "MedSat Variables.xlsx")
master_path   <- file.path(DATA_DIR, "2019_spatial_raw_master.csv")

stopifnot(file.exists(codebook_path), file.exists(master_path))

cat("======================================================================\n")
cat("PART 1: CODEBOOK  (", codebook_path, ")\n", sep = "")
cat("======================================================================\n")
sheets <- readxl::excel_sheets(codebook_path)
cat("Sheets:\n"); print(sheets); cat("\n")

# Read every sheet and stash it
codebook <- lapply(sheets, function(s) {
  readxl::read_excel(codebook_path, sheet = s, .name_repair = "minimal")
})
names(codebook) <- sheets

for (s in sheets) {
  cat("--- Sheet:", s, " | dim =", paste(dim(codebook[[s]]), collapse = " x "), "---\n")
  cat("    cols:", paste(names(codebook[[s]]), collapse = " | "), "\n")
}
cat("\n")

# Concatenate all sheets into one searchable text frame: (sheet, col, row_text)
flatten_sheet <- function(s) {
  df <- codebook[[s]]
  if (nrow(df) == 0) return(NULL)
  txt <- apply(df, 1, function(r) paste(r, collapse = " || "))
  data.frame(sheet = s, row = seq_along(txt), text = txt, stringsAsFactors = FALSE)
}
flat <- do.call(rbind, lapply(sheets, flatten_sheet))

search_codebook <- function(pattern, label) {
  cat("### Search:", label, "  (regex: ", pattern, ")\n", sep = "")
  hits <- flat[grepl(pattern, flat$text, ignore.case = TRUE), ]
  if (nrow(hits) == 0) {
    cat("    (no matches)\n\n")
  } else {
    # Limit so output is readable
    n_show <- min(nrow(hits), 12)
    for (i in seq_len(n_show)) {
      cat(sprintf("    [%s row %d] %s\n",
                  hits$sheet[i], hits$row[i], substr(hits$text[i], 1, 300)))
    }
    if (nrow(hits) > n_show) cat(sprintf("    ... %d more matches\n", nrow(hits) - n_show))
    cat("\n")
  }
}

search_codebook("asthma",                            "asthma prescription rate")
search_codebook("(^|[^a-z])no2|nitrogen dioxide",    "NO2")
search_codebook("ndvi|vegetation index",             "NDVI")
search_codebook("imd|multiple deprivation",          "IMD")
search_codebook("lsoa",                              "LSOA identifier")
search_codebook("easting|northing|longitude|latitude|lon|lat|centroid|coord",
                "coordinates")
search_codebook("borough|local authority|ladcd|lad21|district|lad ",
                "borough / LAD")

cat("======================================================================\n")
cat("PART 2: MASTER FILE  (", master_path, ")\n", sep = "")
cat("======================================================================\n")
fsize <- file.info(master_path)$size
cat(sprintf("File size: %.1f MB\n\n", fsize / 1024^2))

# Read header only, get every column name
hdr <- data.table::fread(master_path, nrows = 0)
all_cols <- names(hdr)
cat("Total columns:", length(all_cols), "\n\n")

# Helper: find columns by regex, report hits
find_cols <- function(pattern, label) {
  hits <- grep(pattern, all_cols, ignore.case = TRUE, value = TRUE)
  cat(sprintf("[%s]  (regex: %s)\n", label, pattern))
  if (length(hits) == 0) {
    cat("    (none)\n\n")
  } else {
    cat("   ", paste(hits, collapse = " | "), "\n\n")
  }
  hits
}

c_asthma  <- find_cols("asthma",                                    "asthma")
c_no2     <- find_cols("(^|[^a-z])no2|nitrogen",                    "NO2")
c_ndvi    <- find_cols("ndvi",                                      "NDVI")
c_imd     <- find_cols("imd|deprivation",                           "IMD")
c_lsoa    <- find_cols("^lsoa|lsoa[0-9]+(cd|nm)|^code$|^geography$", "LSOA id")
c_coords  <- find_cols("east|north|lon|lat|centroid|x_|y_|coord|^x$|^y$",
                       "coordinates")
c_lad     <- find_cols("borough|^lad|ladcd|laduanm|lad21|local_auth|district",
                       "borough / LAD")

# Row count via wc -l minus 1 if header.
# Use data.table::fread with select=1 — faster than reading full file but still slow for 149MB.
# Cheaper: count newlines.
cat("Counting rows (via R.utils::countLines if available, else fallback)...\n")
nrow_est <- tryCatch({
  if (requireNamespace("R.utils", quietly = TRUE)) {
    as.integer(R.utils::countLines(master_path)) - 1L
  } else {
    # Streaming fallback: read in chunks of 1e5 lines
    con <- file(master_path, "r")
    on.exit(close(con))
    total <- 0L
    repeat {
      chunk <- readLines(con, n = 200000L)
      if (length(chunk) == 0) break
      total <- total + length(chunk)
    }
    total - 1L  # subtract header
  }
}, error = function(e) NA_integer_)
cat(sprintf("Data rows (excluding header): %s\n\n", format(nrow_est, big.mark = ",")))

# Read first 5 rows for the columns we care about
want <- unique(c(c_asthma[1], c_no2[1], c_ndvi[1], c_imd[1],
                 c_lsoa[1], c_coords[1:2], c_lad[1]))
want <- want[!is.na(want)]
if (length(want) > 0) {
  head_df <- data.table::fread(master_path, nrows = 5L, select = want)
  cat("First 5 rows of selected columns:\n")
  print(head_df)
  cat("\n")
}

cat("======================================================================\n")
cat("END\n")
