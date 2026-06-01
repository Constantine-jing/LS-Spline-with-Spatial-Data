# ============================================================
# save_munich_data.R
# Run this ONCE on your local machine (where gamlss.data is installed)
# to export rent99 and rent99.polys as standalone .rda files.
#
# Then upload both .rda files to Hellbender alongside run_munich_rent.R
#
# Usage: Rscript save_munich_data.R
# ============================================================

library(gamlss.data)

data(rent99)
data(rent99.polys)

cat("rent99:      ", nrow(rent99), "obs x", ncol(rent99), "cols\n")
cat("rent99.polys:", length(rent99.polys), "districts\n")

save(rent99, file = "rent99.rda")
save(rent99.polys, file = "rent99_polys.rda")

cat("Saved: rent99.rda, rent99_polys.rda\n")
cat("Upload both to Hellbender in the same directory as run_munich_rent.R\n")
