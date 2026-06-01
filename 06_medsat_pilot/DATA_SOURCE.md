# MEDSAT data source

The raw MEDSAT input files under `06_medsat_pilot/data/medsat/` are **not tracked
in git** (≈297 MB; excluded via `.gitignore`). This file documents where they came
from so the pilot is reconstructable.

## Source

- **Repository:** TUM mediaTUM, `https://mediatum.ub.tum.de/1714817`
- **Dataset:** MEDSAT — satellite + sociodemographic dataset of LSOA-level
  prescription rates in England.
- **Citation:** Scepanovic et al., *MedSat: A Public Health Dataset for England
  Featuring Medical Prescriptions and Satellite Imagery*, NeurIPS 2023 Datasets &
  Benchmarks Track. OpenReview ID `CSJYz1Zovj`.
- **License:** Open.
- **Granularity:** Lower Layer Super Output Areas (LSOAs) in England (~33,000
  areas), years 2019 and 2020.
- **Per LSOA:** five prescription-rate outcomes (diabetes, hypertension, asthma,
  depression, anxiety); 43 environmental features (NO2, NDVI, PM, temperature,
  humidity, LST); 111 sociodemographic features (age, ethnicity, IMD, income,
  education); four seasonal Sentinel-2 composite tiles.

## Files expected in `06_medsat_pilot/data/medsat/`

| File | Size | Description |
|------|------|-------------|
| `2019_spatial_raw_master.csv` | ~149 MB | Master table, all LSOAs, year 2019 (pilot uses this). |
| `2020_spatial_raw_master.csv` | ~149 MB | Master table, all LSOAs, year 2020. |
| `IMD_2019_File7.csv` | ~10 MB | English Indices of Multiple Deprivation 2019 (File 7). |
| `LSOA11_to_LSOA21_lookup.csv` | ~3 MB | ONS LSOA 2011→2021 code lookup. |
| `imd_2019_lsoa21.csv` | ~0.6 MB | IMD 2019 re-keyed to LSOA 2021 codes. |
| `MedSat Variables.xlsx` | ~1.8 MB | Variable dictionary / column descriptions. |

The IMD and LSOA-lookup files are standard ONS / MHCLG open releases (English
Indices of Deprivation 2019; ONS geography lookups) and can be re-downloaded from
the UK government open-data portals if the mediaTUM bundle does not include them.

## Pilot subset (derived, not raw)

The pilot uses Greater London (~4,800 LSOAs), year 2019, asthma prescription rate
as outcome, with covariates NO2, NDVI, IMD (p=3) and headline interaction NO2×NDVI,
subsampled to n=500. The derivation is performed in code (`prepare_medsat_london.R`)
from the raw files above; see `06_medsat_pilot/CLAUDE.md` for the full pilot spec.
