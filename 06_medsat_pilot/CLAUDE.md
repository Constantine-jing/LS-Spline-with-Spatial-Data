# CLAUDE.md — medsat-pilot

Project context for Claude Code. Read this first in any new session.

---

## What this project is

This is the **real-data application** for the user's Statistics PhD Project 1: a
Bayesian LS-spline spatial additive model with pairwise tensor-product
interactions, applied to the MEDSAT (Scepanovic et al., 2023) dataset of
LSOA-level prescription rates in England.

This folder is NOT where the methodology lives — that's the larger dissertation
codebase. This folder is the **pilot application** that demonstrates the method
on real public-health data, scoped to a subset of London LSOAs.

The methodology (Bayesian LS-spline + Matérn GP + collapsed Gibbs sampler with
pairwise interactions) is already validated on simulation. Project scope here
is application, not methodology development.

---

## The MEDSAT data

- **Source:** TUM mediaTUM repository at `mediatum.ub.tum.de/1714817`.
- **Citation:** Scepanovic et al., NeurIPS 2023 Datasets & Benchmarks Track,
  OpenReview ID `CSJYz1Zovj`.
- **Granularity:** Lower Layer Super Output Areas (LSOAs) in England,
  ~33,000 areas, years 2019 and 2020.
- **Per LSOA:** five prescription-rate outcomes (diabetes, hypertension,
  asthma, depression, anxiety), 43 environmental features (NO2, NDVI, PM,
  temperature, humidity, LST), 111 sociodemographic features (age, ethnicity,
  IMD, income, education), plus four seasonal Sentinel-2 composite tiles.
- **License:** Open.

**Existing analyses to position against:**
- Scepanovic et al. (2023): RF and gradient-boosting baselines.
- Maitra et al. (arXiv:2501.02111, Jan 2025): XGBoost + interpretable ML.
- **Gap this project fills:** no published analysis uses a Bayesian additive
  spline + spatial GP + pairwise tensor-product interaction model.

---

## Pilot specification (locked)

- **Region:** Greater London (~4,800 LSOAs in the 33 boroughs).
- **Year:** 2019 (pre-COVID).
- **Outcome (Y):** asthma prescription rate, continuous, treated as Gaussian
  after any transform check.
- **Covariates (p=3):** NO2 (environmental, pollution), NDVI (environmental,
  greenness), IMD (sociodemographic, deprivation index).
- **Headline interaction:** NO2 × NDVI. Scientific question: "Does urban
  greenery buffer the effect of pollution on respiratory health?"
- **Sample size for pilot:** n=500 LSOAs, stratified random subsample across
  boroughs.

---

## Why n=500 (subsample rationale)

The full London LSOA set has ~4,800 LSOAs; England has ~33,000.
The sampler's bottleneck is the O(n^3) Cholesky of the n x n Matérn covariance
at every MH proposal of rho. Munich at n=3,082 took ~61 hours.

Strategy is subsample (no code change, cleanest framing) for the pilot, then
decide on full-scale strategy later. At n=500 the pilot runs in roughly the
same regime as the validated 2x2 interaction experiments and Scenario B in
the methodology paper.

If we eventually want full-coverage runs:
- Path (a) [chosen for pilot]: subsample to n ~ 1000-2000.
- Path (b): aggregate to MSOAs (~7,000 nationally, ~1,000 London).
- Path (c): SPDE / Nyström low-rank approximations — separate methodology
  project, explicitly OUT OF SCOPE for Project 1.

---

## Pilot run plan (P0–P4)

| Run     | n    | p | Interactions | Status (2026-05-25)           | Artefact RDS                                              |
|---------|------|---|--------------|-------------------------------|-----------------------------------------------------------|
| P0      | 500  | 3 | mgcv sanity  | **DONE** (~1 min)             | `output/fit_medsat_p0_mgcv_local.rds`                     |
| P1      | 500  | 3 | none (ours additive) | **DONE** (~20 min, M=20) | `output/fit_medsat_p1_additive_n500_seed1_local.rds`      |
| P1-diff | 500  | 3 | none, diffuse rho prior | **DONE** (rho-prior sensitivity check) | `output/fit_medsat_p1_additive_n500_seed1_diffuse_rho_local.rds` |
| P1plus  | 500  | 3 | all-pairs (replaces P2) | **DONE** (~75 min)        | `output/fit_medsat_p1plus_allpairs_n500_seed1_local.rds`  |
| P3      | 1500 | 3 | all-pairs    | **DONE on Hellbender** (7.25 h, 66 MB) | `output/fit_medsat_n1500_seed1_allpairs_hb.rds`           |
| P4      | 2000–3000 | TBD | 2–3 interactions | not started               | —                                                          |

Note: "P2" is no longer a separate run in this folder. The full-interaction
n=500 fit was re-tagged as `p1plus` after the `fit_interactions` flag was
added, since the only thing that distinguishes it from P1 is the flag.

Run P0–P1 IN ORDER. The point is to catch data-prep errors at minutes-scale
before they cost hours of compute. Skipping P0 (the mgcv sanity check) is
how the user lost time on the Munich τ²/σ² imbalance issue — that would
have been caught in P0 at ~1 minute.

**MCMC budget per run is configurable in `settings`:**
- Validated default (Scenario B): `n_iter=10500, n_burn=2500, n_thin=4, n_draws=2000`. Runs ~hours at n=500.
- Pilot-fast (used for P1, P1plus, P3): `n_iter=3000, n_burn=1000, n_thin=1, n_draws=2000`.
- Constraint: `floor((n_iter - n_burn) / n_thin)` must equal `n_draws`
  (hard `stopifnot` in the wrapper).

---

## Key files in `R/`

These are sourced as a dependency chain. **Do not edit anything in R/ during
the pilot** unless the diff is recorded and discussed — these are the
production methodology files. New code goes in `scripts/`.

- `prepare_medsat_london.R` — NEW. Builds a sim-shaped object from a MEDSAT
  CSV/RDS file. Handles London filtering, column-name mapping, NA dropping,
  stratified subsampling by borough, [0,1] scaling of covariates, unit-square
  scaling of coordinates, and optional Y transformation. Returns a list with
  the schema below.
- `fit_ours_interaction_wrapper.R` — PATCHED in May 2026 with "Option A" so it
  tolerates real-data sim objects whose `truth_params$rho`, `$sigma2`,
  `$tau2_s` are NA. Falls back to `rho_init = 0.2`,
  `sigma2_init = tau2_s_init = var(y)/2`. `nu` is still required (it's a
  fixed modelling choice, not a truth quantity).
  Also patched (May 2026) with a `fit_interactions` setting flag — see the
  "fit_interactions flag" section below.
- `ls_basis.R` — Lancaster–Šalkauskas natural cubic spline basis construction
  (Chib & Greenberg 2010 formulation).
- `ls_interaction.R` — tensor-product interaction basis, Khatri–Rao design,
  ANOVA identifiability via T_u ⊗ T_v, 2D Kronecker-sum RW2 penalty. The
  `orthogonalize=TRUE` flag enforces design-level orthogonality of W_uv
  against `[1, W_u, W_v]` (the Lang–Brezger centering fix).
- `ls_interaction_core.cpp` — Rcpp/RcppEigen hot-loop kernels:
  `khatri_rao_cpp`, `build_block_prior_precision_cpp`, `whiten_cpp`,
  `sample_mvn_chol_cpp`, `quad_form_cpp`.
- `gibbs_stage_c_full.R` — main-effects-only collapsed Gibbs sampler.
- `gibbs_interaction.R` — extended sampler with interaction Step 5b. Also has
  diagnostic clamp hooks `tau2_s_main_fixed`, `tau2_s_int_fixed`,
  `sigma2_fixed`, `tau2_fixed`, `rho_fixed` used in the May 2026 undersmoothing
  diagnosis. Step 5b is gated on `fit_interactions`; when FALSE, no
  interaction blocks are drawn even if `int_keys` is non-empty.
- `spatial_utils.R` — `pairdist`, `matern_cor`.
- `fit_spatial_reml.R` — REML init used for MCMC starting values.
- `canonical_schema.R` — schema definition + `validate_fit_result` for the
  canonical fit-result object that downstream metric scripts consume.

---

## The "sim" object schema (real data version, Option A)

`prepare_medsat_london()` returns a list with these fields:

**Real fields (have values):**
- `scenario` = `"medsat_london"`
- `seed`, `n`, `p`
- `data` — data.frame with columns `X1`..`Xp` (scaled to [0,1]),
  `lon`, `lat` (coords scaled to unit square jointly), `y`
- `x_grid_1d`, `x_grid_2d`, `flat_grid_2d` — evaluation grids
- `int_keys` — character vector like `c("1_2", "1_3", "2_3")` for the
  pairs the wrapper will fit (always all-pairs at present)
- `truth_params$nu` — fixed modelling choice (default `1.0`)
- `coords_raw`, `coords_scale` — for back-transforming plots
- `X_scale` — list of `(name, lo, hi)` per covariate, for back-transforming
  fitted f_j curves to original axis units (e.g. NO2 µg/m³)
- `borough` — borough label per row, for stratification diagnostics
- `y_raw` — untransformed response for back-transforming fitted curves
- `meta`, `settings` — provenance and reproducibility

**Real-data NULL/NA fields (no ground truth):**
- `truth_f_grid = NULL`, `truth_s_obs = NULL`,
  `truth_f_int = NULL`, `truth_f12_obs = NULL`
- `truth_params$sigma2 = NA`, `$tau2 = NA`, `$tau2_s = NA`,
  `$rho = NA`, `$c_int = NA`

These NULLs are intentional: any sim-only metric script will fail loudly
rather than silently treating estimates as truth.

---

## Things baked into the prep that may need user input

1. **Column-name map.** Defaults assume MEDSAT files use plain column names
   (`asthma`, `NO2`, `NDVI`, `IMD`, `easting`, `northing`, `borough`,
   `LSOA21CD`). MEDSAT may use suffixed names like `asthma_2019` or
   `no2_mean_2019`. The user can pass `col_map = list(...)` to override.

2. **Year filtering.** Not implemented. If MEDSAT is long-format (one row per
   LSOA × year), filter to `year == 2019` upstream. If wide-format with
   per-year columns, use `col_map` to point at the right column.

3. **Response transform.** Default `identity`. Check the histogram of
   `y_raw` once data is in hand; switch to `"log"` if right-skewed.

4. **Greater London filter mode.** Default `"borough_list"` matches by
   borough name against the 33 ONS-canonical borough names. If the file
   has LAD codes instead, use `london_filter = "lad_code_prefix"`.
   If pre-filtered upstream, use `"none"`.

5. **`int_pairs`.** Recorded but the wrapper currently fits all p*(p-1)/2
   pairs. With p=3 you get all of `1_2`, `1_3`, `2_3`. Fitting only a
   subset would require a non-trivial change in `ls_interaction.R`,
   out of scope for the pilot.

---

## The `fit_interactions` setting flag (May 2026 addition)

Both `fit_ours_interaction_wrapper.R` and `gibbs_interaction.R` now respect a
`fit_interactions` setting (default `TRUE`):

- `TRUE`  — all-pairs interaction blocks are built and Step 5b of the
  collapsed Gibbs sampler is run. This is the production behaviour and was
  used for P1plus and P3.
- `FALSE` — interaction-basis construction is skipped, the design matrix is
  `H = [1 | W_main]` only, and Step 5b is bypassed. `sim$int_keys` is still
  consulted but no posterior is sampled — `fit$f_int` and
  `fit$var_comp$tau2_int` come back as empty lists. This is what produced
  the additive P1 fit.

Use `fit_interactions = FALSE` only when you want a clean additive baseline
for diagnostics. Production runs (and anything used in the headline)
should have it `TRUE`.

---

## Key empirical findings so far (P0 → P3)

1. **mgcv P0:** NO2 × NDVI interaction is non-significant on London asthma
   in 2019 (p ≈ 0.12, ΔAIC ≈ 0). The spatial smooth `te(lon, lat)` has
   EDF ≈ 30 and k-index ≈ 0.72 — strong residual spatial signal that the
   Matérn GP in our model is the appropriate tool for.

2. **Bayesian P1 vs mgcv (re-centered, centering-invariant metric):**
   curve correlation 0.97 (NO2), 0.95 (NDVI), 0.79 (IMD). Methods agree on
   smooth shapes where there is signal. On IMD mgcv shrinks to flat
   (edf ≈ 1); our additive Bayesian fit is consistent with this.

3. **All-pairs vs additive Bayesian at n=500:** `sigma2`, `tau2_s`, `rho`
   differ by < 6 % between the additive and all-pairs fits. The
   interaction blocks absorb negligible variance.

4. **`rho` is DATA-IDENTIFIED at ≈ 0.017** (unit-square coordinates;
   this corresponds to a Matérn range of about a few hundred metres on
   the London scale). A diffuse rho-prior gives a 1.3 % mean shift and
   91 % interval overlap with the original prior bundle; the prior
   `log(0.2), sd=1.0` is kept. See `notes/methods_phrases.md` for the
   write-up paragraph.

5. **n=1500 Hellbender run (P3):** 7.25 h wall, 66 MB RDS. MH acceptance
   `sigma2=0.10` (low — proposal can be widened on the next run),
   `tau2=0.55` (acceptable after tuning), `rho=0.20` (acceptable).
   - `rho` 95 % CI [0.0143, 0.0171] — even tighter than n=500, consistent.
   - `sigma2` shrunk from ≈ 0.125 (n=500) to ≈ 0.077 (n=1500), as expected.
   - Main-effect `tau2_s[j]` shrunk with more data (good).
   - **WARNING — interaction-variance ESS is poor:** the three
     `tau2_int[k]` chains have effective sample sizes 8.4, 12.0, 20.6 out
     of 2000 draws. The headline scientific question depends on these
     blocks. We should not draw conclusions about NO2 × NDVI from this
     fit alone — either rerun with longer chains / better tuning for
     the interaction blocks, or run multiple seeds and pool.

   Tuning that produced these acceptance rates:
   - `mh_sd_log_rho = 0.15` (down from default 0.3; rho is tight)
   - `mh_sd_log_tau2 = 0.6` (up from default 0.3 to break a sticky 0.62)
   - `mh_sd_log_sigma2` left at default 0.5 (accept still ≈ 0.10)

---

## Workflow: Hellbender division of labor

Anything that runs on the MU "Hellbender" HPC cluster goes through the user,
not Claude Code:

- **Claude Code** writes scripts locally on Windows, prepares `.sub`
  SLURM files in `scripts/submit/`, and runs diagnostics / produces plots
  on fit RDS files once they are returned.
- **The user** handles all Hellbender access manually via VPN: copying
  the sim/RDS up, `sbatch`ing the job, monitoring it, and copying the
  fit RDS back into `output/`.

Claude Code MUST NOT attempt `ssh`, `scp`, `sbatch`, or any other
direct Hellbender access. The `.sub` files are templates the user
launches by hand.

See `scripts/HELLBENDER_README.md` for the current submission recipe and
`scripts/submit/run_n1500_allpairs.sub` for the canonical n=1500 example.

---

## Methodological decisions locked for the pilot

These match the dissertation's Scenario B comparison study:

- `nu = 1.0` (Matérn smoothness)
- `M = 20` (knots per smooth)
- `orthogonalize = TRUE` (the Lang–Brezger design-level centering fix)
- Default prior bundle B1 from the prior-sensitivity analysis (a_smooth=1,
  b_smooth=0.005; a_sigma=2, b_sigma=1; a_tau=2, b_tau=0.3;
  log_rho_mu=log(0.2), log_rho_sd=1.0)

**Known unresolved methodology issue (May 2026):**
The conjugate IG update on τ²_{s,j} is amplitude-sensitive — large-amplitude
smooths can push τ²_{s,j} too high and produce visible undersmoothing on
that smooth. In simulation this affects f_1 in Scenario B. The fix
(half-Cauchy or PC prior on σ_{s,j}) is still being validated in the
methodology folder. For the MEDSAT pilot we proceed with the conjugate IG
prior and will note any visible undersmoothing in the fitted curves.

---

## Working principles (from prior conversations)

- **Sanity checks before full multi-seed runs.** Running ablations and
  diagnostic single-seed experiments before committing to expensive
  multi-seed HPC jobs has consistently caught structural problems early.
  This is why the pilot runs are staged P0→P4 with increasing cost.

- **Preserve old results when testing new configurations.** Use separate
  `output/` subdirectories per run config. Do not overwrite an existing
  fit RDS in-place.

- **New code in new files, not edits to old.** Specifically: new pilot
  runners go in `scripts/`; do not edit anything in `R/` unless the diff
  is recorded.

- **Keep implementation details private** until work is published or
  archived. Do not paste full code in poster sessions or external
  presentations.

- **No psychoanalysis, no overconfident claims on real data.** This is
  the application section; no causal claims, only described associations
  with uncertainty quantification.

---

## File-naming conventions

- Prep outputs: `sim_medsat_london_<outcome>_n<N>_seed<S>[_v2].rds`
  (the `_v2` suffix marks the sims rebuilt after the NDVI-ceiling and
  IMD-join fixes — current code always produces v2-equivalent sims).
- Fit outputs (current convention on disk):
  - `fit_medsat_p0_mgcv_local.rds`
  - `fit_medsat_p1_additive_n500_seed1_local.rds`
  - `fit_medsat_p1_additive_n500_seed1_diffuse_rho_local.rds`
  - `fit_medsat_p1plus_allpairs_n500_seed1_local.rds`
  - `fit_medsat_n1500_seed1_allpairs_hb.rds`  (the `_hb` tag = Hellbender)
- SBATCH scripts: `.sub` extension, under `scripts/submit/`.
- Plots: `output/plots/<runtag>/<what>.{pdf,png}` (per-run subfolder).
  Current subfolders: `p0_mgcv/`, `p1/`. New folder for the n=1500
  diagnostics: `p3/`.

---

## What's IN scope for this pilot

- Loading, filtering, scaling, subsampling MEDSAT London data.
- Running P0 (mgcv sanity), P1 (ours additive), P1plus (all-pairs n=500),
  P3 (all-pairs n=1500); eventually P4.
- Generating diagnostic plots: fitted smooths f_j(X_j) on original axes,
  fitted interaction surface f_{1,2}(NO2, NDVI), spatial random effect
  map over London, posterior trace/density for σ², τ², ρ, τ²_int.
- ESS / mixing diagnostics — especially on the interaction-variance
  chains, which are known to mix poorly at n=1500.

## What's OUT of scope for this pilot

- SPDE/Nyström acceleration (separate methodology project).
- Multi-outcome modelling (depression, anxiety, etc.) — pilot is asthma only.
- Multi-year modelling — pilot is 2019 only.
- The methodology's open issues (amplitude-dependent undersmoothing,
  boundary regularization, WAIC replacement) — those live in the
  methodology codebase.
- Refactoring the dissertation codebase. Apply patches if needed via
  surgical str_replace edits with explanation.
