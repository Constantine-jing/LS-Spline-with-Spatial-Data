# CANONICAL_CHECK.md

Read-only canonical-source verification across the project + the new `06_medsat_pilot/` application. **Nothing was moved, renamed, edited, or deleted. Nothing was propagated.**

Generated: 2026-05-29.

## Method
- The medsat application lives at `06_medsat_pilot/` (note: `06_` prefix, not `medsat_pilot/`). Its shared-library files are under `06_medsat_pilot/R/`; its application-specific code is `06_medsat_pilot/R/prepare_medsat_london.R` plus everything in `06_medsat_pilot/scripts/`. **App-specific code is never treated as a canonical source and is excluded from propagation consideration.**
- For every core/shared file I computed MD5 of all project copies and the medsat copy, then diffed where hashes differed.
- medsat's copy is treated as a **comparison point, never an authority**. Canonical among *project* copies is chosen by feature-detection (per your rules), independent of medsat.

## Headline finding (changes your tiering assumption)
You expected `ls_interaction.R`, `gibbs_interaction.R`, and `fit_ours_interaction_wrapper.R` to be **Tier 2 (medsat lacks them)**. They are not. **medsat carries its own copy of all three** — so they are **Tier 1**. medsat is *not* silent on them:
- `ls_interaction.R` — medsat copy is **content-identical** to the project canonical (differs only by CRLF line endings).
- `gibbs_interaction.R` — medsat copy is a **real superset divergence** (adds a `fit_interactions` additive-only switch on top of the project canonical).
- `fit_ours_interaction_wrapper.R` — medsat copy is a **real divergence** (adds real-data "Option A" NA-truth fallback + `fit_interactions`).

Consequently **all 12 shared files are Tier 1**. There are **no Tier 2 files**. The feature-detection rules you supplied were still used — to pick the canonical *project* copy for the three multi-version files.

## Second finding: three "differing" hashes are EOL-only, not code
`canonical_schema.R`, `ls_interaction.R`, and `ls_interaction_core.cpp` have a different MD5 in medsat but **zero differing line content**. Byte-level check confirms the project copies are bare-LF (Unix) and the medsat copies are CRLF (Windows); the byte-size delta exactly equals the newline count (260 / 361 / 380). These are **not** code divergences. (The other 7 identical-hash medsat files are LF, matching the project exactly — so medsat's `R/` folder has mixed line endings internally, a minor cosmetic inconsistency.)

---

## Canonical decision table

| Filename | Tier | Proposed canonical path (project authority) | Stale copies to fix | medsat differences flagged for review | Notes |
|---|---|---|---|---|---|
| `baby_bayes.R` | 1 | `codes/run_codes/munich/baby_bayes.R` | none (also in `codes/BABY_BAYES_codes/` — byte-identical dup, not stale) | **None** — medsat hash `8EDC5C1B` identical | All 3 copies identical. |
| `canonical_schema.R` | 1 | `codes/test and comparison/canonical_schema.R` | none (`comparison_stage1/`, `scenario_A/` copies byte-identical dups) | **EOL only** (medsat CRLF, hash `34011640` vs project `76EAB373`); content identical line-for-line | Not a code divergence. |
| `fit_ours_interaction_wrapper.R` | 1 | `codes/test and comparison/fit_ours_interaction_wrapper.R` (feature: defaults `orthogonalize = TRUE`, records it in `fit$settings`) | `codes/interaction_old/comparison_stage1/fit_ours_interaction_wrapper.R` — **(a) stale/old** (322 lines, lacks the orthogonalize default; pre-patch) | **(b) deliberate app-specific divergence** — see detail below. medsat hash `6D1421FD` differs from BOTH project copies | 3 distinct versions exist. **YOU DECIDE** whether medsat's real-data additions get back-ported. |
| `fit_spatial_reml.R` | 1 | `codes/run_codes/munich/fit_spatial_reml.R` | none (5 byte-identical copies across project) | **None** — medsat hash `0F803D62` identical | Fully consistent everywhere. |
| `gibbs_bayes.R` | 1 | `codes/run_codes/munich/gibbs_bayes.R` | none (also in `codes/BABY_BAYES_codes/` — identical dup) | **None** — medsat hash `B8D2CD8B` identical | All 3 copies identical. |
| `gibbs_interaction.R` | 1 | `codes/run_codes/munich/gibbs_interaction.R` (= `test and comparison/` copy; feature: has clamp args `tau2_s_int_fixed`, `sigma2_fixed`, `rho_fixed`) | `codes/interaction_old/comparison_stage1/gibbs_interaction.R` and `codes/interaction_old/interaction_week1/gibbs_interaction.R` — **(a) stale/old** (452 lines, no clamp args) | **(b) deliberate app-specific divergence** (medsat is a 516-line superset adding `fit_interactions`); hash `98E86DB1` differs from project `9C0A1288` | **YOU DECIDE** whether `fit_interactions` switch belongs in canonical. |
| `gibbs_stage_c_full.R` | 1 | `codes/run_codes/munich/gibbs_stage_c_full.R` | none (5 byte-identical copies) | **None** — medsat hash `D22348BA` identical | Fully consistent everywhere. |
| `ls_basis.R` | 1 | `codes/run_codes/munich/ls_basis.R` | none (4 byte-identical copies) | **None** — medsat hash `6A3132B2` identical | Fully consistent everywhere. |
| `ls_interaction.R` | 1 | `codes/run_codes/munich/ls_interaction.R` (= `test and comparison/` copy; feature: `ls_build_interaction()` has the `orthogonalize` arg) | `codes/interaction_old/comparison_stage1/ls_interaction.R` and `codes/interaction_old/interaction_week1/ls_interaction.R` — **(a) stale/old** (289 lines, no `orthogonalize`) | **EOL only** (medsat CRLF, hash `734CC0D3`); content identical to project canonical | medsat agrees with canonical; no code review needed. |
| `ls_interaction_core.cpp` | 1 | `codes/run_codes/munich/ls_interaction_core.cpp` | none (4 byte-identical project copies) | **EOL only** (medsat CRLF, hash `4BE7F3D8` vs project `81F931C4`); content identical | Not a code divergence. |
| `marginal_utils.R` | 1 | `codes/run_codes/munich/marginal_utils.R` | none (4 byte-identical copies incl. `REML_codes/`) | **None** — medsat hash `04498862` identical | Fully consistent everywhere. |
| `spatial_utils.R` | 1 | `codes/run_codes/munich/spatial_utils.R` | none (4 byte-identical copies) | **None** — medsat hash `F0216DB4` identical | Fully consistent everywhere. |

> "Proposed canonical path" picks `codes/run_codes/munich/` for files that have a copy there (it is the active munich application, the closest analogue to the medsat application) and falls back to `codes/test and comparison/` otherwise. This is a *content* judgment (which version is authoritative); it does **not** assert that a single-shared-library layout should exist — that remains your open question from the prior audit. The medsat copies are evaluated only as comparison points and are excluded from "stale copies to fix."

---

## Detail on the two REAL medsat divergences (for your per-file decision)

### `fit_ours_interaction_wrapper.R` — medsat (319 ln) vs project canonical (`test and comparison`, 275 ln) — 126 differing lines
medsat is the patched-canonical wrapper **plus** application-support additions. Label: **(b) deliberate application-specific divergence**, with one sub-feature (`fit_interactions`) that is arguably general-purpose and currently ahead of the project canonical.

Quoted medsat-only additions:
- Real-data tolerance ("Option A") — falls back to data-driven MCMC inits when truth params are NA:
  ```
  # Use data-driven defaults for the REML init in that case...
  rho_init_val    <- if (is.na(sim$truth_params$rho))    0.2       else sim$truth_params$rho
  sigma2_init_val <- if (is.na(sim$truth_params$sigma2)) var(y) / 2 else sim$truth_params$sigma2
  tau2_s_init_val <- if (is.na(sim$truth_params$tau2_s)) var(y) / 2 else sim$truth_params$tau2_s
  ```
- Requires `nu` as a fixed modelling choice:
  ```
  if (is.null(sim$truth_params) || is.na(sim$truth_params$nu)) {
    stop("sim$truth_params$nu is required (fixed modelling choice).")
  ```
- `fit_interactions` flag gating the whole interaction-basis block:
  ```
  fit_interactions = TRUE
  ...
  if (S$fit_interactions) {
    int_list <- ls_build_all_interactions(full_objs, orthogonalize = S$orthogonalize)
    ...
  } else {
    int_list <- list(); ... ; H <- cbind(1, W_main)
  }
  ```
**Decision needed:** the real-data NA-init fallback is clearly medsat-application logic (do not propagate). But `fit_interactions` is a general capability the project canonical lacks — decide whether to back-port that one piece.

### `gibbs_interaction.R` — medsat (516 ln) vs project canonical (`run_codes/munich` = `test and comparison`, 498 ln) — 54 differing lines
medsat is a **superset** of the project canonical: it keeps the clamp diagnostic hooks and adds an additive-only switch. Label: mixed — driven by medsat's P1 additive runs **(b)**, but the switch itself is general and currently absent from canonical.

Quoted medsat-only additions:
```
# ---- ADDITIVE-ONLY SWITCH (added May 2026 for P1 cleanup) ----
# When FALSE, Step 5b (interaction tau2_s_int IG draw) is skipped...
fit_interactions = TRUE
...
if (!fit_interactions && n_int > 0L) {
  stop("fit_interactions = FALSE but col_map_int has length ", n_int, ...)
}
...
if (fit_interactions) {
  for (k in seq_len(n_int)) { ... tau2_s_int[k] <- 1 / rgamma(...) }
}
```
The project canonical's argument list ends `rho_fixed = NA_real_` (no `fit_interactions`); medsat appends `rho_fixed = NA_real_,` + `fit_interactions = TRUE` and wraps the interaction loop in `if (fit_interactions)`.
**Decision needed:** same `fit_interactions` capability as the wrapper above — decide jointly whether it is promoted into the project canonical.

---

## Files I can't decide for you (flagged, not resolved)
1. **`fit_ours_interaction_wrapper.R`** and **`gibbs_interaction.R`** — both carry the `fit_interactions` capability in medsat that the project canonical lacks. This is *not* purely application-specific (it's a modelling switch), so auto-classifying it as "don't propagate" would be wrong, but auto-promoting it would violate "never auto-prefer medsat." **Your call**, ideally decided for both files together since they're coupled.
2. **Canonical *home* for the 7 fully-identical files** (`baby_bayes`, `fit_spatial_reml`, `gibbs_bayes`, `gibbs_stage_c_full`, `ls_basis`, `marginal_utils`, `spatial_utils`) — content is unambiguous (one version everywhere), but *which directory* should be the single source of truth is still the open layout question from `CODE_AUDIT.md` §6 Q4. I proposed `run_codes/munich` paths above for concreteness only.

## Open items worth noting (no action taken)
- The `interaction_old/` copies of `ls_interaction.R`, `gibbs_interaction.R`, and `fit_ours_interaction_wrapper.R` are the **stale/old** versions per feature-detection — consistent with that folder's "old" name.
- medsat's `R/` folder mixes line endings (7 files LF, 3 files CRLF). Cosmetic; flag only.
- medsat application-specific files intentionally **not** assessed for canonical status: `06_medsat_pilot/R/prepare_medsat_london.R` and all 23 files under `06_medsat_pilot/scripts/`.
