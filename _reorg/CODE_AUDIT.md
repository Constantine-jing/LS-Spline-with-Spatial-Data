# CODE_AUDIT.md

Read-only structure map of the source code in this repository. **No files were modified, moved, renamed, or deleted.**

Scope: only source files (`*.R`, `*.cpp`, `*.h`, `*.md`, `*.sub`, `*.sh`). All `*.pdf`, `*.rds`, `*.csv`, `*.RData`, `*.png` and any `data/ output/ outputs/ results/ figures/` directories were excluded. (No `*.h` or `*.sh` files exist in the tree; `*.cpp`, `*.sub`, `*.md` and `*.R` do.)

Generated: 2026-05-29.

---

## 1. Directory tree (code files only)

```
LS-Spline-with-Spatial-Data/
├── README.md                          # 1-line stub only
├── baby_bayes_guide.md                # narrative guide to the baby-bayes / Gibbs prototypes
├── submit_boundary.sub                # HPC submit script (boundary experiment)
├── submit_exp1.sub                    # HPC submit script (experiment 1)
│
├── codes/
│   ├── hist_codes/                    # OLDEST prototypes (Jan 2026): staged LS-only / spatial-X scripts
│   ├── BABY_BAYES_codes/              # baby-bayes + first Gibbs samplers + their run scripts (prototype lineage)
│   ├── REML_codes/                    # frequentist REML fit + simulation sweeps (Sim1-5) + boundary/variable-selection
│   ├── bayes_n400_codes/              # full Bayesian runs at n=400 (reml-vs-bayes, multi-seed, prior sensitivity)
│   ├── bayes_n1000_mcmc10000_codes/   # full Bayesian runs at n=1000, 10k MCMC (sim1-4, prior sensitivity, sample-size)
│   │   ├── bayes_overlay/             # overlay-plot variants of the sim1-4 n1000 runs
│   │   └── bayes_test_sim1/           # one-off sim1 replotting experiment
│   ├── interaction_old/               # ARCHIVED interaction work (name says "old")
│   │   ├── comparison_stage1/         #   Scenario-A/B comparison harness: ours vs mgcv vs inla + many tests/diagnostics
│   │   │   └── scenario_A/            #     self-contained Scenario-A multiseed comparison (own copies of all deps)
│   │   └── interaction_week1/         #   first-week interaction prototypes (run_interaction_2x2*, null test)
│   ├── poisson_extension/             # Poisson/count-data extension of the Stage-C sampler
│   │   ├── michigan_lung_data/        #   real-data Poisson run (Michigan lung)
│   │   └── sim1_2_code/               #   Poisson simulations 1 & 2
│   ├── realtest_code/munich_rent/codes/ # Munich rent real-data analysis (HPC variants)
│   ├── run_codes/munich/              # NEWEST run-set (May 2026): current Munich pipeline w/ interaction + plotting
│   └── test and comparison/           # NEWEST (May 2026): Scenario-B "h3 / postfix / clamp" testing harness
│
├── interaction/interaction constract/ # math construction notes (markdown) for the tensor-product interaction basis
├── LS-Spline-Dissertation/            # Obsidian dissertation vault (markdown notes, MOCs, methods, literature, daily)
└── presenting files/                  # presentation prep notes (markdown)
```

Note the folder name `test and comparison` contains spaces — see Problem Flags.

---

## 2. File inventory

Line counts and last-modified dates are exact. Descriptions are inferred from headers, `source()` lines, and content.

### codes/hist_codes/ — oldest staged prototypes
| Path | Lines | Modified | Description |
|---|---|---|---|
| stage1_ls_only.R | 127 | 2026-01-29 | Stage 1a: LS-basis fit with no spatial random effect. |
| stage1b_ls_only.R | 254 | 2026-01-29 | Stage 1b LS-only with additive truth; sources `ls_basis.R`/`spatial_utils.R`. |
| stage2_spatialX.R | 121 | 2026-01-29 | Stage 2: adds spatial covariates + spatial residual GP. |

### codes/BABY_BAYES_codes/ — prototype Bayesian lineage
| Path | Lines | Modified | Description |
|---|---|---|---|
| baby_bayes.R | 265 | 2026-03-07 | Minimal Bayesian fit (eta posterior, variance fixed) — earliest sampler. |
| gibbs_bayes.R | 424 | 2026-03-07 | First-pass Gibbs sampler (samples b and variance components). |
| gibbs_bayes_v2.R | 230 | 2026-03-07 | **Collapsed/marginalized** Gibbs (b integrated out, MH on log-variances) — a *different* sampler, not a patch of v1. |
| gibbs_stage_c.R | 304 | 2026-03-07 | Staged Gibbs intermediate; defines RW2 penalty builders. |
| run_baby_bayes.R | 267 | 2026-03-07 | Run script for `baby_bayes.R`. |
| run_gibbs_bayes.R | 252 | 2026-03-07 | Run script for `gibbs_bayes.R`. |
| run_gibbs_v2.R | 194 | 2026-03-07 | Run script comparing v1 vs v2 collapsed sampler. |
| run_stage_c.R | 208 | 2026-03-07 | Run script for `gibbs_stage_c.R`. |
| run_stage_c_all.R | 319 | 2026-03-07 | Driver sourcing all baby-bayes samplers + stage_c. |
| submit_stage_c.sub | 12 | 2026-03-07 | HPC submit for stage_c. |

### codes/REML_codes/ — frequentist REML + simulation sweeps
| Path | Lines | Modified | Description |
|---|---|---|---|
| fit_spatial_reml.R | 131 | 2026-03-07 | REML fit of the spatial LS model (frequentist benchmark). |
| ls_basis.R *(absent here — see note)* | — | — | *(REML scripts source `ls_basis.R` but no copy lives in this folder)* |
| marginal_utils.R | 196 | 2026-03-07 | Marginal-likelihood evaluation helpers. |
| spatial_utils.R *(absent here)* | — | — | *(sourced by REML scripts; no local copy — relies on working dir)* |
| plot_utils.R | 170 | 2026-03-07 | Plotting helpers — **not sourced by any file (orphan).** |
| variable_selection.R | 145 | 2026-03-07 | Variable-selection routine; sources core utils + `sim5_sweep.R`. |
| boundary_experiments.R | 391 | 2026-03-07 | Boundary-behavior experiments driver. |
| run.R | 3 | 2026-03-07 | Tiny driver: sources `run_sim1/2/3_marginal.R`. |
| run_boundary.R | 278 | 2026-03-07 | Boundary experiment runner. |
| run_boundary_exp1.R | 190 | 2026-03-07 | Boundary experiment 1 runner. |
| run_sim1_marginal.R | 86 | 2026-03-07 | Sim1 marginal-likelihood run. |
| run_sim2_marginal.R | 92 | 2026-03-07 | Sim2 marginal-likelihood run. |
| run_sim3_marginal.R | 99 | 2026-03-07 | Sim3 marginal-likelihood run. |
| run_sim3_marginal_n1000.R | 88 | 2026-03-07 | Sim3 marginal run at n=1000. |
| sim1_sweep.R | 119 | 2026-03-07 | Sim1 parameter sweep (no spatial truth, spatial model fit). |
| sim2_sweep.R | 115 | 2026-03-07 | Sim2 parameter sweep. |
| sim3_sweep.R | 133 | 2026-03-07 | Sim3 parameter sweep. |
| sim4_sweep.R | 241 | 2026-03-07 | Sim4 parameter sweep. |
| sim5_sweep.R | 247 | 2026-03-07 | Sim5 parameter sweep (used by variable_selection). |
| Sim3_plot_surfaces.R | 180 | 2026-03-07 | Sim3 surface plotting script. |

> Note: REML run scripts `source("ls_basis.R")` / `source("spatial_utils.R")` but neither file exists in `REML_codes/`. They are present in other folders only. These scripts therefore only run if launched from a directory that has those files — a working-directory coupling, flagged below.

### codes/bayes_n400_codes/ — full Bayes, n=400
| Path | Lines | Modified | Description |
|---|---|---|---|
| gibbs_stage_c_full.R | 347 | 2026-03-07 | **Main Chapter-1 sampler** (full Stage-C Gibbs). |
| run_full_bayes_all.R | 298 | 2026-03-07 | Runs the full Bayes pipeline (all components). |
| run_multi_seed.R | 140 | 2026-03-07 | Multi-seed replication run. |
| run_prior_sensitivity.R | 156 | 2026-03-07 | Prior-sensitivity sweep B1–B5. |
| run_reml_vs_bayes.R | 211 | 2026-03-07 | REML-vs-Bayes comparison run. |
| submit_full_bayes.sub / submit_multi_seed.sub / submit_prior_sensitivity.sub / submit_reml_vs_bayes.sub | 13 each | 2026-03-07 | HPC submit scripts for the above. |

### codes/bayes_n1000_mcmc10000_codes/ — full Bayes, n=1000
| Path | Lines | Modified | Description |
|---|---|---|---|
| run_sim1_n1000.R | 248 | 2026-03-11 | Sim1 at n=1000. |
| run_sim2_n1000.R | 255 | 2026-03-11 | Sim2 at n=1000. |
| run_sim3_n1000.R | 270 | 2026-03-11 | Sim3 at n=1000. |
| run_sim4_n1000.R | 254 | 2026-03-07 | Sim4 at n=1000 (header says `run_sim4_n1000_MH1.R`). |
| run_sim4_plus_n1000.R | 267 | 2026-03-07 | Sim4-plus variant at n=1000. |
| run_bayes_sim4_n1000.R | 289 | 2026-03-07 | Bayes Sim4 n=1000 run. |
| run_prior_sensitivity_sim4.R | 372 | 2026-03-10 | Prior sensitivity for Sim4. |
| run_sample_size_comparison.R | 293 | 2026-03-07 | Sample-size comparison run. |
| submit_*.sub (8 files) | 13 each | 2026-03-07→11 | HPC submit scripts. |
| **bayes_overlay/** run_sim1/2/3/4/4_plus_n1000_overlay.R | 293–309 | 2026-03-13 | Overlay-plot variants of the sim runs. |
| **bayes_overlay/** sim2/3/4/4plus_n1000_overlay.sub | 13 each | 2026-03-13 | Submit scripts for overlays. |
| **bayes_test_sim1/** run_sim1_n1000_plot_new.R | 302 | 2026-03-13 | One-off Sim1 replot experiment. |
| **bayes_test_sim1/** submit_sim1_n1000_overlay.sub / _plot_new.sub | 13 | 2026-03-13 | Submit scripts. |

### codes/interaction_old/comparison_stage1/ — archived comparison harness
| Path | Lines | Modified | Description |
|---|---|---|---|
| canonical_schema.R | 255 | 2026-04-25 | Canonical output schema shared by all fit wrappers. |
| ls_basis.R | 329 | 2026-03-07 | LS-basis construction (1D). |
| ls_interaction.R | 289 | 2026-04-15 | **OLD** tensor-product interaction basis (no orthogonalize option). |
| ls_interaction_core.cpp | 333 | 2026-04-15 | Rcpp/Eigen Khatri-Rao hot loops for interaction basis. |
| gibbs_interaction.R | 452 | 2026-04-15 | **OLD** interaction sampler (no diagnostic clamp hooks). |
| gibbs_stage_c_full.R | 347 | 2026-03-07 | Main Stage-C sampler (copy). |
| fit_spatial_reml.R | 131 | 2026-03-07 | REML fit (copy). |
| marginal_utils.R | 196 | 2026-03-07 | Marginal-likelihood helpers (copy). |
| spatial_utils.R | 49 | 2026-03-07 | Matérn covariance + spatial helpers (copy). |
| simulate_scenario_A.R | 181 | 2026-04-25 | Scenario-A data generator (additive truth). |
| simulate_scenario_B.R | 200 | 2026-04-25 | Scenario-B data generator (interaction truth). |
| fit_ours_interaction_wrapper.R | 322 | 2026-04-25 | **OLD** wrapper for our interaction model (no orthogonalize). |
| fit_mgcv_interaction_wrapper.R | 287 | 2026-04-26 | mgcv interaction comparison wrapper. |
| diagnose_centering_bug.R | 100 | 2026-04-26 | Diagnostic for a centering bug. |
| diagnose_f12_absorption.R | 110 | 2026-04-26 | Diagnostic for f12 absorption issue. |
| test_b_smooth_loose.R | 143 | 2026-04-26 | Test: loose smoothing on b. |
| test_fit_mgcv_interaction.R | 102 | 2026-04-26 | Test harness for mgcv wrapper. |
| test_fit_ours_interaction.R | 133 | 2026-04-25 | Test harness for our wrapper. |
| test_fit_ours_interaction_ablation.R | 117 | 2026-04-25 | Ablation test for our wrapper. |
| test_interaction_on_additive_data.R | 122 | 2026-04-27 | Tests interaction model on additive (null) data. |
| test_simulate_scenario_B.R | 65 | 2026-04-25 | Test for scenario-B generator. |

### codes/interaction_old/comparison_stage1/scenario_A/ — self-contained Scenario-A subset
| Path | Lines | Modified | Description |
|---|---|---|---|
| run_stage1_scenarioA_multiseed.R | 176 | 2026-04-25 | Multiseed Scenario-A driver; conditionally sources ours/mgcv/inla wrappers. |
| simulate_scenario_A.R | 181 | 2026-04-25 | Scenario-A generator (copy). |
| compute_recovery_metrics.R | 159 | 2026-04-25 | Recovery-metric computation. |
| fit_ours_wrapper.R | 270 | 2026-04-25 | Our-model wrapper (Scenario-A, additive). |
| fit_mgcv_wrapper.R | 270 | 2026-04-25 | mgcv wrapper (Scenario-A). |
| fit_inla_wrapper.R | 326 | 2026-04-25 | INLA wrapper (Scenario-A). |
| canonical_schema.R | 255 | 2026-04-25 | Output schema (copy). |
| ls_basis.R / fit_spatial_reml.R / gibbs_stage_c_full.R / marginal_utils.R / spatial_utils.R | (copies) | 2026-03-07 | Copies of core utils. |
| test_fit_inla.R / test_fit_mgcv.R / test_fit_ours.R / test_simulate_scenario_A.R | 60–113 | 2026-04-25 | Per-wrapper test harnesses. |
| scenarioA_settings.md | 103 | 2026-04-25 | Scenario-A settings doc. |

### codes/interaction_old/interaction_week1/ — first-week prototypes
| Path | Lines | Modified | Description |
|---|---|---|---|
| ls_interaction.R | 289 | 2026-04-15 | OLD interaction basis (identical to comparison_stage1 copy). |
| ls_interaction_core.cpp | 333 | 2026-04-15 | Rcpp core (identical to all other copies). |
| gibbs_interaction.R | 452 | 2026-04-15 | OLD interaction sampler (identical to comparison_stage1 copy). |
| run_interaction_2x2.R | 840 | 2026-04-15 | 2×2 interaction experiment driver (large). |
| run_interaction_2x2_v2.R | 857 | 2026-04-16 | Revised 2×2 driver (sourced by p4 script). |
| run_interaction_2x2_p4.R | 316 | 2026-04-22 | p=4 variant; sources `run_interaction_2x2_v2.R`. |
| run_interaction_null_test.R | 264 | 2026-04-15 | Null (no-interaction) test driver. |
| scenarioA_settings.md | 96 | 2026-04-25 | Scenario-A settings doc (shorter copy). |

### codes/poisson_extension/ — count-data extension
| Path | Lines | Modified | Description |
|---|---|---|---|
| gibbs_stage_c_poisson.R | 357 | 2026-03-18 | Poisson-likelihood version of the Stage-C sampler. |
| sim1_2_code/run_sim_poisson_1.R | 231 | 2026-03-18 | Poisson simulation 1. |
| sim1_2_code/run_sim_poisson_2.R | 316 | 2026-03-18 | Poisson simulation 2. |
| sim1_2_code/run_sim_poisson_1.sub / _2.sub | 13 | 2026-03-18 | Submit scripts. |
| michigan_lung_data/run_michigan_lung_poisson.R | 368 | 2026-03-18 | Michigan lung real-data Poisson run. |
| michigan_lung_data/run_michigan_lung_poisson.sub | 13 | 2026-03-18 | Submit script. |

### codes/realtest_code/munich_rent/codes/ — Munich rent real data
| Path | Lines | Modified | Description |
|---|---|---|---|
| save_munich_data.R | 18 | 2026-03-13 | Loads/saves the Munich rent dataset to RDS. |
| run_munich_rent.R | 544 | 2026-03-13 | Munich rent analysis (n=full). |
| run_munich_rent_n1000.R | 559 | 2026-03-18 | Munich rent analysis, n=1000, 10k iter. |
| run_munich_rent_n1000_iter5000_testMH.R | 559 | 2026-03-18 | Same as above but n_iter=5000/burn=1500 (only those 2 lines differ). |
| submit_munich_rent.sub / submit_munich_rent_n1000.sub | 13 | 2026-03-13/15 | Submit scripts. |

### codes/run_codes/munich/ — NEWEST current Munich pipeline
| Path | Lines | Modified | Description |
|---|---|---|---|
| ls_interaction.R | 312 | 2026-05-06 | **NEW** interaction basis with `orthogonalize=TRUE` option. |
| ls_interaction_core.cpp | 333 | 2026-04-15 | Rcpp core (identical to all copies). |
| gibbs_interaction.R | 498 | 2026-05-09 | **NEW** interaction sampler with May-2026 diagnostic clamp hooks. |
| baby_bayes.R / gibbs_bayes.R / gibbs_stage_c_full.R / ls_basis.R / fit_spatial_reml.R / marginal_utils.R / spatial_utils.R | (copies) | 2026-03-07 | Copies of core/sampler utils. |
| run_munich_rent_local.R | 265 | 2026-05-15 | Local (non-HPC) Munich run wiring everything incl. interaction. |
| plot_modelA_checkpoint.R | 310 | 2026-05-18 | Plots a saved Model-A posterior checkpoint. |

### codes/test and comparison/ — NEWEST Scenario-B test harness
| Path | Lines | Modified | Description |
|---|---|---|---|
| ls_interaction.R | 312 | 2026-05-06 | NEW interaction basis (identical to run_codes/munich copy). |
| ls_interaction_core.cpp | 333 | 2026-04-15 | Rcpp core (identical). |
| gibbs_interaction.R | 498 | 2026-05-09 | NEW interaction sampler (identical to run_codes/munich copy). |
| fit_ours_interaction_wrapper.R | 275 | 2026-05-09 | **PATCHED** wrapper (`orthogonalize=TRUE`, clamp args) — differs from interaction_old copy. |
| gibbs_stage_c_full.R / ls_basis.R / fit_spatial_reml.R / spatial_utils.R | (copies) | 2026-03-07 | Copies of core utils. |
| canonical_schema.R | 255 | 2026-04-25 | Output schema (copy). |
| simulate_scenario_A.R | 181 | 2026-04-25 | Scenario-A generator (copy). |
| simulate_scenario_B.R | 200 | 2026-04-25 | Scenario-B generator (copy). |
| test_fit_ours_interaction.R | 133 | 2026-04-25 | Test harness (copy of comparison_stage1 version). |
| run_scenarioB_confirm_h3.R | 104 | 2026-05-07 | Confirm "Hypothesis 3". |
| run_scenarioB_tightbs.R | 112 | 2026-05-07 | Scenario-B with tight b-smoothing. |
| run_scenarioB_clamp_all_h3.R | 141 | 2026-05-08 | Clamp all variance components (H3). |
| run_scenarioB_clamp_main_int_h3.R | 183 | 2026-05-08 | Clamp main+interaction smoothing (H3). |
| run_scenarioB_clamp_full_h3.R | 167 | 2026-05-09 | Full clamp variant (H3). |
| run_scenarioB_postfix.R | 114 | 2026-05-09 | Post-fix Scenario-B run. |
| run_scenarioB_postfix_multiseed.R | 102 | 2026-05-10 | Multiseed post-fix run. |
| run_scenarioB_postfix_orthoFALSE.R | 106 | 2026-05-14 | Post-fix run with orthogonalize=FALSE. |
| run_mgcv_scenarioB_multiseed.R | 94 | 2026-05-14 | mgcv multiseed Scenario-B. |
| compare_h3_metrics.R | 245 | 2026-05-07 | Compare H3 metrics. |
| compare_seed1_metrics.R | 78 | 2026-05-07 | Compare seed-1 metrics. |
| compare_postfix_multiseed.R | 121 | 2026-05-10 | Compare post-fix multiseed metrics. |
| meeting_comparison.R | 197 | 2026-05-07 | Comparison figures/tables for a meeting. |
| plot_seed1_comparison.R | 144 | 2026-05-11 | Plot seed-1 comparison. |
| diagnose_orthogonalization.R | 172 | 2026-05-06 | Diagnostic for orthogonalization behavior. |
| verify_offby1_fix.R | 72 | 2026-05-09 | Verifies an off-by-one fix in the patched sampler. |

### Documentation / markdown
| Path | Lines | Modified | Description |
|---|---|---|---|
| README.md | 1 | 2026-02-09 | Stub (essentially empty). |
| baby_bayes_guide.md | 249 | 2026-03-03 | Narrative guide to the baby-bayes/Gibbs prototype lineage. |
| interaction/interaction constract/ls_interaction_construction.md | 325 | 2026-04-08 | Tensor-product interaction basis math (early). |
| interaction/interaction constract/ls_interaction_bayes_and_experiment.md | 305 | 2026-04-16 | Bayesian interaction + experiment write-up. |
| presenting files/idad_simple_prep.md | 149 | 2026-04-20 | Presentation prep notes. |
| LS-Spline-Dissertation/** (40+ .md) | 13–68 | 2026-04-28 | Obsidian dissertation vault: MOCs, methods, simulations, applications, literature, decisions, daily notes, templates. Includes `60 Code and HPC/R file index.md` — the author's own file-purpose index (a useful ground truth). |

> Two dissertation method-note filenames contain non-ASCII (`Matérn`) which rendered as mojibake in the directory listing (`Mat├⌐rn`) — encoding cosmetics, not a code issue.

---

## 3. Functional grouping

Every source file bucketed by what it **does**. Copies are grouped together.

**(a) REML / non-Bayesian**
- `fit_spatial_reml.R` (×5 identical copies: REML_codes, comparison_stage1, scenario_A, run_codes/munich, test-and-comparison)
- `REML_codes/`: sim1–5_sweep.R, run_sim1/2/3_marginal.R, run_sim3_marginal_n1000.R, run.R, boundary_experiments.R, run_boundary.R, run_boundary_exp1.R, variable_selection.R, Sim3_plot_surfaces.R, plot_utils.R, marginal_utils.R
- `hist_codes/`: stage1_ls_only.R, stage1b_ls_only.R, stage2_spatialX.R (LS-only, no Bayes)

**(b) baby-bayes / prototype**
- `baby_bayes.R` (×2: BABY_BAYES_codes, run_codes/munich)
- `gibbs_bayes.R` (×2), `gibbs_bayes_v2.R`, `gibbs_stage_c.R`
- run_baby_bayes.R, run_gibbs_bayes.R, run_gibbs_v2.R, run_stage_c.R, run_stage_c_all.R

**(c) Bayesian core / full model**
- `gibbs_stage_c_full.R` — **the main sampler** (×5 identical copies)
- `marginal_utils.R` (×4), `ls_basis.R` (×4), `spatial_utils.R` (×4) — shared math (also listed under (g))
- Run scripts: bayes_n400_codes/* and bayes_n1000_mcmc10000_codes/* (sim1–4, prior sensitivity, sample-size, reml-vs-bayes, multi-seed) and their overlay variants

**(d) interaction extension**
- `ls_interaction.R` — OLD 289-line (×2: comparison_stage1, interaction_week1) **and** NEW 312-line (×2: run_codes/munich, test-and-comparison)
- `gibbs_interaction.R` — OLD 452-line (×2) **and** NEW 498-line (×2)
- `ls_interaction_core.cpp` (×4, all identical)
- `fit_ours_interaction_wrapper.R` (comparison_stage1 322-line vs test-and-comparison 275-line patched), `fit_mgcv_interaction_wrapper.R`, `canonical_schema.R` (×3)
- Drivers: interaction_week1/run_interaction_2x2*.R, run_interaction_null_test.R; test-and-comparison/run_scenarioB_*.R
- scenario_A/ comparison harness: run_stage1_scenarioA_multiseed.R, fit_ours/mgcv/inla_wrapper.R, compute_recovery_metrics.R

**(e) diagnostics & testing**
- test-and-comparison: diagnose_orthogonalization.R, verify_offby1_fix.R, compare_h3_metrics.R, compare_seed1_metrics.R, compare_postfix_multiseed.R, plot_seed1_comparison.R, meeting_comparison.R
- comparison_stage1: diagnose_centering_bug.R, diagnose_f12_absorption.R, test_b_smooth_loose.R, test_*_interaction*.R, test_simulate_scenario_B.R
- scenario_A: test_fit_inla/mgcv/ours.R, test_simulate_scenario_A.R

**(f) real-data analysis**
- realtest_code/munich_rent/: run_munich_rent.R, run_munich_rent_n1000.R, run_munich_rent_n1000_iter5000_testMH.R, save_munich_data.R
- run_codes/munich/: run_munich_rent_local.R, plot_modelA_checkpoint.R
- poisson_extension/michigan_lung_data/run_michigan_lung_poisson.R

**(g) shared utilities**
- `ls_basis.R`, `spatial_utils.R`, `marginal_utils.R`, `canonical_schema.R`, `simulate_scenario_A.R`, `simulate_scenario_B.R`, `compute_recovery_metrics.R`, `plot_utils.R` — all duplicated across folders (see Problem Flags)

**(h) run scripts / SBATCH**
- All `*.sub` files (~25): one per experiment, 12–13 lines each, in nearly every folder

**(i) unclear / can't categorize**
- `poisson_extension/sim1_2_code/run_sim_poisson_2.R` vs `_1.R` — both "Poisson simulation", but what distinguishes sim 1 vs 2 (DGP? sample size?) is not stated in the header; needs a glance at the body to be sure.
- `REML_codes/run.R` (3 lines) — trivial driver; harmless but unclear if still the intended entry point.

---

## 4. Dependency map (what each entry script wires together)

`source()` targets are **bare filenames** (no paths), so every script depends on being launched from a directory that contains its dependencies. That is why core files are copied into each folder. Below, each runnable script → files it sources.

### REML_codes/ (deps live partly outside this folder — see flag)
- `run.R` → run_sim1_marginal.R, run_sim2_marginal.R, run_sim3_marginal.R
- `run_sim1/2/3_marginal.R`, `run_sim3_marginal_n1000.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R, marginal_utils.R
- `sim1/2/3_sweep.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R
- `sim4/5_sweep.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R, marginal_utils.R
- `variable_selection.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R, marginal_utils.R; → sim5_sweep.R
- `run_boundary.R`, `run_boundary_exp1.R`, `boundary_experiments.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R, marginal_utils.R, variable_selection.R
- `Sim3_plot_surfaces.R` → ls_basis.R, spatial_utils.R, fit_spatial_reml.R

### BABY_BAYES_codes/
- `run_baby_bayes.R` → ls_basis, spatial_utils, fit_spatial_reml, marginal_utils, baby_bayes
- `run_gibbs_bayes.R` → + gibbs_bayes
- `run_gibbs_v2.R` → + gibbs_bayes, gibbs_bayes_v2
- `run_stage_c.R` → + gibbs_bayes, gibbs_bayes_v2, gibbs_stage_c
- `run_stage_c_all.R` → + gibbs_bayes, gibbs_bayes_v2, gibbs_stage_c

### bayes_n400_codes/ + bayes_n1000_mcmc10000_codes/ (incl. overlay & test_sim1)
- All `run_*` scripts → ls_basis, spatial_utils, fit_spatial_reml, marginal_utils, baby_bayes, gibbs_bayes, gibbs_stage_c_full
  *(uniform 7-file dependency block; these folders do **not** contain local copies of those 7 files — they rely on the working directory)*

### realtest_code/munich_rent/codes/
- `run_munich_rent.R`, `run_munich_rent_n1000.R`, `run_munich_rent_n1000_iter5000_testMH.R` → same 7-file block as above (no local copies present)
- `save_munich_data.R` → (no source; standalone data writer)

### poisson_extension/
- `run_sim_poisson_1.R` / `_2.R` → ls_basis, spatial_utils, fit_spatial_reml, marginal_utils, gibbs_stage_c_poisson
- `run_michigan_lung_poisson.R` → ls_basis, spatial_utils, gibbs_stage_c_poisson
  *(no local copies of these deps in poisson_extension — relies on working dir)*

### interaction_old/interaction_week1/ (self-contained — has local copies)
- `run_interaction_2x2.R`, `run_interaction_2x2_v2.R`, `run_interaction_null_test.R` → spatial_utils, ls_basis, ls_interaction, gibbs_interaction
- `run_interaction_2x2_p4.R` → run_interaction_2x2_v2.R, then spatial_utils, ls_basis, ls_interaction, gibbs_interaction

### interaction_old/comparison_stage1/ (self-contained)
- `fit_ours_interaction_wrapper.R` → ls_basis, ls_interaction, spatial_utils, gibbs_stage_c_full, gibbs_interaction, canonical_schema; (fit_spatial_reml inside a branch)
- `fit_mgcv_interaction_wrapper.R` → canonical_schema
- `test_*` / `diagnose_*` → simulate_scenario_B + (fit_ours_interaction_wrapper or fit_mgcv_interaction_wrapper)
- `test_interaction_on_additive_data.R` → simulate_scenario_A, simulate_scenario_B, fit_ours_interaction_wrapper

### interaction_old/comparison_stage1/scenario_A/ (self-contained)
- `run_stage1_scenarioA_multiseed.R` → simulate_scenario_A, compute_recovery_metrics, then conditionally fit_ours_wrapper / fit_mgcv_wrapper / fit_inla_wrapper
- `fit_ours_wrapper.R` → ls_basis, spatial_utils, gibbs_stage_c_full, canonical_schema (+ fit_spatial_reml in branch)
- `fit_mgcv_wrapper.R`, `fit_inla_wrapper.R` → canonical_schema
- `test_fit_*` → simulate_scenario_A + corresponding wrapper

### run_codes/munich/ (self-contained, NEWEST)
- `run_munich_rent_local.R` → ls_basis, spatial_utils, fit_spatial_reml, marginal_utils, baby_bayes, gibbs_bayes, gibbs_stage_c_full, **ls_interaction, gibbs_interaction** (NEW versions)
- `plot_modelA_checkpoint.R` → ls_basis, spatial_utils, ls_interaction

### test and comparison/ (self-contained, NEWEST)
- `fit_ours_interaction_wrapper.R` (patched) → ls_basis, ls_interaction, spatial_utils, gibbs_stage_c_full, gibbs_interaction, canonical_schema (+ fit_spatial_reml in branch)
- `run_scenarioB_*.R` → simulate_scenario_B, fit_ours_interaction_wrapper (+ spatial_utils/fit_spatial_reml in the clamp_full variant)
- `run_mgcv_scenarioB_multiseed.R`, `compare_*`, `plot_seed1_comparison.R` → simulate_scenario_B
- `verify_offby1_fix.R` → simulate_scenario_B, gibbs_stage_c_full, ls_basis, ls_interaction, spatial_utils, fit_spatial_reml, gibbs_interaction, fit_ours_interaction_wrapper
- `diagnose_orthogonalization.R` → ls_basis, ls_interaction, simulate_scenario_B

---

## 5. Problem flags

### 5.1 Duplicate files (exact, identical content by MD5)
These are byte-identical copies scattered across folders:

| File | # copies | Locations |
|---|---|---|
| `fit_spatial_reml.R` | 5 | REML_codes, comparison_stage1, scenario_A, run_codes/munich, test-and-comparison |
| `gibbs_stage_c_full.R` | 5 | bayes_n400, comparison_stage1, scenario_A, run_codes/munich, test-and-comparison |
| `ls_basis.R` | 4 | comparison_stage1, scenario_A, run_codes/munich, test-and-comparison |
| `spatial_utils.R` | 4 | comparison_stage1, scenario_A, run_codes/munich, test-and-comparison |
| `marginal_utils.R` | 4 | REML_codes, comparison_stage1, scenario_A, run_codes/munich |
| `ls_interaction_core.cpp` | 4 | comparison_stage1, interaction_week1, run_codes/munich, test-and-comparison |
| `canonical_schema.R` | 3 | comparison_stage1, scenario_A, test-and-comparison |
| `simulate_scenario_A.R` | 3 | comparison_stage1, scenario_A, test-and-comparison |
| `simulate_scenario_B.R` | 2 | comparison_stage1, test-and-comparison |
| `baby_bayes.R` | 2 | BABY_BAYES_codes, run_codes/munich |
| `gibbs_bayes.R` | 2 | BABY_BAYES_codes, run_codes/munich |
| `test_fit_ours_interaction.R` | 2 | comparison_stage1, test-and-comparison |

These exist because `source()` uses bare filenames, so each runnable folder needs its own copy. The risk: editing one copy of `gibbs_stage_c_full.R` silently leaves four stale copies.

### 5.2 Near-duplicate files (same name, DIFFERENT content — diffed)

- **`gibbs_bayes.R` vs `gibbs_bayes_v2.R`** (BABY_BAYES_codes): **NOT versions of the same sampler.** v1 (424 lines) samples `b` and the variance components directly. v2 (230 lines) is a *collapsed/marginalized* Gibbs: it integrates `b` out (`y ~ N(Hη, σ²R + τ²I)`), samples η from the marginal, and uses random-walk Metropolis on `log σ²` / `log τ²`, deriving `b` afterward. The v2 header documents that v1's `σ²` could "explode" due to b–σ² coupling — v2 is the fix. They serve different roles and both are still sourced (`run_stage_c.R`, `run_gibbs_v2.R`).

- **`ls_interaction.R` — OLD (289 lines) vs NEW (312 lines):** the NEW version (run_codes/munich, test-and-comparison) adds an `orthogonalize = TRUE` option to `ls_build_interaction()` that residualizes the interaction design `W_uv` against `[1, W_u, W_v]` at the observations (data-level ANOVA orthogonality). `orthogonalize = FALSE` reproduces the OLD behavior. ~323 lines differ (largely new functions + header). The OLD copies live only under `interaction_old/`.

- **`gibbs_interaction.R` — OLD (452 lines) vs NEW (498 lines):** the NEW version adds May-2026 "diagnostic hook" arguments — `tau2_s_main_fixed`, `tau2_s_int_fixed`, `sigma2_fixed`, `tau2_fixed`, `rho_fixed` — that let you clamp variance components and skip their sampling steps (used for the "Hypothesis 3" undersmoothing investigation). 76 lines differ; otherwise identical. OLD copies only under `interaction_old/`.

- **`fit_ours_interaction_wrapper.R` — comparison_stage1 (322 lines) vs test-and-comparison (275 lines):** the test-and-comparison copy is the **patched** one — header states `orthogonalize = TRUE` default, forwards the new clamp args, and records `orthogonalize` in `fit$settings`. 199 lines differ. These are two generations of the same wrapper living in two folders.

- **`run_munich_rent_n1000.R` vs `run_munich_rent_n1000_iter5000_testMH.R`:** identical except **two lines** — `n_iter` (10000 → 5000) and `n_burn` (2500 → 1500). This is a copy-to-tweak-MCMC-length pattern; the "testMH" suffix is the only signal of intent.

### 5.3 Orphans / dead code (nothing sources them; not obviously a run entry point)
- **`REML_codes/plot_utils.R`** — defines plotting helpers but **no file sources it.** Either dead, or meant to be sourced interactively. Strongest orphan candidate.
- The entire **`interaction_old/`** subtree (week1 + comparison_stage1 + scenario_A) is superseded by the NEW interaction code in `run_codes/munich` and `test and comparison`. The folder name signals this, but the old code is still present and self-consistent, so it's "archived" rather than strictly dead — worth confirming nothing live still points here.
- **`hist_codes/`** (stage1/stage1b/stage2) — Jan-2026 prototypes, nothing sources them; effectively historical/dead.
- Top-level **`submit_boundary.sub` / `submit_exp1.sub`** and **`README.md`** (1 line) — the README is an empty stub; the two root `.sub` files presumably target REML boundary scripts but live far from them.

### 5.4 Inconsistent naming
- Folder **`test and comparison/`** has **spaces in its name** — fragile in shell/HPC contexts and inconsistent with every other underscore-named folder.
- Header-vs-filename mismatches: several files' top comment names a *different* file than the file itself, e.g. `run_sim4_n1000.R` header says `run_sim4_n1000_MH1.R`; the `bayes_overlay/run_sim*_overlay.R` headers say the non-overlay names; all three `run_munich_rent*` headers just say `run_munich_rent.R`. Headers were copied and not updated.
- Casing: `Sim3_plot_surfaces.R` (capital S) vs `sim3_sweep.R` (lowercase); commented self-references like `# source("Sim4_sweep.R")` vs the actual `sim4_sweep.R` — case-sensitive `source()` would break on Linux/HPC.
- `interaction_old/` contains two parallel `scenarioA_settings.md` files (96 vs 103 lines) in `interaction_week1/` and `scenario_A/`.
- `gibbs_bayes_v2.R` uses a `_v2` suffix while the *newer* `ls_interaction.R` / `gibbs_interaction.R` revisions carry **no** suffix and are distinguished only by folder — two different conventions for "this is the revised one."

### 5.5 Working-directory coupling (structural fragility, not a per-file bug)
`source()` everywhere uses bare filenames. Folders that contain their own copies (interaction_old/*, run_codes/munich, test-and-comparison, BABY_BAYES_codes) are self-runnable. But **REML_codes, bayes_n400_codes, bayes_n1000_mcmc10000_codes (+ overlay/test_sim1), poisson_extension, and realtest_code/munich_rent do NOT contain local copies of the core files they source** (`ls_basis.R`, `spatial_utils.R`, etc.). They only run if launched from a directory where those files are present — likely the HPC layout flattens everything, but in this repo layout they would fail. This is the single biggest "what's actually wired together" hazard.

### 5.6 Files whose purpose I genuinely can't pin down without reading the body
- `poisson_extension/sim1_2_code/run_sim_poisson_1.R` vs `run_sim_poisson_2.R` — what differs between Poisson sim 1 and 2 isn't in the header.
- `REML_codes/run.R` — 3-line driver; unclear if it's the current entry point or a leftover.
- The exact intended-canonical home: `run_codes/munich/` and `test and comparison/` both hold the NEW interaction code (identical), but `fit_ours_interaction_wrapper.R` exists *only* in test-and-comparison (patched) and comparison_stage1 (old) — not in run_codes/munich. Which folder is meant to be "the" current interaction pipeline is ambiguous.

---

## 6. Questions I'd need answered before any cleanup

1. **Which folder is the canonical "current" code?** `run_codes/munich/` and `test and comparison/` both carry the newest (May-2026) interaction code. Is one the source of truth and the other a working copy, or do they serve different purposes (real-data vs simulation)?
2. **Is `interaction_old/` truly archived?** Can the OLD 289-line `ls_interaction.R` / 452-line `gibbs_interaction.R` and the whole `comparison_stage1` + `scenario_A` + `interaction_week1` tree be moved to an `archive/` (or deleted), or are results in your dissertation still tied to that exact code?
3. **How do you actually launch jobs on HPC?** Do the `.sub` scripts `cd` into a flat directory where all sourced files coexist, or copy files in? This determines whether the duplicate copies are load-bearing or safe to collapse into one shared location.
4. **Do you want a single shared library folder?** i.e. one copy each of `ls_basis.R`, `spatial_utils.R`, `gibbs_stage_c_full.R`, etc., sourced via a path — or do you require each experiment folder to stay self-contained (current behavior) for reproducibility/HPC reasons?
5. **`gibbs_bayes_v2.R` and `baby_bayes.R` — still needed?** They're sourced by current run scripts (e.g. `run_munich_rent_local.R` sources `baby_bayes.R` and `gibbs_bayes.R`). Are these used for warm-starts/diagnostics, or vestigial sources you'd like to drop?
6. **`REML_codes/plot_utils.R` orphan** — keep (sourced interactively) or remove?
7. **`hist_codes/` (Jan-2026 stages)** — safe to archive/delete, or referenced anywhere in the writeup?
8. **The `run_munich_rent_n1000` vs `..._iter5000_testMH` pair** (differ only in MCMC length) — should MCMC length be a parameter/CLI arg instead of a forked file? Same question for the overlay vs non-overlay run families.
9. **Naming**: OK to rename `test and comparison/` to `test_and_comparison/` and fix the case-mismatched `source("Sim4_sweep.R")` references (Linux-case-sensitivity risk)?
10. **`poisson_extension`** — is this an active dissertation chapter or an exploratory side-branch? Affects whether it gets first-class treatment in any reorg.
