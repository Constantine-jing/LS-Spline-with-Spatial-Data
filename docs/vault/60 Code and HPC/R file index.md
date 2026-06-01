# R file index

What each R / C++ source file does and how it connects to the rest of the project.

## Core / utility

| File | Purpose |
|---|---|
| `ls_basis.R` | [[LS basis]] construction (1D) |
| `spatial_utils.R` | [[Matérn covariance function]] construction; spatial helpers |
| `marginal_utils.R` | Marginal likelihood evaluations for [[Collapsed Gibbs sampler]] |
| `fit_spatial_reml.R` | REML fit (frequentist comparison; [[REML vs Bayes]]) |

## Sampler implementations (1D additive, no interactions)

| File | Purpose |
|---|---|
| `baby_bayes.R` | Earliest / minimal Bayesian fit (likely deprecated) |
| `gibbs_bayes.R` | First-pass Bayesian Gibbs |
| `gibbs_bayes_v2.R` | Revised Gibbs |
| `gibbs_stage_c.R` | Staged Gibbs intermediate |
| `gibbs_stage_c_full.R` | **Main Chapter 1 sampler** |

## Sampler — interaction extension

| File | Purpose |
|---|---|
| `ls_interaction.R` | [[Tensor-product LS basis]] basis construction (R) |
| `gibbs_interaction.R` | Sampler with interaction terms |
| `ls_interaction_core.cpp` | Rcpp/RcppEigen hot loops for interaction sampler |

## Run scripts (one per experiment)

| File | What it runs |
|---|---|
| `run_sim1_n1000.R` | [[Sim1 n1000]] |
| `run_sim2_n1000.R` | [[Sim2 n1000]] |
| `run_sim3_n1000.R` | [[Sim3 n1000]] |
| `run_sim4_n1000_MH1.R` | [[Sim4 n1000 MH1]] |
| `run_sim4_plus_n1000.R` | [[Sim4 plus n1000]] |
| `run_bayes_sim4_n1000.R` | [[Bayes Sim4 n1000]] |
| `run_reml_vs_bayes.R` | [[REML vs Bayes]] |
| `run_full_bayes_all.R` | [[Full Bayes all]] |
| `run_multi_seed.R` | [[Multi-seed run]] |
| `run_prior_sensitivity.R` | [[Prior sensitivity B1-B5]] |
| `run_prior_sensitivity_sim4.R` | [[Prior sensitivity Sim4]] |
| `run_sample_size_comparison.R` | [[Sample size comparison]] |
| `run_interaction_null_test.R` | [[Interaction null test]] |
| `run_interaction_2x2_v2.R` | [[Interaction 2x2 v2]] |

## Documentation files (not code)

| File | Content |
|---|---|
| `ls_interaction_construction.md` | Tensor-product math (early version) |
| `ls_interaction_bayes_and_experiment.md` | Bayes + experiment description |
| `ls_tensor_product_construction.md` | **27-step construction document** |
| `ls_tensor_product_construction.pdf` | Same, typeset |

## Conventions

- Interactive Windows R 4.5 development; HPC runs use r/4.4.0
- Outputs saved as both **RDS** (full posterior) and **CSV** (tidy summary)
- Plots generated in the same loop immediately after each experiment
- Math-first, then code: new extensions get new files, not edits to existing ones
- `vapply` over `sapply` — see [[Bug log]]

## Connects to

- [[Hellbender HPC notes]] (running these on the cluster)
- [[Bug log]] (gotchas)
