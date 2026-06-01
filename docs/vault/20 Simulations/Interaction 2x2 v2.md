# Interaction 2x2 v2

The main validation experiment for the [[Tensor-product LS basis]] extension. ✅ **Phase 1 complete.**

## Purpose

Validate that the interaction extension:
1. Recovers $f_{u,v}$ when present (DGP-INT)
2. Doesn't hallucinate interaction when absent (DGP-NULL)
3. Doesn't degrade main-effect estimation
4. Has acceptable computational overhead

## Design

A 2×2 factorial:

| | DGP-NULL ($f_{12} \equiv 0$) | DGP-INT ($f_{12}$ active) |
|---|---|---|
| **M0** (main effects only) | sanity | misspecification check |
| **M1** (with interaction) | should not over-fit | should recover $f_{12}$ |

Crossed with $p \in \{2, 3, 4\}$ predictors → **6 configurations × 4 cells = 24 cells total**.

## Settings

- $n = 1000$
- M = [fill in knots per dimension]
- MCMC iterations: [fill in]
- Burn-in: [fill in]
- Prior bundle: [fill in — likely B1 default]
- Random seed: [fill in]

Script: `run_interaction_2x2_v2.R`

## Results

After the [[Centering for ANOVA identifiability]] fix:

- **RMSE_$f_{12}$ improved 3–5×** vs. uncentered version
- **Coverage ≥ 90%** in all six configurations (DGP-INT × $p$ ∈ {2,3,4} and DGP-NULL × $p$ ∈ {2,3,4})
- **M1 cost ≤ 10% over M0** (computational overhead is small)

[Add full RMSE / coverage table here when you port over the CSV summary]

| Config | RMSE_f1 | RMSE_f12 | Coverage_f1 | Coverage_f12 | Time |
|---|---|---|---|---|---|
| DGP-INT, M1, p=2 | | | | | |
| DGP-INT, M1, p=3 | | | | | |
| DGP-INT, M1, p=4 | | | | | |
| DGP-NULL, M1, p=2 | | | | | |
| DGP-NULL, M1, p=3 | | | | | |
| DGP-NULL, M1, p=4 | | | | | |

## Bugs encountered & fixed

- **Off-by-one in `col_map_main`** — fixed
- **`vapply` over `sapply`** — `sapply` on empty list returns list, breaking type assumption in `build_interaction_prior_precision`. See [[Bug log]].
- **Centering** — DGP $f_1 = 2\sin(\pi x)$ has nonzero mean on $[0,1]$; replaced with $2(\sin(\pi x) - 2/\pi)$. See [[Centering for ANOVA identifiability]].

## Compares to

- [[Sim1 n1000]] — pre-interaction baseline
- [[Lang and Brezger 2004]] — field benchmark for tensor-product interactions

## Status

- ✅ Phase 1 complete (X1×X2 only)
- ⏳ [[Interaction Phase 2]] — all $p(p-1)/2$ pairwise interactions

## Files

- Script: `run_interaction_2x2_v2.R`
- Sampler: `gibbs_interaction.R`
- Basis: `ls_interaction.R`
- C++: `ls_interaction_core.cpp`
- Outputs: [fill in path on disk to RDS / CSV / plots]

## Connects to

- [[Tensor-product LS basis]]
- [[Centering for ANOVA identifiability]]
- [[ANOVA identifiability]]
