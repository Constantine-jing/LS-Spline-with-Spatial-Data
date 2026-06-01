# Nandy Lim Maiti 2017

Additive model building for spatial regression.

## TL;DR

Additive models with **Matérn covariance** spatial component, fit via a **frequentist penalized framework**. Closest analog on the spatial side, but not Bayesian.

## Why it matters here

This is the **comparison target** for the spatial component of the dissertation. Same Matérn-on-continuous-locations idea. Different inferential framework.

## Key differences from this dissertation

| | Nandy, Lim & Maiti 2017 | This dissertation |
|---|---|---|
| Inference | Frequentist (penalized likelihood) | Fully Bayesian |
| Spatial | Matérn on continuous locations ✓ | Matérn on continuous locations ✓ (same) |
| Basis | [fill in: B-spline / smoother of choice] | [[LS basis]] |
| Uncertainty | Plug-in / asymptotic | Full posterior |
| Component selection | [fill in] | [[RW2 prior]] now; spike-and-slab in Chapter 2 |
| Interactions | [fill in: do they have any?] | [[Tensor-product LS basis]] |

## What this paper establishes for our purposes

- Matérn-on-continuous-locations with additive smoothing is a viable, published modeling strategy
- The dissertation can position itself as "the Bayesian counterpart" — fully posterior, with hierarchical priors that enable Chapter 2 spike-and-slab without additional approximation

## Citation

[Fill in full citation when checking the PDF on disk]
File: `Additive_model_building_for_spatial_regression.pdf`

## Connects to

- [[Matérn GP random effect]]
- [[REML vs Bayes]]
- [[Novelty stack]] (contribution 2: Bayesian + LS + Matérn)
