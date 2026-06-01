# REML vs Bayes

Comparison of frequentist REML estimation vs the full Bayesian [[Collapsed Gibbs sampler]]. ✅ Complete.

## Purpose

Show that the Bayesian approach matches REML in point estimates but additionally provides coherent uncertainty quantification, including for nonlinear functionals and small-sample regimes.

## Setup

- REML implementation: `fit_spatial_reml.R`
- Bayesian: `gibbs_stage_c_full.R`
- Same data, same model structure
- Base design: [fill in which simulation — Sim1? Sim4?]

Script: `run_reml_vs_bayes.R`

## Headline findings

[Fill in from your actual results]

- Point estimates: REML and posterior mean agree closely
- Coverage: Bayesian credible intervals show [better / similar] frequentist coverage
- Uncertainty propagation: Bayesian naturally handles uncertainty in $\tau^2$, $\sigma^2_b$, $\rho$; REML plug-in approach underestimates this
- Computational cost: REML faster, Bayes more informative

## Why include this in the dissertation

- Demonstrates the model can be fit by either framework
- The Bayesian extras (full posterior, hierarchical priors → Chapter 2 spike-and-slab, prior sensitivity) are dissertation-justifying
- Provides a sanity check for both implementations

## Connects to

- [[Sim1 n1000]] (base design)
- [[REML estimation]]
- [[Collapsed Gibbs sampler]]
