# Lang and Brezger 2004

Bayesian P-Splines. Journal of Computational and Graphical Statistics.

## TL;DR

The classic Bayesian P-spline paper. Establishes Bayesian P-splines with RW2 priors, discrete-region geoadditive spatial models, and tensor-product P-spline interactions — all fit by MCMC.

## Why it matters here

This is the **primary blueprint for [[Tensor-product LS basis|the interaction extension]]** and the **field benchmark** the dissertation is in dialogue with. Their simulation templates remain valid for the basic Chapter 1 experiments.

If a reader asks "what is the standard approach this dissertation departs from?", the answer is "Lang and Brezger 2004 with extensions for spatial interactions."

## Key differences from this dissertation

| | Lang & Brezger 2004 | This dissertation |
|---|---|---|
| Basis | P-splines (B-splines + difference penalty) | [[LS basis]] (exact natural cubic spline) |
| Smoothing prior | RW2 on B-spline coefficients | [[RW2 prior]] on LS coefficients |
| Spatial | Markov random field on **discrete regions** | [[Matérn GP random effect]] on **continuous locations** |
| Estimation | MCMC, sample $b$ explicitly | [[Collapsed Gibbs sampler]], marginalize $b$ |
| Interactions | Tensor-product P-splines | [[Tensor-product LS basis]] (Khatri-Rao) |
| Model selection | DIC | WAIC / LOO-CV ([[Vehtari et al 2017]]) |

## Methodological inheritance

The dissertation borrows from L&B:
- Geoadditive structure: $y = \beta_0 + \sum f_j(X_j) + s(\text{loc}) + \varepsilon$
- RW2 prior on smoothing coefficients (with different basis)
- Tensor-product construction philosophy (with different basis)
- Simulation templates

The dissertation departs in the four ways listed in [[Novelty stack]].

## Citation

Lang, S., & Brezger, A. (2004). Bayesian P-splines. *Journal of Computational and Graphical Statistics*, 13(1), 183–212.

## Connects to

- [[Tensor-product LS basis]]
- [[RW2 prior]]
- [[Why LS over P-splines]]
- [[INLA as computational benchmark]]
- [[Novelty stack]]
