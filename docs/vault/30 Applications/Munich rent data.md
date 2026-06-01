# Munich rent data

The real-data application for Chapter 1. ✅ Complete.

## Source

`gamlss.data::rent99` — Munich rent index data, $n = 3082$ observations.

Subsampled to $n = 1000$ for $O(n^3)$ Cholesky feasibility (see [[Collapsed Gibbs sampler]] cost discussion).

## Variables

[Fill in based on actual analysis]

- Response: rent (per square meter, or total)
- Continuous predictors: living area, year of construction, ...
- Spatial: district / location coordinates

## Why this dataset

- Standard benchmark in the geoadditive Bayesian literature
- Used in [[Lang and Brezger 2004]] and successors — directly comparable
- Covers the model's intended use: continuous covariates + spatial dependence

## Model fit

[Fill in headline findings]

- Posterior summaries for each smooth component
- Spatial component visualization (posterior mean of $b$ over locations)
- Hyperparameter posteriors: $\sigma^2$, $\sigma^2_b$, $\rho$, $\tau^2_{s,j}$
- WAIC / LOO-CV scores (see [[Vehtari et al 2017]])

## Comparison to baselines

- vs. linear model
- vs. additive model without spatial term
- vs. REML implementation (see [[REML vs Bayes]])

## Notes

- Subsampling to $n=1000$ is a real limitation. Phase 3 (paused) was about removing this cap.
- The interaction extension hasn't yet been applied to Munich rent — that's an open direction.

## Connects to

- [[MOC - Dissertation]]
- [[Collapsed Gibbs sampler]]
- [[REML vs Bayes]]
- [[Lang and Brezger 2004]]
