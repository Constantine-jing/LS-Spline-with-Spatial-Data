# Vehtari et al 2017

Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.

## TL;DR

Establishes WAIC (widely applicable information criterion) and PSIS-LOO-CV (Pareto-smoothed importance sampling leave-one-out) as the **modern replacements for DIC** for Bayesian model comparison.

## Why it matters here

[[Lang and Brezger 2004]] — the dissertation's primary blueprint — used **DIC** for model comparison. DIC has known problems (not invariant to reparameterization, can fail for hierarchical models, can be negative). It's no longer the right default.

Use WAIC and/or PSIS-LOO-CV instead. Both are:
- Computable from posterior MCMC samples (no rerun needed)
- Better-justified theoretically
- Implemented in R via the `loo` package

## What this means for the dissertation

When reporting model comparison (M0 vs M1, with vs without component, prior bundle effects), use:
- **WAIC** as the primary criterion
- **PSIS-LOO-CV** as a check (with Pareto $\hat{k}$ diagnostics — a high $\hat{k}$ flags an observation where importance sampling is unreliable, i.e. a high-influence point)

Don't use DIC unless explicitly asked for legacy comparability.

## Code

The R `loo` package handles both. Inputs: pointwise log-likelihood matrix from MCMC samples ($n \times S$ where $S$ is post-burn-in iterations).

## Citation

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413–1432.

## Connects to

- [[Lang and Brezger 2004]] (uses DIC, which we replace)
- [[INLA as computational benchmark]]
- [[WAIC and LOO-CV]] (method note)
