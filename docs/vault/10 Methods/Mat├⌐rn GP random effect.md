# Matérn GP random effect

The spatial component of the model. Captures dependence between observations at nearby locations.

## What it is

For observed locations $\mathbf{loc}_1, \ldots, \mathbf{loc}_n \in \mathbb{R}^2$:
$$
b \sim \mathcal{N}(\mathbf{0}, \sigma^2_b \, R(\rho))
$$
where $R(\rho)_{ij} = C_\nu(\|\mathbf{loc}_i - \mathbf{loc}_j\|; \rho)$ is the [[Matérn covariance function]] with smoothness $\nu$ and range $\rho$.

## Why Matérn

- Standard, defensible, well-studied family
- Smoothness $\nu$ controls differentiability of the spatial surface
- Range $\rho$ has interpretable meaning (effective spatial scale)
- Common in geostatistics — accessible to applied readers

## Why this is novelty contribution 2

Combining a Matérn GP on **continuous** locations with a **fully Bayesian LS-spline framework** has not appeared in the literature.

- [[Nandy Lim Maiti 2017]] uses Matérn with additive models, but **frequentist** penalized framework, not Bayesian
- [[Chib and Greenberg 2010]] uses LS basis in Bayesian regression, but **no spatial component**
- [[Lang and Brezger 2004]] does Bayesian geoadditive, but with **MRF on discrete regions**, not continuous-location Matérn

## How it appears in the model

$$
y_i = \beta_0 + \sum_{j=1}^{p} f_j(X_{ij}) + b_i + \varepsilon_i
$$
where $b_i$ is the $i$-th element of $b$.

Crucially, $b$ is **never sampled directly** — see [[Collapsed Gibbs sampler]].

## Hyperparameters

- $\sigma^2_b$: spatial signal variance
- $\rho$: range
- $\nu$: smoothness — typically fixed (often $\nu = 3/2$ or $\nu = 5/2$). Check current setting in `spatial_utils.R`.

## Code

- `spatial_utils.R` — covariance construction
- `marginal_utils.R` — marginal likelihood evaluations
- See [[R file index]]

## Connects to

- [[Collapsed Gibbs sampler]]
- [[Matérn covariance function]]
- [[Why collapsed Gibbs]]
