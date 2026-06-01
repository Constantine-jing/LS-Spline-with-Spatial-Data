# Collapsed Gibbs sampler

The MCMC strategy used throughout. Marginalizes out the spatial random effect $b$ analytically.

## The core idea

The full hierarchy includes $b \sim \mathcal{N}(\mathbf{0}, \sigma^2_b \, R(\rho))$ where $R(\rho)$ is the [[Matérn covariance function]] evaluated at observed locations.

A naive Gibbs sampler updates $b$ explicitly. A **collapsed** Gibbs integrates $b$ out:
$$
\mathbf{y} \mid \boldsymbol{\theta}, \sigma^2, \sigma^2_b, \rho \sim \mathcal{N}\!\left(Z\boldsymbol{\theta}, \; \sigma^2 I + \sigma^2_b R(\rho)\right)
$$

This means the sampler only updates $(\boldsymbol{\theta}, \sigma^2, \sigma^2_b, \rho, \tau^2)$ — **never $b$ itself**.

## Why this is essential

**Empirical finding:** sampling $b$ explicitly causes $\sigma^2$–$b$ coupling instability — the chain mixes badly, posterior concentrates poorly.

This is a **core design decision**, not a convenience. See [[Why collapsed Gibbs]].

## Why this enables the LS approach

P-splines dominate the field partly because sparse banded penalty matrices fit GMRF / INLA. LS produces dense penalty matrices, which would be a problem under naive Gibbs. The collapsed Gibbs sidesteps that issue entirely. See [[Why LS over P-splines]].

## What gets sampled

In each iteration:
1. $\boldsymbol{\theta}_j$ — spline coefficients (Gaussian conditional, marginalized over $b$)
2. $\tau^2_{s,j}$ — smoothing parameter for each component (inverse-gamma)
3. $\sigma^2$ — error variance
4. $\sigma^2_b$ — spatial signal variance
5. $\rho$ — Matérn range parameter (Metropolis–Hastings step, see [[fill in MH details]])

After sampling, $b$ can be **post-hoc reconstructed** from posterior samples if needed (kriging-like formula).

## Cost

- Per iteration: $O(n^3)$ Cholesky of $\sigma^2 I + \sigma^2_b R(\rho)$
- This is why $n=1000$ is the upper bound for laptop runs — see [[Why n is capped at 1000]] [TODO]
- Chapter 3 (paused) was about scaling this beyond $n=1000$

## Code

- Main implementation: `gibbs_stage_c_full.R`
- Interaction extension: `gibbs_interaction.R`
- Hot loops in C++: `ls_interaction_core.cpp` via Rcpp/RcppEigen

See [[R file index]].

## Connects to

- [[Marginal likelihood derivation]] [TODO]
- [[Matérn GP random effect]]
- [[Why collapsed Gibbs]]
