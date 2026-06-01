# INLA as computational benchmark

Why INLA is the natural comparison and how to position the dissertation against it.

## The context

INLA (Integrated Nested Laplace Approximation) is the de facto tool for Bayesian geoadditive models with smooth additive components. It's fast, reliable, and built on a P-spline + GMRF foundation.

Anyone reading this dissertation will think: "Why not just use INLA?"

## The honest answer

INLA is excellent for what it does. It's not a competitor in capability; it's a competitor in *defaults*.

This dissertation's pitch is:

1. INLA's speed comes from sparse precision matrices → P-splines → an approximate basis
2. By using the [[LS basis]] + [[Collapsed Gibbs sampler]], we trade speed for the *exact* natural cubic spline solution under a Bayesian framework
3. We retain access to full MCMC tooling (component-wise spike-and-slab in Chapter 2, etc.) that INLA's approximations don't easily support

## What to actually compare

For the dissertation's Gibbs vs. INLA framing:

- **Point estimates:** should agree closely (same posterior, different computation)
- **Uncertainty:** should agree at $n=1000$; possibly differ at small $n$
- **Speed:** INLA wins by a lot
- **Flexibility for extensions:** Gibbs wins (spike-and-slab, custom priors, nonstandard likelihoods like Poisson)

## Updates to the original blueprint

[[Lang and Brezger 2004]] is the classic geoadditive Bayesian P-spline paper, and its simulation templates are still valid for the basic comparisons. But:

- DIC has been **superseded** by WAIC / LOO-CV — see [[Vehtari et al 2017]]. Use those, not DIC.
- INLA wasn't widely available when Lang & Brezger published; it now *is* the natural computational benchmark, more so than vanilla MCMC

## What to claim

- ✅ "INLA is the natural computational benchmark for this class of models. Our Gibbs sampler is slower per iteration but accommodates the dense LS penalty matrix and supports the spike-and-slab extension in Chapter 2 without further approximation."
- ✅ "Where direct comparison is feasible, posterior summaries agree closely."

## What to NOT claim

- ❌ "Our method is more accurate than INLA." (INLA is often very accurate)
- ❌ "Our method scales better than INLA." (it doesn't — see [[Why n is capped at 1000]] [TODO])

## Connects to

- [[Why LS over P-splines]]
- [[Lang and Brezger 2004]]
- [[Vehtari et al 2017]]
- [[Collapsed Gibbs sampler]]
