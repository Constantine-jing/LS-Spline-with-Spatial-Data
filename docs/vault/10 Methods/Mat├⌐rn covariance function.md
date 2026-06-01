# Matérn covariance function

The kernel underlying the [[Matérn GP random effect]].

## The form

For two locations at distance $d = \|\mathbf{loc}_i - \mathbf{loc}_j\|$:
$$
C_\nu(d; \rho) = \frac{2^{1-\nu}}{\Gamma(\nu)} \left(\frac{\sqrt{2\nu}\, d}{\rho}\right)^\nu K_\nu\!\left(\frac{\sqrt{2\nu}\, d}{\rho}\right)
$$
where $K_\nu$ is the modified Bessel function of the second kind, $\nu > 0$ is the smoothness, and $\rho > 0$ is the range parameter.

## Useful special cases

- $\nu = 1/2$: exponential kernel $C(d) = \exp(-d/\rho)$ — paths are continuous but not differentiable
- $\nu = 3/2$: $C(d) = (1 + \sqrt{3}d/\rho)\exp(-\sqrt{3}d/\rho)$ — paths are once differentiable
- $\nu = 5/2$: $C(d) = (1 + \sqrt{5}d/\rho + 5d^2/(3\rho^2))\exp(-\sqrt{5}d/\rho)$ — twice differentiable
- $\nu \to \infty$: squared exponential / Gaussian kernel — infinitely differentiable

## What we use

[Fill in: typically $\nu$ is fixed — most likely $\nu = 3/2$ or $5/2$. Check `spatial_utils.R`.]

$\rho$ is sampled (Metropolis–Hastings step in [[Collapsed Gibbs sampler]]).

## Why fix $\nu$

Smoothness is weakly identified from typical spatial data. Conventionally fixed; a sensitivity analysis across $\nu \in \{1/2, 3/2, 5/2\}$ is sometimes reported.

## Code

`spatial_utils.R` constructs $R(\rho)$ given locations and the chosen $\nu$.

## Connects to

- [[Matérn GP random effect]]
- [[Collapsed Gibbs sampler]]
