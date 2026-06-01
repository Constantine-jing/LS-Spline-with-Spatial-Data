# RW2 prior

Second-order random walk prior on spline coefficients. Provides smoothness regularization.

## What it is

For coefficients $\boldsymbol{\theta}_j = (\theta_{j,1}, \ldots, \theta_{j,M_j})^\top$, the RW2 prior is
$$
\theta_{j,m} - 2\theta_{j,m-1} + \theta_{j,m-2} \stackrel{iid}{\sim} \mathcal{N}(0, \tau^2_{s,j})
$$
for $m = 3, \ldots, M_j$. This penalizes the **second differences** — i.e., curvature.

In matrix form: $\boldsymbol{\theta}_j \mid \tau^2_{s,j} \sim \mathcal{N}(\mathbf{0}, \tau^2_{s,j} \, K_j^{-})$ where $K_j$ is the rank-deficient RW2 precision matrix.

## Why RW2 over RW1

- RW1 penalizes first differences → favors **constant** functions
- RW2 penalizes second differences → favors **linear** functions
- The null space of RW2 is constants + linear trends — exactly what we want for smooth nonparametric components

## Hyperparameter

The smoothing parameter is $\tau^2_{s,j}$. We place an inverse-gamma prior on it. See [[Prior sensitivity B1-B5]] for robustness across choices.

## Combined with

- [[LS basis]] — the novel pairing (contribution 1 of [[Novelty stack]])
- [[Kronecker-sum RW2 penalty]] — the 2D extension for [[Tensor-product LS basis]]

## Why never combined with LS before

Field defaulted to P-splines (see [[Why LS over P-splines]]). With P-splines the RW2 penalty matrix is sparse and banded, which fits GMRF/INLA. With LS, the penalty matrix is dense — only made tractable here by the [[Collapsed Gibbs sampler]].

## Future: spike-and-slab on $\tau^2_{s,j}$

Chapter 2 will extend this to a spike-and-slab prior for component selection. See [[MOC - Dissertation]].
