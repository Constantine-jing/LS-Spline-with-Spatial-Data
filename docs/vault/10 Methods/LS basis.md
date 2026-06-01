# LS basis

Lancaster–Šalkauskas natural cubic spline basis. The core spline used throughout the dissertation.

## What it is

A basis for natural cubic splines parameterized so that **basis coefficients equal function values at knots**:
$$
f(\tau_m) = \theta_m
$$
This *knot-ordinate property* is what makes the [[Tensor-product LS basis]] construction clean.

## Why we use it (key claim)

LS yields the **exact** natural cubic spline solution, not an approximation.

Compare to P-splines (the field default): P-splines approximate via B-spline basis with difference penalty. They dominate because their sparse banded penalties fit GMRF / INLA infrastructure — a **computational** advantage, not a theoretical one. See [[Why LS over P-splines]].

The LS approach is made tractable here by [[Collapsed Gibbs sampler]], which sidesteps the sparsity advantage P-splines normally provide.

## How it appears in the model

Each smooth $f_j(X_{ij})$ is expanded as
$$
f_j(X_{ij}) = \sum_{m=1}^{M_j} Z_{j,im} \, \theta_{j,m}
$$
where $Z_j$ is the LS design matrix at the knots $\tau_{j,1}, \ldots, \tau_{j,M_j}$ and $\boldsymbol{\theta}_j$ is the coefficient vector.

## Combined with

- [[RW2 prior]] on $\boldsymbol{\theta}_j$ — never previously combined with LS basis (novelty contribution 1)
- [[Tensor-product LS basis]] for pairwise interactions
- [[ANOVA identifiability]] — sum-to-zero constraint for component identifiability
- [[Collapsed Gibbs sampler]] for posterior estimation

## Sources

- Theoretical: [[Lancaster and Šalkauskas 1986]]
- In Bayesian regression (1D, no spatial): [[Chib and Greenberg 2010]]

## Code

- Implementation: `ls_basis.R` — see [[R file index]]
