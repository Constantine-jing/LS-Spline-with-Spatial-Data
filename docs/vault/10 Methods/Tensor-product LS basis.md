# Tensor-product LS basis

The interaction-term construction. **Current focus of the dissertation (April 2026).**

## What it is

For two predictors $X_u$ and $X_v$ with [[LS basis]] design matrices $Z_u$ ($n \times M_u$) and $Z_v$ ($n \times M_v$), the bivariate smooth $f_{u,v}(X_u, X_v)$ uses the **row-wise Khatri-Rao product**:
$$
Z_{uv} = Z_u \odot Z_v \quad (n \times M_u M_v)
$$

The coefficient vector $\boldsymbol{\theta}_{uv}$ has dimension $M_u M_v$, with the **knot-ordinate property** preserved:
$$
\theta_{uv,(m,l)} = f_{u,v}(\tau_{u,m}, \tau_{v,l})
$$

## Why this is novelty contribution 4

Tensor-product **LS-spline** pairwise interactions in a **spatial** Bayesian setting — territory the field never pursued because it defaulted to P-splines for computational reasons.

[[Lang and Brezger 2004]] is the closest analog (P-spline tensor products in geoadditive Bayesian models) but uses P-splines and discrete-region spatial.

## The penalty: Option C

The 2D smoothness penalty is the [[Kronecker-sum RW2 penalty]]:
$$
K_{uv}^{(\theta)} = K_u^{(\theta)} \otimes I + I \otimes K_v^{(\theta)}
$$

**"Option C"** — the chosen approach — derives this directly from second differences of surface heights at knot grid points. Motivated by the LS knot-ordinate property: because $\theta_{ml}$ literally **is** $f$ evaluated at the knot grid, second differences in coefficient space equal second differences of the function on the grid.

Full derivation: `ls_tensor_product_construction.md` (27-step document) and the PDF `ls_tensor_product_construction.pdf`.

## Identifiability

The interaction must be orthogonal to the corresponding main effects, otherwise the ANOVA decomposition is non-identifiable.

Achieved via contrast matrices $T_u$ and $T_v$:
$$
\tilde{Z}_{uv} = (Z_u T_u) \odot (Z_v T_v)
$$
which removes the directions corresponding to $f_u$ and $f_v$ from the interaction space.

See [[ANOVA identifiability]] and [[Centering for ANOVA identifiability]].

## Validation

- [[Interaction null test]] — sanity check under $f_{u,v} \equiv 0$
- [[Interaction 2x2 v2]] — main 2×2 factorial validation, ✅ complete
- [[Interaction Phase 2]] — all $p(p-1)/2$ pairs, in progress

Key result from 2×2: after the centering fix (see [[Centering for ANOVA identifiability]]), RMSE_f12 improved 3–5× and coverage ≥90% in all six configurations. M1 cost ≤10% over M0.

## Code

- `ls_interaction.R` — basis construction (R)
- `gibbs_interaction.R` — sampler with interaction terms
- `ls_interaction_core.cpp` — Rcpp/RcppEigen hot loops
- Documentation: `ls_interaction_construction.md`, `ls_interaction_bayes_and_experiment.md`

## Open work

- Dissertation notation alignment (translate construction into consistent notation, like was done for [[Chib and Greenberg 2010]] in the 1D case)
- Phase 2 simulations
- Novelty/contribution framing relative to [[Lang and Brezger 2004]] for writeup

## Connects to

- [[LS basis]]
- [[Kronecker-sum RW2 penalty]]
- [[ANOVA identifiability]]
- [[Centering for ANOVA identifiability]]
