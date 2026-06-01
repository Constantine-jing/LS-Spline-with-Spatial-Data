# Eilers and Marx 2003

Multidimensional penalized signal regression / tensor-product P-spline penalty.

## TL;DR

Establishes the **Kronecker-sum penalty structure** for tensor-product P-spline interactions:
$$
K_{uv} = K_u \otimes I + I \otimes K_v
$$

This is the standard penalty form for 2D smoothing on a tensor-product basis.

## Why it matters here

Direct precedent for the [[Kronecker-sum RW2 penalty]] in this dissertation. The form is the same; the basis underneath is different (P-spline → LS).

## What we inherit

- The Kronecker-sum penalty form
- The justification (penalize curvature along each axis independently)
- The general tensor-product modeling philosophy

## What's different in this dissertation

The "Option C" derivation directly motivates the same penalty form from the [[LS basis]] knot-ordinate property — i.e., the penalty is on second differences of *function values* at the knot grid, not just on coefficient second differences. With LS this is the same thing (because $\theta_{(m,l)} = f(\tau_{u,m}, \tau_{v,l})$); with P-splines it isn't, and the connection is less clean.

This is documented in `ls_tensor_product_construction.md` and the PDF `ls_tensor_product_construction.pdf`.

## Citation

Eilers, P. H. C., & Marx, B. D. (2003). Multidimensional penalized signal regression. *Technometrics*, 45(3), 263–273.

[Verify exact title and year against the actual paper.]

## Connects to

- [[Kronecker-sum RW2 penalty]]
- [[Tensor-product LS basis]]
- [[Lang and Brezger 2004]] (which uses this penalty in Bayesian P-spline tensors)
