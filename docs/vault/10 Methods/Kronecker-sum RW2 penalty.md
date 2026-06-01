# Kronecker-sum RW2 penalty

The 2D smoothness penalty for [[Tensor-product LS basis]].

## The construction

For 1D RW2 precision matrices $K_u$ (size $M_u \times M_u$) and $K_v$ (size $M_v \times M_v$):
$$
K_{uv}^{(\theta)} = K_u \otimes I_{M_v} + I_{M_u} \otimes K_v
$$

This penalizes second differences along **both** axes of the knot grid. A surface flat in either direction has zero penalty contribution from that direction.

## Why Kronecker sum (not Kronecker product)

A Kronecker **product** $K_u \otimes K_v$ would penalize **mixed** second differences only — the surface could be wildly curved in one axis as long as it was curved-and-cancelling in the other. That's not what we want.

The Kronecker **sum** penalizes marginal curvature along each axis independently. This is the standard choice in tensor-product P-splines (see [[Eilers and Marx 2003]]).

## Why this fits LS naturally — "Option C"

The LS [[knot-ordinate property]] says $\theta_{(m,l)} = f(\tau_{u,m}, \tau_{v,l})$. So second differences of $\theta$ along the $u$ axis are *literally* second differences of the function values at the knot grid — there is no approximation step.

This is the cleanest theoretical motivation for the Kronecker-sum form when the basis is LS. Documented in detail in `ls_tensor_product_construction.md`.

## Computational note

This penalty matrix is **dense** in the LS parameterization (unlike the banded P-spline case). The [[Collapsed Gibbs sampler]] makes this tractable. See [[Why LS over P-splines]].

## Connects to

- [[Tensor-product LS basis]]
- [[RW2 prior]] (1D version)
- [[Eilers and Marx 2003]] — original Kronecker-sum penalty for P-splines
