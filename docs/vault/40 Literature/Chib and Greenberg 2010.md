# Chib and Greenberg 2010

Additive cubic spline regression with Dirichlet process mixture errors.

## TL;DR

Uses the **LS basis** in a Bayesian additive regression framework. Demonstrates that natural cubic spline coefficients can be inferred under flexible (DP-mixture) error models.

## Why it matters here

This is the **closest precedent for using LS in Bayesian regression** — and the comparison Jimmy uses for **dissertation notation alignment** (translating the 1D construction into consistent dissertation notation).

But: **no spatial component**, **no RW2 smoothing prior**, **no interactions**. The pairing of LS + RW2 + spatial in this dissertation does not appear in C&G.

## Key differences from this dissertation

| | C&G 2010 | This dissertation |
|---|---|---|
| Basis | LS (same family) | LS (same family) |
| Smoothing | DP mixture on errors, no RW2 prior on coefficients | [[RW2 prior]] on coefficients |
| Spatial | None | [[Matérn GP random effect]] |
| Interactions | None | [[Tensor-product LS basis]] |

## Why this is useful for the writeup

C&G is the paper to cite when establishing **"the LS basis has been used in Bayesian regression before, but not with these other ingredients."** That framing is what makes [[Novelty stack]] contribution 1 (LS + RW2) credible — there's a precedent for half of the pairing, but not the full combination.

## Notation alignment

For the 1D case, the dissertation's notation follows C&G closely. The plan is to extend the same alignment to the [[Tensor-product LS basis]] / 2D case for consistency.

## Citation

Chib, S., & Greenberg, E. (2010). Additive cubic spline regression with Dirichlet process mixture errors. *Journal of Econometrics*, 156(2), 322–336.

(File on disk: `Additive_cubic_spline_regression_with_Dirichlet_process_mixture_errors.pdf`)

## Connects to

- [[LS basis]]
- [[Novelty stack]]
- [[Tensor-product LS basis]] (notation alignment open task)
