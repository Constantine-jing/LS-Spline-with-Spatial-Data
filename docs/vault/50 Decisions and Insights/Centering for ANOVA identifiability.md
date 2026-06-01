# Centering for ANOVA identifiability

A bug-and-fix story that became a methodological insight. Worth re-reading before any talk on the interaction extension.

## The bug

In the simulation DGPs for [[Interaction 2x2 v2]], the true $f_1$ was specified as
$$
f_1(x) = 2\sin(\pi x), \qquad x \in [0,1]
$$
This function has a **nonzero mean** on $[0,1]$:
$$
\int_0^1 2\sin(\pi x) \, dx = \frac{4}{\pi} \approx 1.273
$$

The fitted model imposes [[ANOVA identifiability]] — every component is centered. So the data-generating $f_1$ has a nonzero mean baked in, but the fitted $f_1$ is forced to be mean-zero.

**Consequence:** the $4/\pi$ mass got absorbed by the intercept (or worse, leaked into other components), and the comparison between fitted and true curves had a systematic offset → inflated RMSE, broken coverage.

## The fix

Replace the DGP with a **centered** version:
$$
f_1(x) = 2\!\left(\sin(\pi x) - \frac{2}{\pi}\right)
$$
This integrates to zero on $[0,1]$ and matches the constraint imposed during fitting.

After this fix, RMSE_$f_{12}$ improved 3–5× and coverage reached ≥90% in [[Interaction 2x2 v2]].

## The general principle

> **Whenever you simulate from an additive model and fit one with sum-to-zero constraints, your DGP components must also be mean-zero.**

If they're not, you're not actually testing what you think you're testing — you're testing the model's ability to recover a function that violates its own identifiability constraints, which is impossible by construction.

## Why this matters beyond the bug

1. It's a **lesson about reading the contract of your own model**. ANOVA identifiability is not just a constraint on inference — it's a constraint on what counts as "the truth" in any simulation.
2. It generalizes to interactions: $f_{u,v}$ in the DGP must satisfy the same row-sum-zero / column-sum-zero constraints that the fitted model imposes. Otherwise the "interaction" includes hidden main-effect mass.
3. It's the kind of thing that's easy to miss for months, then obvious in retrospect.

## What to do going forward

When designing any new DGP for this dissertation:

- For each main effect $f_j$: center it explicitly: $f_j(x) - \int f_j$ on the support
- For each interaction $f_{u,v}$: also remove the marginal main-effect projections (this is what the contrast matrices do for the *fitted* model — the truth must match)

## Connects to

- [[ANOVA identifiability]]
- [[Tensor-product LS basis]]
- [[Interaction 2x2 v2]]
- [[Interaction Phase 2]] — apply the same care
