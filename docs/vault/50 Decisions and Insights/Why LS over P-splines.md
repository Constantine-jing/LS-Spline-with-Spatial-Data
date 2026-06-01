# Why LS over P-splines

The single most important framing decision in the dissertation. Reviewers, talk audiences, and the defense committee will ask. This is the answer.

## The question

P-splines (Eilers & Marx, [[Lang and Brezger 2004]], etc.) dominate the field for smooth additive Bayesian regression. Why use the [[LS basis]] instead?

## The wrong answer

"LS is theoretically superior because it gives the exact natural cubic spline."

This is **technically true but misleading**. P-splines approximate the same target very well in practice. Saying "exact" implies P-splines are inadequate, which is not the case.

## The right answer (the framing)

**P-splines dominate for computational reasons, not theoretical superiority.**

Specifically:
- P-splines yield **sparse banded penalty matrices**, which fit GMRF / INLA infrastructure
- This makes them fast and the natural choice when you're plugging into existing software stacks
- LS in contrast yields **dense penalty matrices**, which are computationally awkward in standard MCMC

The LS approach in this dissertation is a **deliberate methodological alternative** made tractable by the [[Collapsed Gibbs sampler]], which sidesteps the sparsity advantage P-splines normally provide (the cost is dominated by the $n \times n$ Matérn covariance Cholesky, not the penalty matrix).

## Why this framing matters

1. **It's accurate.** It doesn't over-claim relative to what the code and methods actually do.
2. **It anticipates the obvious objection.** Anyone who has worked with smoothing in Bayesian models will think "but P-splines work fine — why bother?" — this answers them directly.
3. **It opens the door to the novelty stack.** Once you grant that LS is a viable alternative (not strictly better, but defensibly different), the contributions ([[LS basis]] + [[RW2 prior]] + [[Matérn GP random effect]] + collapsed Gibbs + tensor-product) line up naturally as new territory.

## Things to NOT say

- ❌ "LS is more accurate than P-splines"
- ❌ "P-splines are an approximation; ours is exact"
- ❌ "We chose LS for accuracy reasons"

## Things to say

- ✅ "We use the LS basis as a deliberate alternative to the P-spline default. The collapsed Gibbs sampler makes the dense LS penalty tractable, which is what enables the rest of the contribution stack."
- ✅ "The LS basis yields the exact natural cubic spline solution. This was inaccessible in standard MCMC because of the dense penalty; the collapsed Gibbs is what unlocks it."

## Connects to

- [[LS basis]]
- [[Collapsed Gibbs sampler]]
- [[Why collapsed Gibbs]]
- [[Novelty stack]]
- [[Lang and Brezger 2004]] (P-spline blueprint)
- [[INLA as computational benchmark]]
