# Why collapsed Gibbs

Why we marginalize out $b$ rather than sampling it directly.

## The decision

In the [[Collapsed Gibbs sampler]], the spatial random effect $b \sim \mathcal{N}(\mathbf{0}, \sigma^2_b R(\rho))$ is **integrated out analytically** so the sampler operates on the marginal likelihood
$$
\mathbf{y} \mid \boldsymbol{\theta}, \sigma^2, \sigma^2_b, \rho \sim \mathcal{N}(Z\boldsymbol{\theta}, \, \sigma^2 I + \sigma^2_b R(\rho))
$$

A naive Gibbs sampler would update $b$ explicitly in each iteration. We do **not** do this.

## Why

**Empirical:** sampling $b$ explicitly produces $\sigma^2$–$b$ coupling instability. The chain mixes badly. Posterior summaries are unreliable. Confirmed in early experiments.

**Theoretical:** when components of the model are strongly correlated in the conditional updates, blocked or collapsed sampling is the textbook fix (Liu, 1994; van Dyk & Park, 2008). Here, $\sigma^2$ and $b$ are exactly such a pair — the model can absorb residual variation into either one.

## Why this is more than a convenience

This is a **core design decision**, on par with the basis choice itself:

1. It's what makes the [[LS basis]] practical. Without collapsing $b$, the dense LS penalty matrix would have to be paired with a Gibbs update on $b$ that *also* requires $n \times n$ work — no win over P-splines.
2. It enables [[Why LS over P-splines]] as a defensible framing.
3. $b$ can still be reconstructed post-hoc from posterior samples via the conditional mean formula (kriging-like), so nothing of inferential value is lost.

## What gets sampled instead

Only $(\boldsymbol{\theta}, \sigma^2, \sigma^2_b, \rho, \tau^2)$. See [[Collapsed Gibbs sampler]] for the full update list.

## Cost

The trade-off: each collapsed iteration requires an $n \times n$ Cholesky of $\sigma^2 I + \sigma^2_b R(\rho)$, which is $O(n^3)$. This caps tractable sample sizes around $n \approx 1000$ on a laptop. Chapter 3 (paused) was about scaling beyond this.

## Things to say

- ✅ "We use a collapsed Gibbs sampler that integrates out the spatial random effect. This avoids the $\sigma^2$–$b$ coupling instability we observed when sampling $b$ explicitly."
- ✅ "Marginalizing $b$ is essential for stable mixing and is what makes the dense LS penalty matrix computationally acceptable."

## Connects to

- [[Collapsed Gibbs sampler]]
- [[Why LS over P-splines]]
- [[Matérn GP random effect]]
