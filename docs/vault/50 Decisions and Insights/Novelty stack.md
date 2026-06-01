# Novelty stack

The four contributions of the dissertation. Useful as a one-page summary for the introduction, talks, and the defense.

## The four contributions

### 1. [[LS basis]] + [[RW2 prior]]

This pairing has not appeared before. [[Chib and Greenberg 2010]] uses LS basis in Bayesian regression but with different smoothing priors and no spatial component. The field has defaulted to P-splines + RW2.

**Why it's possible here:** the [[Collapsed Gibbs sampler]] handles the dense LS penalty matrix that would otherwise be awkward.

### 2. [[Matérn GP random effect]] in a fully Bayesian LS-spline framework

Existing work either uses:
- Continuous-location Matérn but **frequentist** ([[Nandy Lim Maiti 2017]]), or
- Bayesian spatial smoothing but with **discrete-region MRFs** ([[Lang and Brezger 2004]])

Neither is "Bayesian + LS basis + continuous Matérn GP" — that's the territory this dissertation occupies.

### 3. [[Collapsed Gibbs sampler]] marginalizing $b$

Not novel as a sampling technique (collapsed Gibbs is textbook). What's novel is:
- It's the enabling computational ingredient for contributions 1 and 2
- The empirical motivation ($\sigma^2$–$b$ coupling instability) is documented and clean
- $b$ is reconstructed post-hoc — no inferential loss

See [[Why collapsed Gibbs]].

### 4. [[Tensor-product LS basis]] pairwise interactions in a spatial setting

The most distinctive contribution. Tensor-product P-spline interactions exist ([[Eilers and Marx 2003]], [[Lang and Brezger 2004]]), but tensor-product **LS-spline** interactions in a **spatial Bayesian** setting do not — because the field hasn't pursued LS bases due to their density.

Validated in [[Interaction 2x2 v2]].

## How they reinforce each other

Contributions 1, 3, and 4 are linked:
- (1) is **enabled by** (3)
- (4) is **enabled by** (1) — once you have LS + RW2 working, the Khatri-Rao construction with [[Kronecker-sum RW2 penalty]] follows naturally

(2) is more independent — it's a "we put the right spatial component in this framework."

## Ordering for the writeup

For Chapter 1 introduction, the suggested order is:

1. Set up the methodological gap (LS basis underused; spatial Bayesian dominated by P-splines)
2. Present the framework (1 + 2 + 3 together)
3. Present the extension (4)
4. Validation: simulations + Munich rent

## The one-line pitch

> "Fully Bayesian additive geostatistical regression using an exact natural cubic spline basis and continuous-location Matérn dependence, made tractable by a collapsed Gibbs sampler, and extended to pairwise interactions via tensor products — none of which has appeared together in the literature."

## Connects to

- [[MOC - Dissertation]]
- [[Why LS over P-splines]]
- [[Why collapsed Gibbs]]
