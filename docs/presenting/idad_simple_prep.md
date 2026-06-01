# IDAD 2026 — Simple Poster Prep Guide
## Jimmy Jing · April 21, 2026

---

## Where you are right now

Your poster (from Revvity Summit) has 3 columns:
- **Col 1**: Intro, Model Architecture diagram, Model Specification, Collapsed Gibbs
- **Col 2**: Simulation Design, RMSE tables, Key Findings, simulation figures
- **Col 3**: Munich Rent Application, Lang & Brezger comparison, Conclusions, References

You have REML (M=6 knots) and Full Bayes (M=20 knots, RW2 prior) implemented and validated across 5 simulation scenarios + Munich rent data. The interaction extension (Chapter 1 focus) is in progress but NOT on this poster.

---

## The basics: explain like they've never seen this

### What is a spline?

**Simple version:**
"You know how in a regression you fit a straight line? A spline is just a smooth curve instead. We break the covariate range into pieces and fit a smooth polynomial on each piece, and they join up smoothly at the boundaries."

**Health example:**
"Think about air pollution and asthma risk. At low pollution, maybe there's almost no effect. Then above some threshold the risk climbs fast, then maybe it plateaus. A straight line can't capture that shape. A spline can."

**If they know restricted cubic splines (RCS) from Cox models:**
"It's the same idea as the restricted cubic splines you might use in a Cox model — like with the rms package in R. We're using a specific type called the Lancaster–Šalkauskas basis that has a nice property: the spline coefficients directly equal the function values at the knot points."

### What are knots?

"Knots are the anchor points along the covariate range where the curve pieces join. More knots = more flexible curve. We use 20 knots per covariate."

"Think of it like pushpins on a flexible ruler — the ruler bends between the pins. 6 pins gives you a rough shape. 20 pins can capture detailed curves."

### Why LS basis instead of B-splines?

"In B-splines, the coefficients are abstract numbers — β₃ doesn't mean anything you can point to. In the LS basis, β₃ IS the function value at knot 3. So if we want to put a Bayesian prior on the coefficients, we can reason about it: 'adjacent function values should be similar' — that's our RW2 prior."

---

## What is REML?

**Simple version:**
"REML stands for Restricted Maximum Likelihood. It's a standard way to estimate variance components in mixed models. If you've ever fit a mixed model in SAS or R (like lmer), it probably used REML."

**What it does in our model:**
"We have three unknowns to estimate: how much spatial variation there is (σ²), how much noise there is (τ²), and how far the spatial correlation reaches (ρ). REML finds the values of these three parameters that maximize a special version of the likelihood that accounts for the fact that we also estimated the fixed effects."

**Our algorithm in plain English:**
1. Pick candidate values of ρ (range) and λ = τ²/σ² (noise-to-signal ratio)
2. For those values, build the covariance matrix Σ = σ²R(ρ) + τ²I
3. "Whiten" the data — transform it so the spatial correlation is removed
4. Do ordinary least squares on the whitened data to get β̂
5. The leftover residual variance gives us σ̂²
6. Search over (ρ, λ) to find the best fit
7. Once we have (ρ̂, σ̂², τ̂²), do a final GLS to get the best β̂

"It's like two-step GLS, but we search for the covariance parameters using REML instead of just plugging in moment estimates."

**Limitation:**
"REML gives you point estimates. It doesn't give you uncertainty bands on the curves or on the variance parameters. That's why we also do Full Bayes."

---

## What is MCMC?

**If they don't know:**
"MCMC = Markov Chain Monte Carlo. It's a way to explore a complicated probability distribution by taking random walks through it."

"Imagine you want to know the shape of a mountain range but you can't see it — you're in fog. MCMC is like wandering around the mountains, always tending to walk uphill. After enough walking, the places you visited most often trace out the shape of the peaks. Those peaks are where the most probable parameter values are."

**What it does in our model:**
"We want the posterior distribution — the probability of each parameter value given the data. The posterior is too complicated to write down exactly. So we use MCMC: at each iteration we update one parameter at a time from its conditional distribution, holding the others fixed. After thousands of iterations, the samples give us the posterior."

**Our sampler (6 steps each iteration):**
1. Update μ and β (the spline coefficients) — Gaussian, so we can sample directly
2. Update b (spatial effect) — also Gaussian
3. Update each τ²ₛ,ⱼ (smoothing variance per covariate) — Inverse-Gamma, direct
4. Update τ² (noise variance) — Inverse-Gamma, direct
5. Update σ² (spatial variance) — Inverse-Gamma, direct
6. Update ρ (spatial range) — this one is NOT a standard distribution, so we use Metropolis-Hastings: propose a new value, accept or reject based on the likelihood ratio

"Steps 1–5 are Gibbs sampling (we can sample directly). Step 6 is Metropolis-Hastings (we propose and accept/reject). The whole thing is called a Gibbs/MH hybrid."

**What is "collapsed"?**
"In a collapsed sampler, we integrate out (marginalize) one of the parameters before sampling. We integrate out the spatial effect b. Why? Because if you sample b and σ² separately in a naive Gibbs sampler, they get coupled — b grows large, σ² grows to compensate, and the chain gets stuck. By integrating out b, we break that coupling."

**Health analogy:**
"It's like a clinical trial where you have patient-level random effects and a treatment effect. If you try to separately estimate each patient's random effect AND the between-patient variance, they chase each other. It's better to integrate out the patient effects and estimate the variance from the marginal likelihood."

---

## How and why we set the priors

### RW2 prior on spline coefficients

"RW2 = second-order random walk. It says: adjacent spline coefficients should change smoothly. Specifically, the second differences β_{k} - 2β_{k-1} + β_{k-2} should be small."

"This is the Bayesian version of the penalized spline (P-spline) idea from Eilers and Marx. In the frequentist version you add a penalty λ∫f''(x)²dx. In the Bayesian version, that penalty becomes a prior. The smoothing parameter λ becomes a variance parameter τ²ₛ that we estimate from the data."

**Why RW2 not RW1?**
"RW1 penalizes first differences (slope changes). RW2 penalizes second differences (curvature changes). RW2 leaves linear functions completely unpenalized, which is usually what we want — we only want to shrink the wiggly parts."

### Per-covariate smoothing variance τ²ₛ,ⱼ

"Each covariate gets its own smoothing variance. Large τ²ₛ,ⱼ = wiggly curve (the prior is loose). Small τ²ₛ,ⱼ = smooth/flat curve (the prior pulls it toward a line)."

"The beautiful thing: for 'garbage' covariates with no true effect, the model learns τ²ₛ,ⱼ ≈ 0, pulling the curve flat. For real signals, τ²ₛ,ⱼ is larger. It's automatic variable selection without any penalty."

### Priors on variance components

| Parameter | Prior | Why |
|-----------|-------|-----|
| μ (intercept) | N(0, 10⁶) | Vague — let the data decide |
| τ²ₛ,ⱼ (smoothing) | IG(1, 0.005) | Weakly informative, allows both smooth and wiggly |
| σ² (spatial variance) | IG(2, 1) | Proper but not too tight |
| τ² (noise variance) | IG(2, 0.3) | Similar — weakly informative |
| log ρ (spatial range) | N(log 0.2, 1) | Centered on a reasonable range for the unit square |

"We tested five different prior specifications including a deliberately wrong one (prior mean ρ = 0.5 when true ρ = 0.2). Results are nearly identical — at n=1000 the data overwhelm the prior."

---

## How we set the knots

**REML model: M = 6 knots per covariate**
- Equally spaced from min to max of each covariate
- Few knots = less flexible, faster computation
- Problem: can't capture high-frequency oscillations (like 2sin(πx))

**Bayes model: M = 20 knots per covariate**
- Also equally spaced
- Many knots = very flexible basis
- Overfitting prevented by the RW2 prior (not by limiting knots)
- This is the P-spline philosophy: use plenty of knots, let the prior control smoothness

"Why 20? It's a standard choice in the P-spline literature. Ruppert (2002) showed that beyond about 20–30 knots the results are insensitive to the exact number, as long as you have a proper smoothing penalty."

"The jump from 6 to 20 knots is the single biggest reason our Bayes model beats REML for X₁ = 2sin(πx). Six knots simply can't represent that curve well enough."

---

## Why Bayes? The uncertainty argument

### The professional answer

"REML produces point estimates of the variance components (σ², τ², ρ) and then constructs confidence intervals for the curves conditional on those estimates. This ignores the estimation uncertainty in the covariance parameters. Our Bayesian approach samples all parameters jointly — the credible bands for f_j(·) integrate over the posterior uncertainty in σ², τ², ρ, τ²ₛ,ⱼ, and b simultaneously. The result is properly calibrated uncertainty that reflects all sources of variability in the model."

### The concrete example (use this — it lands well)

"Say you're studying PM2.5 and childhood asthma. Your model estimates an exposure-response curve. A policy maker asks: at what pollution level does risk start increasing?

With REML, you get a curve and a narrow band. It looks like risk jumps at 12 µg/m³. The band is tight, so the policy maker sets the threshold at 12.

But that narrow band was lying — it ignored uncertainty in the spatial correlation, the noise level, the smoothing. If you'd accounted for all that, the band would be wider, and the threshold could plausibly be anywhere from 8 to 15. A more conservative threshold at 8 might be warranted.

Narrow bands that ignore parameter uncertainty give false confidence — and in health policy, false confidence means setting the wrong threshold."

### The technical detail (if a statistician pushes)

"Formally, the REML-based interval for f_j(x) conditions on (σ̂², τ̂², ρ̂) and uses the sandwich or observed Fisher information. The Bayesian credible band marginalizes:

  p(f_j(x) | y) = ∫ p(f_j(x) | β, ...) p(β, σ², τ², ρ, τ²ₛ | y) dβ dσ² dτ² dρ dτ²ₛ

This integral is approximated by the MCMC samples. Each posterior draw of β gives a curve; the envelope of those curves is the credible band. No plug-in, no conditioning."

---

## The health/biostats connection

When someone asks "how is this relevant to health research?":

**Disease mapping:**
"We observe disease rates at clinic or county locations. The rate depends on risk factors like poverty, air quality, smoking prevalence. The effects might be nonlinear — a threshold effect for pollution, a U-shaped effect for BMI. And nearby locations have correlated rates even after adjusting for risk factors, due to unmeasured shared exposures or healthcare access patterns. Our model handles both: nonlinear covariate effects via splines, residual spatial dependence via the Matérn GP."

**Multi-pollutant exposure:**
"In environmental epi, people are exposed to multiple pollutants simultaneously. PM2.5 and ozone together might be worse than the sum of their individual effects. That's exactly what our interaction extension — LS tensor product surfaces — is designed for. This is directly relevant to multi-pollutant mixture analysis in EPA-funded research."

**Spatial confounding (if someone brings it up):**
"Yes, spatial confounding — where spatially-varying covariates are correlated with the spatial random effect — is a known issue. Our collapsed Gibbs helps with identifiability of the variance components by working on the marginal likelihood. But we don't claim to fully resolve it. Our Sim3 results show some σ² underestimation under confounding, consistent with the theoretical literature (Hodges & Reich 2010). We provide principled decomposition, not automatic separation."

---

## Quick answers to common questions

**Q: Why not just use a GAM (mgcv)?**
"mgcv is excellent for additive models and has spatial extensions (e.g., s(x,y) smooths). Three differences: (1) our credible bands integrate over all parameter uncertainty rather than conditioning on estimated smoothing parameters — this matters when the exposure-response curve informs a threshold decision; (2) per-covariate smoothing is fully Bayesian with adaptive τ²ₛ,ⱼ; (3) we use a continuous-location Matérn GP for the spatial effect, not a basis expansion, which is natural for irregularly-spaced point-referenced health data."

**Q: Why not INLA?**
"INLA requires the model to be expressible as a latent Gaussian model with sparse precision. B-splines + RW priors produce banded precision matrices, which is why they're the standard INLA choice. The LS basis produces a dense penalty matrix — it doesn't fit the sparse GMRF framework. Our collapsed Gibbs is the computational strategy that makes the LS basis tractable with a dense Matérn covariance."

**Q: How does this scale? What if n = 50,000?**
"The bottleneck is the O(n³) Cholesky of the n×n Matérn correlation matrix. At n=1000 each MCMC iteration takes ~0.1s; at n=3000 it's feasible on a cluster. For truly large n, approximate GP methods (nearest-neighbor GPs, predictive processes) would be needed. That's a natural future extension but not the focus of the current work."

**Q: What about binary or count outcomes?**
"Current implementation handles continuous responses. A Poisson log-link extension for count data is planned — we have SEER lung cancer mortality data ready for that application. The LS basis itself extends to non-Gaussian outcomes; Chib & Greenberg (2010) demonstrated it for ordinal data via a latent variable approach."

**Q: How do you choose the number of knots?**
"We follow the P-spline philosophy: use a generous number of knots (M=20) and let the RW2 prior control smoothness. Ruppert (2002) showed that beyond ~20 knots, results are insensitive to the exact count. The key insight is that in our framework, overfitting is controlled by the prior, not by limiting the basis dimension."

**Q: What software did you use?**
"R for the model framework, with the MCMC hot loops (Cholesky, matrix multiplications, MVN sampling) written in C++ via Rcpp/RcppEigen. Jobs run on our university HPC cluster (Hellbender). A typical run at n=1000, p=6, M=20, 10k iterations takes 30–45 minutes."

**Q: How do you check convergence?**
"Standard MCMC diagnostics: trace plots for all variance components and a subset of β coefficients, effective sample sizes, and Geweke's convergence diagnostic. The trace plots are on the poster — you can see the chains mix well for τ² and σ², with ρ showing slightly more autocorrelation due to the MH step."

**Q: Why not penalized likelihood / cross-validation instead of Bayes?**
"Three things. First, full posterior uncertainty — the credible bands propagate all sources of uncertainty, which matters for health policy decisions. Second, the per-covariate smoothing variance τ²ₛ,ⱼ is estimated jointly with all other parameters, no cross-validation needed. Third, the Bayesian framework provides a natural path to model comparison via marginal likelihoods and to extensions like spike-and-slab variable selection."

**Q: What is the Matérn GP?**
"The Matérn is a family of covariance functions for spatial data. It has a smoothness parameter ν that controls how rough the spatial surface is, and a range parameter ρ that controls how far the correlation reaches. We fix ν=1.5 (once differentiable, which is standard for environmental data) and estimate ρ. The key property: Matérn correlation decreases with distance, and the rate is controlled by ρ. It's the standard choice in geostatistics for continuous spatial locations."

---

## Your 30-second pitch

"I work on Bayesian spatial regression for health data. The problem: you observe a health outcome at geographic locations with multiple risk factors. The effects might be nonlinear, and nearby locations are correlated. My model recovers those nonlinear exposure-response curves while accounting for spatial dependence, with automatic adaptive smoothing and full uncertainty quantification. The novelty is combining the Lancaster–Šalkauskas spline basis with a Matérn Gaussian process spatial effect in a fully Bayesian framework, estimated via a collapsed Gibbs sampler."

---

## What to point to on the poster

1. **Model architecture diagram** (Col 1) — gives the big picture in 30 seconds
2. **Marginal curve plots** (Col 2) — "black dashed = truth, blue = our estimate, shaded = 95% credible band"
3. **RMSE table** (Col 2) — "3–4× improvement for oscillating functions"
4. **Munich rent curves** (Col 3) — "matches Lang & Brezger's benchmark, with full posterior uncertainty bands"

Start with the figures, not the math. The figures sell the method. The math is there for people who lean in and want the details.

---

## Things to keep in your back pocket (don't volunteer, but ready if asked)

- **Prior sensitivity results (B1–B5):** "Five prior bundles, including deliberately misspecified ρ, give near-identical posteriors. The data dominate at n=1000."
- **X₁ undercoverage (~74%):** "Approximation bias from basis resolution, not model failure. All other components achieve 95–100% coverage."
- **Interaction extension:** "We're building LS tensor product interaction surfaces for multi-pollutant mixture analysis — Khatri–Rao basis construction with Kronecker sum penalty." (Don't go deeper unless they're a potential collaborator.)
- **σ²–b coupling divergence:** "In a naive Gibbs sampler, b and σ² chase each other and the chain diverges. Marginalizing out b breaks the coupling. This is a well-known problem in spatial mixed models."
