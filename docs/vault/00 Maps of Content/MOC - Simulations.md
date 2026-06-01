# MOC — Simulations

Index of every simulation experiment. Status icons:
- ✅ complete
- 🔄 in progress
- ⏳ queued
- ❄️ paused

---

## Chapter 1 baseline (no interactions)

- ✅ [[Sim1 n1000]] — basic additive validation
- ✅ [[Sim2 n1000]] — [fill in purpose]
- ✅ [[Sim3 n1000]] — [fill in purpose]
- ✅ [[Sim4 n1000 MH1]] — [fill in purpose]
- ✅ [[Sim4 plus n1000]] — extended Sim4
- ✅ [[Bayes Sim4 n1000]] — full Bayesian on Sim4 setup

## Comparisons & robustness

- ✅ [[REML vs Bayes]] — frequentist vs Bayesian comparison
- ✅ [[Prior sensitivity B1-B5]] — robustness across prior bundles
- ✅ [[Prior sensitivity Sim4]] — same on Sim4 design
- ✅ [[Sample size comparison]] — n=400 vs n=1000 diagnostic

## Multi-seed & full Bayes runs

- ✅ [[Multi-seed run]]
- ✅ [[Full Bayes all]]

## Chapter 1 interaction extension (current focus)

- ✅ [[Interaction null test]] — null DGP sanity check
- ✅ [[Interaction 2x2 v2]] — 2×2 factorial design, the main validation
- 🔄 [[Interaction Phase 2]] — all p(p-1)/2 pairwise interactions

---

## Settings cheatsheet

| Knob | Typical value | Notes |
|---|---|---|
| n | 1000 | O(n³) Cholesky cap |
| M (knots, 1D) | 6, 20 | M=6 used in [[X1 undercoverage is basis resolution]] |
| MCMC iterations | [fill in] | |
| Burn-in | [fill in] | |
| Thinning | [fill in] | |
| Prior bundle | B1 (default) | See [[Prior sensitivity B1-B5]] |
