# Prior sensitivity B1-B5

Robustness study across five prior bundles. ✅ Complete. Result is publishable robustness evidence.

## Purpose

Show that posterior inference is stable across reasonable prior choices for the smoothing variance $\tau^2_{s,j}$ and spatial hyperparameters $(\sigma^2_b, \rho)$.

## The bundles

| Bundle | Description |
|---|---|
| B1 | Default (used everywhere else) |
| B2 | [fill in: e.g. tighter $\tau^2$ prior] |
| B3 | [fill in: looser] |
| B4 | [fill in] |
| B5 | **Deliberately misspecified** $\rho$ prior — stress test |

## Headline result

At $n=1000$, posterior summaries (RMSE, coverage, posterior mean of each component) are **near-identical across B1–B5, including B5**.

This is the core finding: **the model is robust to prior misspecification when the data are informative**.

## Why this matters for the dissertation

See [[Prior sensitivity is a strength]]. Reviewers will inevitably ask "are these results sensitive to your prior choices?" — this analysis is the answer.

## Settings

- Base simulation: [[Sim1 n1000]] design, or [fill in which]
- $n = 1000$
- All MCMC settings as in [[Sim1 n1000]]
- See also [[Prior sensitivity Sim4]] for the same analysis on the Sim4 design

Script: `run_prior_sensitivity.R` (and `run_prior_sensitivity_sim4.R`)

## Tables / plots to include

- Posterior mean curves across B1–B5 (visual overlay)
- RMSE table across B1–B5
- Coverage table across B1–B5
- Posterior summary of $\rho$ across bundles (most affected by B5)

[Insert specific numbers when you port them over]

## Connects to

- [[Prior sensitivity is a strength]]
- [[Prior sensitivity Sim4]]
- [[Sim1 n1000]]
