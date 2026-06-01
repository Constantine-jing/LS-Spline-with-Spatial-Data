# Scenario A — Settings & Calibration Notes

This document records every setting used in the Scenario A
recovery comparison so that any cross-method difference can
be traced to a deliberate choice rather than an accident.

Update as new methods/wrappers are added.

---

## 1. Data-generating process (locked)

| Quantity      | Value                              | Notes |
|---------------|------------------------------------|-------|
| n             | 500                                | per the comparison plan |
| p (smooths)   | 3                                  | per the comparison plan |
| Covariates    | X_j iid Uniform(0,1)               | identical across methods |
| Locations     | (lon, lat) iid Uniform([0,1]^2)    | identical across methods |
| sigma^2       | 1.0                                | residual variance |
| tau^2_s       | 1.0                                | spatial signal variance |
| rho           | 0.06                               | Matern range parameter |
| nu            | 1.0                                | Matern smoothness — fixed |
| Eval grid 1D  | seq(0,1, length.out=101)           | identical across methods |
| Eval grid 2D  | 30 x 30 on [0,1]^2 (unused in A)   | kept for API symmetry |

### Truth functions (centered on [0,1])

| j | f_j(x) | sd on [0,1] |
|---|--------|-------------|
| 1 | 2 (sin(pi x) - 2/pi)                              | 0.616 |
| 2 | 1.5 (exp(x - 0.5) - (exp(0.5) - exp(-0.5)))       | 0.448 |
| 3 | 0.7 (x^2 - 1/3)                                   | 0.209 |

All three integrate to ~0 on [0,1] (verified to <1e-3 by trapezoidal
quadrature). This matches the centered-truth convention used in your
Sim1-Sim4+ work and is required for fair RMSE comparison across
methods that internally enforce sum-to-zero identifiability.

### Spatial calibration rationale

- **nu = 1** chosen for clean cross-method matching:
  - Our method: free to fix nu at any value via `matern_cor`.
  - mgcv: `s(lon, lat, bs="gp", m=c(3, ...))` -> nu = 1 in 2D.
  - INLA-SPDE: alpha = 2 -> nu = 1 in 2D (default).
  - BayesX: kriging with explicit nu argument.
  
  nu = 1.5 was the historical default in your work but cannot be
  exactly matched by INLA-SPDE; nu = 1 is the cleanest cross-method
  choice. Documented in conversation, April 2026.

- **rho = 0.06** chosen so the Matern effective range
  (correlation drops to 0.05) is approximately 0.24, i.e. ~25% of
  the spatial domain. This is the "modest spatial dependence"
  setting from the comparison plan. Verified numerically:
  `rho=0.06 -> effective range = 0.240`.

- **tau^2_s = 1.0** with sigma^2 = 1.0 gives spatial-signal-to-
  noise ratio of 1, which is "modest" relative to the overall
  signal sd ~ 0.79 (sum of three centered smooths). All three
  variance components are O(1), keeping the problem
  numerically stable.

---

## 2. Posterior representation (locked)

All four method wrappers must produce **joint posterior draws**
of every smooth, evaluated on the canonical grid.

- Our method, INLA: real MCMC / posterior draws.
- mgcv: `mvrnorm(n_draws, mu_grid, Vp_grid)` from the empirical
  Bayes posterior covariance Vp returned by `predict(..., type="lpmatrix")`.
- BayesX: real MCMC samples (when added).

`n_draws` (post-burn, post-thin) is fixed across methods at
**2000**.

For pointwise RMSE and pointwise 95% coverage these joint draws
give numerically identical results to direct interval-based
computation. Joint structure costs nothing extra and keeps the
schema clean. Documented in conversation, April 2026.

---

## 3. Knot count and basis settings

To be filled in once each wrapper is written.

| Method | Smooth term type             | Knot/basis count | Notes |
|--------|------------------------------|------------------|-------|
| ours   | LS natural cubic + RW2       | M = 20           | matches Sim1-Sim4+ default; n_iter=10500, n_burn=2500, n_thin=4 -> 2000 draws |
| mgcv   | thin-plate / P-spline (`bs="ps"`) | TBD          | will use `k=20` to match knot count |
| inla   | RW2 (`f(..., model="rw2")`)  | TBD              | needs grid for RW2 |
| bayesx | P-spline RW2                 | TBD              | optional |

For the spatial term:

| Method | Spatial term         | Notes |
|--------|----------------------|-------|
| ours   | Matern GP, collapsed Gibbs | nu = 1 fixed |
| mgcv   | `s(lon, lat, bs="gp", m=c(3, rho, nu))` | nu = 1 fixed via m argument |
| inla   | SPDE alpha=2          | nu = 1 in 2D |
| bayesx | kriging               | nu = 1 fixed |

---

## 4. Known unavoidable differences

**Naming convention flip in the `ours` wrapper.** `gibbs_stage_c_full.R` uses the convention `b ~ N(0, sigma2*R)` and `eps ~ N(0, tau2*I)` -- i.e. its `sigma2` is the spatial variance and its `tau2` is the residual variance. The canonical schema (and the comparison plan) follows the more standard convention `b ~ N(0, tau2_s*R)` and `eps ~ N(0, sigma2*I)`. `fit_ours_wrapper.R` performs the rename when assembling `var_comp`. Documented in the wrapper's header comment.

---

## 5. Convergence criteria

For MCMC methods (ours, INLA when in MCMC mode, BayesX):

- R-hat < 1.01 for sigma^2, tau^2_s, rho, and the tau^2_j vector.
- ESS > 400 for the same parameters.
- Re-run with longer chain if either is violated.

mgcv: REML convergence reported by `gam.fit`'s convergence flag.

---

## 6. Hardware / software

- Stage 1 (development): local Windows R 4.5, single core.
- Stage 3 (50 seeds): Hellbender (R 4.4 or 4.5), one seed per
  job, parallel across the cluster.
- Document R version, package versions, and BLAS library at
  paper-writing time (not now).
