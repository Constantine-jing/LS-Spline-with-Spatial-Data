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
| mgcv   | thin-plate / P-spline (`bs="ps"`) | k = 20          | bs="ps", m=c(2,2) -> cubic B-spline + 2nd-order penalty (RW2 analog); REML smoothing |
| inla   | RW2 (`f(..., model="rw2")`)  | bins = 50         | x binned onto 50-grid of [0,1], constr=TRUE, scale.model=TRUE; PC prior on prec |
| bayesx | P-spline RW2                 | TBD              | optional |

For the spatial term:

| Method | Spatial term         | Notes |
|--------|----------------------|-------|
| ours   | Matern GP, collapsed Gibbs | nu = 1 fixed |
| mgcv   | `s(lon, lat, bs="gp", m=c(2))` | nu = 1 fixed via m[1]=2; rho estimated by REML |
| inla   | SPDE alpha=2 via `inla.spde2.pcmatern` | nu = alpha - d/2 = 1 in 2D; PC priors on (range, sigma); inner mesh max.edge = 0.10 |
| bayesx | kriging               | nu = 1 fixed |

---

## 4. Known unavoidable differences

**Naming convention flip in the `ours` wrapper.** `gibbs_stage_c_full.R` uses the convention `b ~ N(0, sigma2*R)` and `eps ~ N(0, tau2*I)` -- i.e. its `sigma2` is the spatial variance and its `tau2` is the residual variance. The canonical schema (and the comparison plan) follows the more standard convention `b ~ N(0, tau2_s*R)` and `eps ~ N(0, sigma2*I)`. `fit_ours_wrapper.R` performs the rename when assembling `var_comp`. Documented in the wrapper's header comment.

**rho parameterization mismatch (ours vs mgcv vs INLA).** The "spatial range parameter" means a different scalar quantity in each method:

- *ours* (`spatial_utils.R::matern_cor`): `cor(d) = (2^(1-nu)/Gamma(nu)) * (d/rho)^nu * K_nu(d/rho)`. Here `rho` is the inverse-rate scale; effective range (correlation drops to 0.05) at nu=1 is approximately `4.0 * rho`. Truth `rho = 0.06` -> effective range ~ 0.24.
- *mgcv*: `bs="gp"` reports a range parameter that mgcv estimates internally; format is version-dependent and not always exposed.
- *INLA-SPDE*: reports the *practical range* `r = sqrt(8*nu)/kappa`, defined as the distance where the Matern correlation drops to ~0.13. For nu=1 this gives `r = sqrt(8) / kappa` ~ `2.83 / kappa`. To compare to ours: at nu=1 the relation `our_rho <-> INLA_range` is approximately `INLA_range ~= 2*our_rho * sqrt(8*nu)/2` ~ `4*our_rho * sqrt(2)` for the same effective spatial decay -- this conversion needs to be done carefully when reporting and is best handled by reporting BOTH numbers and their effective ranges, NOT a single "rho" number.

For Stage 1, each method reports its native parameterization. The variance-component-recovery row of the final comparison table will show effective range (a method-agnostic quantity) rather than the raw rho.

**INLA tau2_s vs ours tau2_s.** INLA's `inla.spde2.pcmatern` reports `Stdev for s_field` which is the marginal sd of the SPDE Gaussian field; squaring it gives a quantity *approximately* equivalent to our `tau2_s`. The mismatch is small for stationary models on bounded domains but not exactly zero. Reported as-is.

**INLA RW2 vs ours LS+RW2.** INLA's RW2 is a discrete random walk on integer indices; ours is LS basis weights with an RW2 prior on the coefficient sequence. Both penalize 2nd differences but the basis is different by construction. This is the SCIENTIFIC POINT of the comparison, not a bug.

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
