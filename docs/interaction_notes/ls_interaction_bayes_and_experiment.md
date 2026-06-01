# The Bayesian Layer and Experimental Setup for LS Tensor Product Interactions

## Purpose

This document sits on top of `ls_tensor_product_construction.md` (the purely mathematical construction) and adds:

1. **Priors** — how each parameter gets a prior distribution
2. **Full conditionals** — the Gibbs sampler update rules
3. **Implementation map** — which equation corresponds to which line of code in `ls_interaction.R`, `gibbs_interaction.R`, and `run_interaction_2x2_v2.R`
4. **Experiment design** — the 2×2 factorial setup we built to validate the construction

No new mathematics of the spline construction itself — that's all in the construction document. Here we're specifying the probabilistic model and describing what we actually run.

---

## Part V: The Bayesian Hierarchical Model

### Step 28. The Full Likelihood

Starting from Step 26 of the construction document, the sampling model is:

$$
\mathbf{y} \mid \boldsymbol{\eta}, \mathbf{b}, \tau^2 \;\sim\; N\big(H\,\boldsymbol{\eta} + \mathbf{b},\; \tau^2 I_n\big)
$$

The spatial random effect is:

$$
\mathbf{b} \mid \sigma^2, \rho \;\sim\; N\big(\mathbf{0},\; \sigma^2 R(\rho, \nu)\big)
$$

where $R(\rho, \nu)$ is the $n \times n$ Matérn correlation matrix with range $\rho$ and smoothness $\nu$ (fixed, typically $\nu = 3/2$). Marginalizing out $\mathbf{b}$ gives the **collapsed likelihood** used throughout the sampler:

$$
\mathbf{y} \mid \boldsymbol{\eta}, \sigma^2, \tau^2, \rho \;\sim\; N\big(H\,\boldsymbol{\eta},\; \Sigma\big), \qquad \Sigma = \sigma^2 R(\rho, \nu) + \tau^2 I_n
$$

This is the same collapsed setup from `gibbs_stage_c_full.R`. Adding interactions changes $H$ (more columns) and changes the prior on $\boldsymbol{\eta}$ (more blocks). Nothing else changes.

### Step 29. Priors on $\boldsymbol{\eta}$ — Block Structure

The coefficient vector partitions into blocks:

$$
\boldsymbol{\eta} = \big(\mu,\; \beta_1,\; \ldots,\; \beta_p,\; \beta_{12},\; \beta_{13},\; \ldots,\; \beta_{p-1,p}\big)
$$

Each block gets its own prior:

**Intercept** (vague normal):

$$
\mu \sim N(0, \kappa^2), \qquad \kappa^2 = 10^6
$$

**Main effect $j$** (Gaussian with RW2 precision, smoothing variance $\tau^2_{s,j}$):

$$
\beta_j \mid \tau^2_{s,j} \;\sim\; N\!\left(\mathbf{0},\; \tau^2_{s,j}\, \big(K_j^{(\beta)}\big)^{-}\right)
$$

The prior precision is $\tau^{-2}_{s,j} K_j^{(\beta)}$. Because $K_j^{(\beta)}$ is rank-deficient (rank $M_j - 2$ out of $M_j - 1$), we use the **partially improper** prior — the rank-1 null space of $K_j^{(\beta)}$ corresponds to the single linear direction that isn't penalized, and that direction is given an implicit flat prior.

**Interaction block $(u,v)$** (same form, new 2D penalty, separate smoothing variance):

$$
\beta_{uv} \mid \tau^2_{s,uv} \;\sim\; N\!\left(\mathbf{0},\; \tau^2_{s,uv}\, \big(K_{uv}^{(\beta)}\big)^{-}\right)
$$

where $K_{uv}^{(\beta)} = K_u^{(\beta)} \otimes G_v + G_u \otimes K_v^{(\beta)}$ from Step 21 of the construction document.

**Full prior precision matrix on $\boldsymbol{\eta}$** (block diagonal):

$$
Q_0 = \text{blockdiag}\!\left(\frac{1}{\kappa^2},\; \frac{1}{\tau^2_{s,1}} K_1^{(\beta)},\; \ldots,\; \frac{1}{\tau^2_{s,p}} K_p^{(\beta)},\; \frac{1}{\tau^2_{s,12}} K_{12}^{(\beta)},\; \ldots\right) + \epsilon_{\text{ridge}}\, I
$$

The tiny ridge $\epsilon_{\text{ridge}} = 10^{-6}$ ensures $Q_0$ is strictly positive definite (needed for Cholesky in the sampler).

**Code reference:** `build_interaction_prior_precision()` in `gibbs_interaction.R` lines 68–95.

### Step 30. Priors on Variance Parameters

**Smoothing variances** (conjugate inverse gamma, same hyperpriors for main effects and interactions):

$$
\tau^2_{s,j} \;\sim\; \text{IG}(a_s, b_s), \qquad \tau^2_{s,uv} \;\sim\; \text{IG}(a_s, b_s)
$$

with $a_s = 1$ and $b_s = 0.005$ (weakly informative, same as main-effect smoothing hyperpriors in `gibbs_stage_c_full.R`).

**Nugget variance:**

$$
\tau^2 \;\sim\; \text{IG}(a_\tau, b_\tau), \qquad a_\tau = 2,\; b_\tau = 0.3
$$

**Spatial variance:**

$$
\sigma^2 \;\sim\; \text{IG}(a_\sigma, b_\sigma), \qquad a_\sigma = 2,\; b_\sigma = 1
$$

**Matérn range** (log-normal):

$$
\log \rho \;\sim\; N(\log \rho_\mu, \log \rho_\sigma^2), \qquad \log \rho_\mu = -1.6,\; \log \rho_\sigma = 1.0
$$

**Matérn smoothness:** $\nu$ is fixed at 1.5 (standard practice — hard to identify $\nu$ from spatial data).

### Step 31. The Joint Posterior

Bringing everything together, the joint posterior is proportional to:

$$
p(\boldsymbol{\eta}, \sigma^2, \tau^2, \rho, \{\tau^2_{s,j}\}, \{\tau^2_{s,uv}\} \mid \mathbf{y}) \;\propto\; p(\mathbf{y} \mid \boldsymbol{\eta}, \sigma^2, \tau^2, \rho) \cdot p(\boldsymbol{\eta} \mid \{\tau^2_s\}) \cdot p(\{\tau^2_s\}) \cdot p(\sigma^2) \cdot p(\tau^2) \cdot p(\rho)
$$

$$
\propto \underbrace{|\Sigma|^{-1/2} \exp\!\left(-\tfrac{1}{2}(\mathbf{y} - H\boldsymbol{\eta})^\top \Sigma^{-1} (\mathbf{y} - H\boldsymbol{\eta})\right)}_{\text{collapsed likelihood}} \times \underbrace{\exp\!\left(-\tfrac{1}{2}\boldsymbol{\eta}^\top Q_0\, \boldsymbol{\eta}\right)}_{\text{prior on } \boldsymbol{\eta}} \times \prod \text{IG and logN factors}
$$

---

## Part VI: Full Conditionals and the Gibbs Sampler

### Step 32. Update for $\boldsymbol{\eta}$ — Multivariate Normal (conjugate)

The collapsed likelihood is Gaussian in $\boldsymbol{\eta}$, the prior is Gaussian in $\boldsymbol{\eta}$, so the full conditional is Gaussian:

$$
\boldsymbol{\eta} \mid \mathbf{y}, \sigma^2, \tau^2, \rho, \{\tau^2_s\} \;\sim\; N\big(\mathbf{m}_\eta,\; Q_\eta^{-1}\big)
$$

with **posterior precision**:

$$
Q_\eta = H^\top \Sigma^{-1} H + Q_0
$$

and **posterior mean**:

$$
\mathbf{m}_\eta = Q_\eta^{-1}\, H^\top \Sigma^{-1} \mathbf{y}
$$

**Sampling in practice.** Directly inverting $Q_\eta$ is unstable. Instead:

1. Cholesky $L_\Sigma L_\Sigma^\top = \Sigma$. Let $\tilde{\mathbf{y}} = L_\Sigma^{-1} \mathbf{y}$, $\tilde H = L_\Sigma^{-1} H$ (whitening).
2. Then $Q_\eta = \tilde H^\top \tilde H + Q_0$.
3. Cholesky $U_\eta^\top U_\eta = Q_\eta$. Solve $U_\eta^\top U_\eta\, \mathbf{m}_\eta = \tilde H^\top \tilde{\mathbf{y}}$ by forward-then-back substitution.
4. Draw $\mathbf{z} \sim N(0, I_p)$ and set $\boldsymbol{\eta} = \mathbf{m}_\eta + U_\eta^{-1} \mathbf{z}$.

**Code reference:** `gibbs_interaction_sampler()` lines 210–227 in `gibbs_interaction.R`. The logic is identical to the main-effects-only sampler in `gibbs_stage_c_full.R` — only the dimensions of $H$, $Q_0$, and $\boldsymbol{\eta}$ grow.

### Step 33. Update for $\tau^2_{s,j}$ (main effect) — Conjugate IG

With $\beta_j \mid \tau^2_{s,j} \sim N(0, \tau^2_{s,j}\, K_j^{-})$ and $\tau^2_{s,j} \sim \text{IG}(a_s, b_s)$, the full conditional is inverse gamma:

$$
\tau^2_{s,j} \mid \beta_j \;\sim\; \text{IG}\!\left(a_s + \frac{\text{rank}(K_j^{(\beta)})}{2},\;\; b_s + \frac{\beta_j^\top K_j^{(\beta)}\, \beta_j}{2}\right)
$$

With $\text{rank}(K_j^{(\beta)}) = M_j - 2$.

**Code reference:** `gibbs_interaction.R` lines 256–265.

### Step 34. Update for $\tau^2_{s,uv}$ (interaction) — Conjugate IG

Identical form, different matrix and different rank:

$$
\tau^2_{s,uv} \mid \beta_{uv} \;\sim\; \text{IG}\!\left(a_s + \frac{r_{uv}}{2},\;\; b_s + \frac{\beta_{uv}^\top K_{uv}^{(\beta)}\, \beta_{uv}}{2}\right)
$$

where $r_{uv} = \text{rank}(K_{uv}^{(\beta)}) = (M_u - 1)(M_v - 1) - 1$ in our implementation (the nullity-1 result from Step 23 of the construction document, verified numerically).

**Code reference:** `gibbs_interaction.R` lines 270–280. Note line 276: `rank_k <- nrow(K_k) - 1`.

### Step 35. Update for $\sigma^2, \tau^2, \rho$ — Metropolis-Hastings on Collapsed Likelihood

These three parameters enter the likelihood nonlinearly through $\Sigma = \sigma^2 R + \tau^2 I$ and through $R$'s dependence on $\rho$. No closed-form full conditional — use random-walk MH on the log scale.

For each parameter $\phi \in \{\sigma^2, \tau^2, \rho\}$:

1. Propose $\phi' = \phi \cdot \exp(w)$ where $w \sim N(0, s_\phi^2)$ (log-scale random walk).
2. Compute the acceptance log-ratio:

$$
\log \alpha = \log p(\mathbf{y} \mid \phi', \cdot) + \log p(\phi') - \log p(\mathbf{y} \mid \phi, \cdot) - \log p(\phi)
$$

3. Accept $\phi'$ with probability $\min(1, \alpha)$, otherwise keep $\phi$.

For $\rho$ specifically, when the proposal is accepted we must recompute $R$ (which changes every entry since Matérn correlation depends on $\rho$).

**Proposal SDs (tuning):** $s_{\log \sigma^2} = s_{\log \tau^2} = 0.3$, $s_{\log \rho} = 0.2$.

**Code reference:** `gibbs_interaction.R` lines 235–251.

### Step 36. Derived Posterior Mean of $\mathbf{b}$

After sampling $\boldsymbol{\eta}, \sigma^2, \tau^2, \rho$, the conditional mean of $\mathbf{b}$ given these is a derived quantity:

$$
\mathbb{E}[\mathbf{b} \mid \mathbf{y}, \boldsymbol{\eta}, \sigma^2, \tau^2, \rho] = \sigma^2 R\, \Sigma^{-1}\, (\mathbf{y} - H\boldsymbol{\eta})
$$

This is stored at each retained iteration so we have posterior samples of the spatial random effect for diagnostics.

**Code reference:** `gibbs_interaction.R` lines 161–165, called on line 285.

### Step 37. Complete Gibbs Sweep

Each iteration does the following in order:

| Step | Parameter | Type | Where in code |
|------|-----------|------|---------------|
| 1 | $\boldsymbol{\eta}$ | MVN (conjugate) | Lines 210–227 |
| 2 | $\sigma^2$ | MH on log scale | Lines 234–238 |
| 3 | $\tau^2$ | MH on log scale | Lines 240–244 |
| 4 | $\rho$ | MH on log scale (recompute $R$) | Lines 246–251 |
| 5a | $\tau^2_{s,j}$, $j = 1, \ldots, p$ | IG (conjugate) | Lines 256–265 |
| 5b | $\tau^2_{s,uv}$, all pairs | IG (conjugate) | Lines 270–280 |
| 6 | $\mathbf{b}$ (derived) | Deterministic | Line 285 |

The sweep's structure is identical to the main-effects-only sampler; Step 5b is the only genuinely new step.

---

## Part VII: The LS Interaction Construction in Code

This section maps the construction document's equations to specific functions in `ls_interaction.R`.

### Step 38. Building the 1D Identified RW2 Penalty $K_j^{(\beta)}$

Construction equation: $K_j^{(\beta)} = T_j^\top K_j^{(\theta)} T_j$ where $K_j^{(\theta)} = (D_2^{(j)})^\top D_2^{(j)}$.

Code (from `run_interaction_2x2_v2.R` lines 132–138):

```r
K_raw <- build_rw2_penalty_1d(M)                      # K_j^{(theta)}
K_main_list[[j]] <- t(objs[[j]]$T) %*% K_raw %*% objs[[j]]$T   # K_j^{(beta)}
```

`build_rw2_penalty_1d()` in `ls_interaction.R` lines 31–40 builds $D_2$ explicitly and returns $D_2^\top D_2$.

### Step 39. Building the 2D Penalty $K_{uv}^{(\beta)}$

Construction equation: $K_{uv}^{(\theta)} = K_u^{(\theta)} \otimes I_{M_v} + I_{M_u} \otimes K_v^{(\theta)}$, then $K_{uv}^{(\beta)} = (T_u \otimes T_v)^\top K_{uv}^{(\theta)} (T_u \otimes T_v)$.

Code (`build_rw2_2d_penalty()` in `ls_interaction.R` lines 65–77):

```r
K_uv_raw <- kronecker(K_u_raw, diag(M_v)) + kronecker(diag(M_u), K_v_raw)
T_uv <- kronecker(T_u, T_v)
K_uv_beta <- t(T_uv) %*% K_uv_raw %*% T_uv
```

This is the "slow but obvious" route — it forms the large $T_u \otimes T_v$ explicitly and does the triple multiplication. For $M = 8$ this is a $49 \times 49$ matrix, no issue. For larger $M$ you'd use the algebraic shortcut $K_{uv}^{(\beta)} = K_u^{(\beta)} \otimes G_v + G_u \otimes K_v^{(\beta)}$ from Step 21, but the current code takes the direct route for clarity.

### Step 40. Building the Interaction Design Matrix $W_{uv}$

Construction shortcut: $W_{uv} = W_u \odot W_v$ (Khatri-Rao, from Step 15 of the construction document).

Code (`khatri_rao_rowwise_R()` in `ls_interaction.R` lines 91–103):

```r
W_u_rep <- W_u[, rep(seq_len(d_u), each  = d_v), drop = FALSE]
W_v_rep <- W_v[, rep(seq_len(d_v), times = d_u), drop = FALSE]
W_u_rep * W_v_rep   # element-wise product
```

Verifying this is the Khatri-Rao: row $i$ of the result has entry at column $(a-1)d_v + b$ equal to $W_u[i,a] \cdot W_v[i,b]$ — exactly the Kronecker of row $i$ of $W_u$ with row $i$ of $W_v$. Matches Step 15.

An Rcpp/Eigen version `khatri_rao_cpp()` exists in `ls_interaction_core.cpp` for speed on larger problems.

### Step 41. Assembling Everything: `ls_build_interaction()`

The function `ls_build_interaction(obj_u, obj_v)` in `ls_interaction.R` lines 125–174 is the main entry point. Given the two 1D LS objects (from `ls_build_one_full`), it returns a list containing:

- `W_uv`: the $n \times (M_u-1)(M_v-1)$ identified interaction design matrix (Step 15 of construction)
- `K_uv`: the $(M_u-1)(M_v-1) \times (M_u-1)(M_v-1)$ identified 2D penalty (Step 21 of construction)
- `T_uv`: the $M_u M_v \times (M_u-1)(M_v-1)$ combined contrast (kept for reference)
- `recipe`: stored so we can build $W_{uv}$ at new $(x_u, x_v)$ points for prediction

### Step 42. Sanity Checks (`ls_interaction_tests`)

`ls_interaction.R` lines 280–332 runs a battery of checks:

1. Dimension check: `W_uv` should be $n \times (M-1)^2$
2. Dimension check: `K_uv` should be $(M-1)^2 \times (M-1)^2$
3. PSD check: all eigenvalues of `K_uv` should be $\geq -10^{-8}$
4. Rank check: `qr(K_uv)$rank` should equal $(M-1)^2 - 1$ (confirming Step 23 of the construction numerically)
5. Contrast check: $\mathbf{1}^\top T_u = \mathbf{0}^\top$ (the key 1D identifiability property)
6. Prediction check: design at new $(x_u, x_v)$ has the correct dimension

Running `ls_interaction_tests(n=200, M=8)` passes all checks.

---

## Part VIII: The 2×2 Experimental Design

This section documents what we actually ran — the setup that validates the interaction construction end-to-end.

### Step 43. Data-Generating Process (DGP)

True functions (centered to have approximate mean zero on $[0,1]$):

$$
f_1(x) = 2.0\,(\sin(\pi x) - 2/\pi)
$$

$$
f_2(x) = 1.5\,(e^{x - 0.5} - C_2), \qquad C_2 = e^{0.5} - e^{-0.5}
$$

$$
f_3(x) = 0.7\,(x^2 - 1/3)
$$

$$
f_4(x) = 0.5\,\sin(2\pi x)
$$

**True interaction** (only $f_{12}$ nonzero — all other pairs are zero):

$$
f_{1,2}(x_1, x_2) = c_{\text{int}}\, (\sin(\pi x_1) - 2/\pi)\,(e^{x_2 - 0.5} - C_2)
$$

**Spatial + noise:**

$$
b(s) \sim \text{Matérn GP}(\rho = 0.3,\, \nu = 1.5),\quad \sigma^2 = 0.5
$$

$$
\varepsilon \sim N(0, \tau^2),\quad \tau^2 = 0.25
$$

Covariates $X_1, \ldots, X_p \sim \text{Uniform}(0,1)$ independently; locations uniform on $[0,1]^2$.

**Code reference:** `simulate_unified()` in `run_interaction_2x2_v2.R` lines 78–111.

### Step 44. The 2×2 Factorial Structure

For each setting $(p, \text{phase})$, we run four cells in a $2 \times 2$ design:

**Truth axis** (vertical):

- **NULL**: $c_{\text{int}} = 0$ (no true interaction — $f_{12}$ is actually zero)
- **INT**: $c_{\text{int}} = 1.5$ (strong true interaction)

**Model axis** (horizontal):

- **$M_0$**: additive main-effects-only model (no interaction terms fit)
- **$M_1$**: model with interaction terms fit (either just $X_1 \times X_2$ in phase 1, or all pairwise in phase 2)

| | $M_0$ (no int) | $M_1$ (with int) |
|---|---|---|
| **NULL** ($c=0$) | Cell A: correctly specified | Cell B: overfitting check |
| **INT** ($c=1.5$) | Cell C: underfitting check | Cell D: correctly specified |

**What each cell tests:**

- **Cell A** (NULL + $M_0$): Baseline — the main-effects-only sampler under its correct assumption. Confirms the existing machinery still works.
- **Cell B** (NULL + $M_1$): When the truth has no interaction but we fit one anyway, does the interaction term shrink toward zero? Does its 95% credible band cover zero?
- **Cell C** (INT + $M_0$): When the truth has an interaction but we ignore it, where does the missing signal go? Usually it leaks into main effects (RMSE inflates) or into $\sigma^2$ (the spatial nugget absorbs it).
- **Cell D** (INT + $M_1$): The target scenario. Does the fitted interaction surface recover the true surface? Does the credible band cover the truth?

**Code reference:** `run_2x2()` in `run_interaction_2x2_v2.R`.

### Step 45. Two-Phase Rollout

For each $p \in \{2, 3, 4\}$:

- **Phase 1**: $M_1$ includes only the $X_1 \times X_2$ interaction. Tests the simplest case — the interaction model exactly matches the single true nonzero interaction.
- **Phase 2**: $M_1$ includes all $\binom{p}{2}$ pairwise interactions. Tests whether the "extra" interaction terms (which are truly zero) correctly shrink to zero while the $X_1 \times X_2$ term recovers the truth.

This gives $3 \times 2 \times 4 = 24$ total cells across all $(p, \text{phase})$ combinations.

### Step 46. Metrics Computed Per Cell

For each fitted cell, `summarise_cell_unified()` (lines 215–341) computes:

**Main effect recovery:**

- `RMSE_f_j` = $\sqrt{\text{mean}((f_j(\hat{\boldsymbol{\beta}}_j) - f_j^{\text{true}})^2)}$ for $j = 1, \ldots, p$

**Interaction recovery** (Cell D only):

- `RMSE_f12` = $\sqrt{\text{mean}((\hat{f}_{12} - f_{12}^{\text{true}})^2)}$ at the training points
- `cov_f12` = fraction of training points where the 95% posterior band contains the truth

**Variance parameter recovery:**

- Posterior mean of $\sigma^2$, $\tau^2$, $\rho$
- Trace plots, posterior densities, autocorrelation, effective sample size

**Model comparison:**

- WAIC (to compare $M_0$ vs $M_1$; lower = better predictive performance)

**Posterior interaction surface:**

- Evaluated on a $50 \times 50$ regular grid over $[0,1]^2$
- Posterior mean surface
- Pointwise 95% credible band (and its width, for mapping uncertainty)

### Step 47. What We're Looking For — Expected Results

If the interaction construction is correct, we expect:

1. **Cell A**: matches existing main-effects results (regression test).
2. **Cell B**: `cov_f12` ≈ 0.95 with $\hat{f}_{12}$ ≈ 0. The credible band covers zero because the posterior correctly identifies no interaction signal. WAIC should slightly favor $M_0$ (no interaction) due to Occam's factor — the penalty for unnecessary complexity.
3. **Cell C**: `RMSE_f_j` inflated compared to Cell A; some estimated main effects absorb structure from the true interaction. $\sigma^2$ may be biased upward.
4. **Cell D**: `RMSE_f12` small, `cov_f12` ≈ 0.95, estimated interaction surface visually close to truth. WAIC should favor $M_1$.
5. **Across $p$**: As $p$ grows, the interaction-blocks-to-observations ratio grows. By phase 2 with $p = 4$, the model has $\binom{4}{2}(M-1)^2 = 6(M-1)^2$ interaction parameters; with $M = 8$ that's 294 interaction parameters on $n = 300$ observations — quite aggressive. Regularization from the RW2 prior matters a lot here.

### Step 48. Outputs

Per $(p, \text{phase})$ combination:

- `interaction_2x2_v2_p{p}_phase{phase}_marginals.pdf`: fitted main effects with truth overlay, for all four cells
- `interaction_2x2_v2_p{p}_phase{phase}_surfaces.pdf`: fitted $f_{12}$ surface image plots + uncertainty + residual-from-truth, for cells B and D
- `interaction_2x2_v2_p{p}_phase{phase}_diagnostics.pdf`: trace plots, posterior densities, ACFs for all variance parameters

Summary across all 24 cells:

- `interaction_2x2_v2_summary.csv`: one row per cell with RMSE, coverage, WAIC, variance estimates, timing
- `interaction_2x2_v2_results.rds`: full R object with chains for reanalysis

---

## Part IX: Summary of What's Inherited vs. New

| Component | Implementation status | File/location |
|-----------|---------------------|---------------|
| 1D LS basis (Part I of construction) | Existing (unchanged) | `ls_basis.R` |
| Collapsed likelihood, $\boldsymbol{\eta}$ MVN update | Existing (unchanged structure) | `gibbs_stage_c_full.R` |
| MH updates for $\sigma^2, \tau^2, \rho$ | Existing (unchanged) | `gibbs_stage_c_full.R` |
| IG update for main-effect $\tau^2_{s,j}$ | Existing (unchanged) | `gibbs_stage_c_full.R` |
| Spatial random effect marginalization | Existing (unchanged) | `spatial_utils.R` |
| **Khatri-Rao $W_{uv} = W_u \odot W_v$** | **New** | `ls_interaction.R` |
| **2D Kronecker-sum penalty $K_{uv}^{(\beta)}$** | **New** | `ls_interaction.R` |
| **Interaction contrast $T_u \otimes T_v$** | **New** | `ls_interaction.R` |
| **IG update for interaction $\tau^2_{s,uv}$** | **New** | `gibbs_interaction.R` |
| **Block prior precision with interaction blocks** | **New** | `gibbs_interaction.R` |
| **2×2 experiment driver** | **New** | `run_interaction_2x2_v2.R` |
| **Surface grid prediction at new $(x_u, x_v)$** | **New** | `ls_interaction_design_new()` |
| Proof-of-concept simulator | **New** | `run_interaction_poc()` |
| Centered Sim4 DGP with interaction | **New** | `run_interaction_2x2_v2.R` |

---

## Notation Additions (beyond construction document)

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| $\tau^2_{s,j}$ | Main effect smoothing variance (covariate $j$) | $\sim 1$ at init |
| $\tau^2_{s,uv}$ | Interaction smoothing variance (pair $u,v$) | $\sim 1$ at init |
| $a_s, b_s$ | IG hyperparameters for smoothing variances | $a_s = 1, b_s = 0.005$ |
| $a_\sigma, b_\sigma$ | IG hyperparameters for spatial variance | $a_\sigma = 2, b_\sigma = 1$ |
| $a_\tau, b_\tau$ | IG hyperparameters for nugget | $a_\tau = 2, b_\tau = 0.3$ |
| $\log \rho_\mu, \log \rho_\sigma$ | Log-normal prior on Matérn range | $-1.6, 1.0$ |
| $\kappa^2$ | Vague intercept prior variance | $10^6$ |
| $\epsilon_{\text{ridge}}$ | Tiny ridge on $Q_0$ for numerical stability | $10^{-6}$ |
| $c_{\text{int}}$ | Interaction strength in DGP | $0$ or $1.5$ |
| $p$ | Number of covariates | $\{2, 3, 4\}$ |
| $K$ (implicit) | Number of interaction blocks $= \binom{p}{2}$ | $\{1, 3, 6\}$ |
