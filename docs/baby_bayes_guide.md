# Baby Bayes: Stage A Guide

## Where You Are Now

Your current model (from `fit_spatial_reml.R`) is:

$$
y = \mu \mathbf{1} + W\beta + b + \varepsilon
$$

where:
- $W$ is the LS-basis design (after identifiability reparameterization $\theta_j = T_j \beta_j$, sum-to-zero)
- $b \sim N(0, \sigma^2 R(\rho, \nu))$ is the spatial random effect
- $\varepsilon \sim N(0, \tau^2 I_n)$ is the nugget

Stacking the intercept and spline into one design matrix $H = [\mathbf{1} \mid W]$ and parameter $\eta = (\mu, \beta^\top)^\top$:

$$
y \sim N(H\eta, \; \Sigma), \quad \Sigma = \sigma^2 R(\rho, \nu) + \tau^2 I_n
$$

Your REML code estimates $(\rho, \sigma^2, \tau^2)$ by optimization, then computes $\hat{\eta}$ by GLS:

$$
\hat{\eta} = (H^\top \Sigma^{-1} H)^{-1} H^\top \Sigma^{-1} y
$$

This is already **Empirical Bayes** — you're plugging in point estimates of variance parameters and computing the conditional mode of $\eta$.

---

## What Baby Bayes Does (and Doesn't Do)

**Baby Bayes = your current GLS, rewritten as a Bayesian posterior with explicit priors.**

- Fix $(\rho, \sigma^2, \tau^2)$ at truth (in simulation) or REML estimates
- Put a vague Gaussian prior on $\eta$: $\eta \sim N(0, \kappa^2 I_p)$ with $\kappa^2$ very large (e.g., $10^6$)
- Compute the closed-form posterior — no MCMC needed

**What you gain:** Bayesian bookkeeping (prior + likelihood → posterior), credible intervals, and a sanity check that everything is set up correctly before adding MCMC complexity.

**What you don't gain (yet):** uncertainty about $\sigma^2, \tau^2, \rho$. Those stay fixed.

---

## The Math (All Closed Form)

### Likelihood

Given fixed $\Sigma$:

$$
p(y \mid \eta) \propto \exp\!\Big(-\tfrac{1}{2}(y - H\eta)^\top \Sigma^{-1}(y - H\eta)\Big)
$$

### Prior

$$
p(\eta) \propto \exp\!\Big(-\tfrac{1}{2}\eta^\top \Lambda_0^{-1} \eta\Big)
$$

where $\Lambda_0 = \kappa^2 I_p$ (prior covariance). For Bundle A, $\kappa^2 = 10^6$.

### Posterior

The posterior $p(\eta \mid y)$ is Gaussian with:

$$
\boxed{
\begin{aligned}
\text{Posterior precision:} \quad Q &= H^\top \Sigma^{-1} H + \Lambda_0^{-1} \\[6pt]
\text{Posterior mean:} \quad \bar{\eta} &= Q^{-1} H^\top \Sigma^{-1} y \\[6pt]
\text{Posterior covariance:} \quad V &= Q^{-1}
\end{aligned}
}
$$

When $\kappa^2 \to \infty$, $\Lambda_0^{-1} \to 0$, so $Q \to H^\top \Sigma^{-1} H$ and $\bar\eta \to \hat\eta_{\text{GLS}}$. This is exactly your current code.

### Credible Intervals

For any component $\eta_k$:

$$
\eta_k \mid y \;\sim\; N(\bar{\eta}_k, \; V_{kk})
$$

95% credible interval: $\bar{\eta}_k \pm 1.96\sqrt{V_{kk}}$

### For Marginal Curves

To get a credible band for $f_j(x_g) = w_j(x_g)^\top \beta_j$ at a grid point $x_g$:

Let $c_g$ be the row of the full design $H$ corresponding to grid point $x_g$ (with zeros for other covariates). Then:

$$
f_j(x_g) \mid y \;\sim\; N(c_g^\top \bar\eta, \; c_g^\top V c_g)
$$

This gives you pointwise credible bands on each marginal curve.

---

## Mapping to Your Code

Here's how each math object maps to your existing R code:

| Math | Your Code | Where |
|------|-----------|-------|
| $H$ | `X_fix = cbind(1, W)` | `fit_ls_spatial()`, line 115 |
| $\eta$ | `fit$beta` (intercept + spline coefficients) | `spatial_reml_fit()`, line 82 |
| $\Sigma$ | `sigma2 * R + tau2 * I` = `sigma2 * Sigma0` | `spatial_reml_fit()`, line 75 |
| $\Sigma^{-1} y$ | You already compute this via Cholesky: `forwardsolve(t(U0), y)` | line 78 |
| $H^\top \Sigma^{-1} H$ | `crossprod(X_t)` where `X_t = L^{-1} H` | line 79, 81 |
| $\hat\eta_{\text{GLS}}$ | `solve(XtX, crossprod(X_t, y_t))` | line 82 |

The key insight: **your code already computes the GLS estimate, which IS the posterior mean when the prior is flat.** Baby Bayes just makes this explicit and adds the posterior covariance.

---

## The Baby Bayes Function (R Code)

```r
# ============================================================
# baby_bayes_fit.R
# Stage A: Closed-form Bayesian posterior for eta = (mu, beta)
# with FIXED (rho, sigma2, tau2)
# ============================================================

baby_bayes_fit <- function(y, X_fix, Sigma, kappa2 = 1e6, jitter = 1e-8) {
  # y      : n-vector of responses
  # X_fix  : n x p design matrix (= cbind(1, W))
  # Sigma  : n x n covariance matrix (= sigma2 * R + tau2 * I)
  # kappa2 : prior variance for eta (large = vague prior)
  
  n <- length(y)
  p <- ncol(X_fix)
  
  # --- Cholesky of Sigma for stable inversion ---
  L <- chol(Sigma + diag(jitter, n))    # upper triangular: Sigma = t(L) %*% L
  
  # Whiten: multiply by Sigma^{-1/2}
  y_w <- forwardsolve(t(L), y)           # L^{-T} y
  H_w <- forwardsolve(t(L), X_fix)       # L^{-T} H
  
  # --- Prior precision ---
  Lambda0_inv <- diag(1 / kappa2, p)     # (kappa^2 I)^{-1}
  
  # --- Posterior ---
  # Q = H' Sigma^{-1} H + Lambda0^{-1}
  HtSinvH <- crossprod(H_w)              # H' Sigma^{-1} H
  Q <- HtSinvH + Lambda0_inv
  
  # Posterior covariance V = Q^{-1}
  U_Q <- chol(Q)
  V <- chol2inv(U_Q)                     # stable inverse via Cholesky
  
  # Posterior mean: eta_bar = V * H' Sigma^{-1} y
  HtSinv_y <- crossprod(H_w, y_w)        # H' Sigma^{-1} y
  eta_bar <- as.vector(V %*% HtSinv_y)
  
  # --- Standard errors and credible intervals ---
  eta_se <- sqrt(diag(V))
  ci_lower <- eta_bar - 1.96 * eta_se
  ci_upper <- eta_bar + 1.96 * eta_se
  
  list(
    eta_bar  = eta_bar,       # posterior mean (should match GLS)
    V        = V,             # posterior covariance
    eta_se   = eta_se,        # posterior std dev
    ci_lower = ci_lower,      # 95% credible interval lower
    ci_upper = ci_upper,      # 95% credible interval upper
    Q        = Q,             # posterior precision
    HtSinvH  = HtSinvH       # useful for diagnostics
  )
}
```

---

## How to Use It With Your Existing Pipeline

```r
source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")

# --- Step 1: Simulate data (same as your Sim2) ---
# ... your simulate_sim2() call ...

# --- Step 2: Fit your current REML model ---
obj <- fit_ls_spatial(y = y, X_raw = X_raw, coords = coords,
                      M_vec = rep(6, 4), nu = 1.5)

# --- Step 3: Build Sigma from fitted (or true) parameters ---
# Option A: use REML estimates
sigma2 <- obj$fit$sigma2
tau2   <- obj$fit$tau2
R      <- obj$fit$R
Sigma  <- sigma2 * R + tau2 * diag(length(y))

# Option B (recommended for sanity check): use TRUE values
# sigma2_true <- 0.8; tau2_true <- 0.15; rho_true <- 0.2
# D <- pairdist(coords)
# R_true <- matern_cor(D, rho = rho_true, nu = 1.5)
# Sigma <- sigma2_true * R_true + tau2_true * diag(length(y))

# --- Step 4: Run Baby Bayes ---
bb <- baby_bayes_fit(y = y, X_fix = obj$X_fix, Sigma = Sigma)

# --- Step 5: Sanity check — posterior mean should match GLS ---
cat("Max |eta_bar - eta_GLS| =", max(abs(bb$eta_bar - obj$fit$beta)), "\n")
# Should be ~1e-6 or smaller when kappa2 is large
```

---

## Getting Credible Bands on Marginal Curves

```r
# For covariate j, on a grid x_g:
bayes_marginal_band <- function(bb, obj, j, x_grid, clip = TRUE) {
  # bb  : output from baby_bayes_fit()
  # obj : output from fit_ls_spatial()
  # j   : covariate index (1..p)
  
  # Build W_grid for covariate j
  W_grid_j <- obj$des$objs[[j]]$design_new(x_grid, type = "W", clip = clip)
  # n_grid x (M_j - 1)
  
  # Which columns in the full eta correspond to covariate j?
  col_idx <- 1 + obj$des$col_map[[j]]   # +1 for intercept
  
  # Posterior mean of f_j(x_grid)
  beta_j_bar <- bb$eta_bar[col_idx]
  f_hat <- as.vector(W_grid_j %*% beta_j_bar)
  
  # Posterior variance of f_j(x_grid) at each grid point
  V_j <- bb$V[col_idx, col_idx]         # sub-block of posterior covariance
  f_var <- rowSums((W_grid_j %*% V_j) * W_grid_j)  # diag(W V W')
  f_se  <- sqrt(f_var)
  
  data.frame(
    x      = x_grid,
    f_hat  = f_hat,
    f_se   = f_se,
    lower  = f_hat - 1.96 * f_se,
    upper  = f_hat + 1.96 * f_se
  )
}

# Example usage:
x_grid <- seq(0, 1, length.out = 101)
band_j1 <- bayes_marginal_band(bb, obj, j = 1, x_grid = x_grid)

# Plot with credible band
plot(band_j1$x, band_j1$f_hat, type = "l", ylim = range(band_j1$lower, band_j1$upper),
     main = "X1 marginal with 95% credible band", xlab = "x", ylab = "f1(x)")
polygon(c(band_j1$x, rev(band_j1$x)),
        c(band_j1$lower, rev(band_j1$upper)),
        col = rgb(0, 0, 1, 0.15), border = NA)
lines(band_j1$x, band_j1$f_hat, lwd = 2)
# Add truth:
lines(x_grid, 2*sin(pi*x_grid) - mean(2*sin(pi*runif(400))), lty = 2, col = "red")
```

---

## What to Check (Your Sanity Checklist)

### Check 1: Posterior mean ≈ GLS estimate

```r
max(abs(bb$eta_bar - obj$fit$beta))
# Should be < 1e-5 when kappa2 = 1e6
```

If this is large, something is wrong in the matrix algebra. Debug by checking that `Sigma` matches what your REML code uses internally.

### Check 2: Credible intervals are reasonable

```r
# Coverage: does the true eta fall inside the 95% CI?
# (Only possible in simulation where you know the truth)
eta_true <- c(mu_true, beta_true)  # you'd need to extract these from your sim
mean(eta_true >= bb$ci_lower & eta_true <= bb$ci_upper)
# Should be roughly 0.95
```

### Check 3: Marginal bands contain the truth

Plot $f_j(x)$ with the credible band from `bayes_marginal_band()` and overlay the true function. The true curve should mostly fall inside the band.

### Check 4: Prior doesn't matter (yet)

Run with `kappa2 = 1e6` and `kappa2 = 1e3`. Results should be nearly identical (because the data dominates the prior at this stage). If they differ substantially, your data is too weak for the model or there's a coding bug.

---

## What Comes Next (Stage B Preview)

Once Baby Bayes is working, Stage B adds:

- **Priors on $\sigma^2$ and $\tau^2$** (e.g., inverse-gamma or half-t)
- **Gibbs sampler** to alternate between:
  - Draw $\eta \mid y, \sigma^2, \tau^2$ (this is just your Baby Bayes posterior — you already have it!)
  - Draw $\sigma^2 \mid y, \eta, \tau^2$ (conditionally conjugate if using inverse-gamma prior)
  - Draw $\tau^2 \mid y, \eta, \sigma^2$ (same)
- $\rho$ still fixed

The key point: **the Baby Bayes posterior for $\eta$ becomes one step of the Gibbs sampler in Stage B.** That's why getting Stage A right matters — you'll reuse this code directly.

---

## Summary

| What | Baby Bayes (Stage A) |
|------|---------------------|
| Parameters with priors | $\eta = (\mu, \beta)$ only |
| Fixed | $\rho, \sigma^2, \tau^2$ |
| Prior | $\eta \sim N(0, 10^6 I)$ |
| Posterior | Closed-form Gaussian |
| MCMC needed? | No |
| Result should match | Your current GLS $\hat\eta$ |
| New output | Credible intervals + bands |
| Test on | Sim1 and Sim2 |
| Paper reference | Crainiceanu et al. (2005), Section 2 |
