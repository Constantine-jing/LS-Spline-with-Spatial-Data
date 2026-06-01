# LS-Spline Tensor Product Interaction: Mathematical Construction

## 1. Recap: 1D LS-Spline Basis for One Covariate

### 1.1 Setup

Fix covariate $j$ with $M_j \geq 4$ knots:

$$
\tau_{1j} < \tau_{2j} < \cdots < \tau_{M_j, j}
$$

Knot spacings: $h_{mj} = \tau_{mj} - \tau_{m-1,j}$ for $m = 2, \ldots, M_j$.

### 1.2 LS Basis Functions

The LS construction provides two families of compactly supported cubic spline basis functions:

- $\{\Phi_{mj}(x)\}_{m=1}^{M_j}$ — the **value basis** (cardinal for function values)
- $\{\Psi_{mj}(x)\}_{m=1}^{M_j}$ — the **slope basis** (cardinal for derivatives)

**Cardinal properties:**

$$
\Phi_{mj}(\tau_{\ell j}) = \delta_{m\ell}, \qquad
\Psi'_{mj}(\tau_{\ell j}) = \delta_{m\ell}
$$

$$
\Phi'_{\ell j}(\tau_{mj}) = 0 \;\;\forall \ell, \qquad
\Psi_{\ell j}(\tau_{mj}) = 0 \;\;\forall \ell
$$

Any natural cubic spline with these knots can be written as:

$$
f_j(x) = \sum_{m=1}^{M_j} \theta_{mj}\, \Phi_{mj}(x) + \sum_{m=1}^{M_j} \gamma_{mj}\, \Psi_{mj}(x)
$$

where:
- $\theta_{mj} = f_j(\tau_{mj})$ are the **knot ordinates** (function values at knots)
- $\gamma_{mj} = f'_j(\tau_{mj})$ are the **knot slopes** (derivative values at knots)

### 1.3 Eliminating Slopes: The A and C Matrices

The natural spline conditions $f''_j(\tau_{1j}) = 0$, $f''_j(\tau_{M_j,j}) = 0$, plus continuity of $f''_j$ at interior knots, impose the linear relation:

$$
C_j\, \theta_j = A_j\, \gamma_j
\qquad \Longrightarrow \qquad
\gamma_j = A_j^{-1} C_j\, \theta_j
$$

where $A_j \in \mathbb{R}^{M_j \times M_j}$ is tridiagonal and $C_j \in \mathbb{R}^{M_j \times M_j}$ is sparse (both defined by knot spacings and $\omega$-weights; see Appendix 5 of the dissertation draft).

### 1.4 Reduced Basis and Design Matrix

Substituting $\gamma_j = A_j^{-1} C_j\, \theta_j$ into the spline representation eliminates the slopes:

$$
f_j(x) = z_j(x)^\top \theta_j
$$

where the **reduced LS basis vector** is:

$$
z_j(x)^\top = \Phi_j(x)^\top + \Psi_j(x)^\top A_j^{-1} C_j,
\qquad z_j(x) \in \mathbb{R}^{M_j}
$$

**Key property:** $z_{m'j}(\tau_{mj}) = \delta_{m'm}$ (cardinal), so $f_j(\tau_{mj}) = \theta_{mj}$.

The **design matrix** for $n$ observations is:

$$
Z_j = \begin{pmatrix} z_j(x_{1j})^\top \\ \vdots \\ z_j(x_{nj})^\top \end{pmatrix} \in \mathbb{R}^{n \times M_j}
$$

so that the vector of fitted values is $\mathbf{f}_j = Z_j\, \theta_j$.

### 1.5 Identifiability: The Contrast Matrix T

The model includes an intercept $\mu$ and each $Z_j$ contains constants in its column space. To separate $f_j$ from $\mu$, we impose **sum-to-zero** on the knot ordinates:

$$
\mathbf{1}^\top \theta_j = \sum_{m=1}^{M_j} \theta_{mj} = 0
$$

This means $\theta_{1j} = -(\theta_{2j} + \cdots + \theta_{M_j,j})$.

Define the **contrast matrix** $T_j \in \mathbb{R}^{M_j \times (M_j - 1)}$:

$$
T_j = \begin{pmatrix} -1 & -1 & \cdots & -1 \\ 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & 1 \end{pmatrix}
$$

Then $\theta_j = T_j\, \beta_j$ enforces $\mathbf{1}^\top \theta_j = 0$ automatically.

**Critical property of $T_j$:**

$$
\mathbf{1}^\top T_j = \mathbf{0}^\top
$$

*Proof:* The first row of $T_j$ is all $-1$'s, and rows $2, \ldots, M_j$ form $I_{M_j - 1}$. So $\mathbf{1}^\top T_j = (-1, -1, \ldots, -1) + (1, 0, \ldots, 0) + (0, 1, \ldots, 0) + \cdots + (0, 0, \ldots, 1) = \mathbf{0}^\top$. $\square$

The **identified design matrix** is:

$$
W_j = Z_j\, T_j \in \mathbb{R}^{n \times (M_j - 1)}
$$

and the **identified coefficients** are $\beta_j \in \mathbb{R}^{M_j - 1}$, so $\mathbf{f}_j = W_j\, \beta_j$.

### 1.6 The RW2 Penalty

The second-order random walk penalty on $\beta_j$ is:

$$
\beta_j^\top K_j^{(\beta)}\, \beta_j,
\qquad K_j^{(\beta)} = D_2^\top D_2
$$

where $D_2 \in \mathbb{R}^{(d_j - 2) \times d_j}$ is the second-difference matrix ($d_j = M_j - 1$):

$$
D_2 = \begin{pmatrix}
1 & -2 & 1 & 0 & \cdots \\
0 & 1 & -2 & 1 & \cdots \\
  &   & \ddots & \ddots & \ddots
\end{pmatrix}
$$

Equivalently, this can be expressed on the unidentified $\theta_j$:

$$
K_j^{(\theta)} = (D_2^{(\theta)})^\top D_2^{(\theta)}
$$

where $D_2^{(\theta)} \in \mathbb{R}^{(M_j - 2) \times M_j}$, and $K_j^{(\beta)} = T_j^\top K_j^{(\theta)} T_j$.


---

## 2. The 2D Tensor Product Interaction

Consider two covariates $u$ and $v$ with $M_u$ and $M_v$ knots respectively. We construct the pairwise interaction surface $f_{u,v}(x_u, x_v)$.

### 2.1 Tensor Product of the Reduced Bases

Define the interaction surface as:

$$
f_{u,v}(x_u, x_v) = \sum_{m=1}^{M_u} \sum_{l=1}^{M_v} \theta_{ml}\; z_{um}(x_u)\; z_{vl}(x_v)
$$

where $z_{um}(x_u)$ is the $m$-th element of $z_u(x_u)$ and $z_{vl}(x_v)$ is the $l$-th element of $z_v(x_v)$.

**Knot-ordinate property in 2D.** Evaluate at the grid point $(\tau_{mu}, \tau_{lv})$:

$$
f_{u,v}(\tau_{mu}, \tau_{lv})
= \sum_{m'=1}^{M_u} \sum_{l'=1}^{M_v} \theta_{m'l'}\; z_{um'}(\tau_{mu})\; z_{vl'}(\tau_{lv})
$$

By the 1D cardinal property, $z_{um'}(\tau_{mu}) = \delta_{m'm}$ and $z_{vl'}(\tau_{lv}) = \delta_{l'l}$. Therefore:

$$
\boxed{f_{u,v}(\tau_{mu}, \tau_{lv}) = \theta_{ml}}
$$

**Each coefficient $\theta_{ml}$ is the surface height at the grid point $(\tau_{mu}, \tau_{lv})$.**

This is the **2D knot-ordinate property**, unique to the LS basis. B-spline tensor product coefficients do not have this interpretation.

### 2.2 Vectorization and Design Matrix

Stack the $M_u \times M_v$ coefficient matrix into a vector using the convention that the $v$-index varies fastest:

$$
\theta_{uv} = (\theta_{11}, \theta_{12}, \ldots, \theta_{1,M_v},\; \theta_{21}, \ldots, \theta_{M_u, M_v})^\top
\in \mathbb{R}^{M_u \cdot M_v}
$$

For observation $r$, the bivariate basis vector is the **Kronecker product**:

$$
z_{r,uv} = z_u(x_{ru}) \otimes z_v(x_{rv}) \in \mathbb{R}^{M_u \cdot M_v}
$$

so that:

$$
f_{u,v}(x_{ru}, x_{rv}) = z_{r,uv}^\top\, \theta_{uv}
= \big(z_u(x_{ru}) \otimes z_v(x_{rv})\big)^\top \theta_{uv}
$$

The **interaction design matrix** is:

$$
Z_{uv} \in \mathbb{R}^{n \times M_u M_v},
\qquad \text{row } r: \quad z_{r,uv}^\top = z_u(x_{ru})^\top \otimes z_v(x_{rv})^\top
$$

Using the mixed-product property of Kronecker products:

$$
Z_{uv} = Z_u \odot Z_v
$$

where $\odot$ denotes the **row-wise Kronecker product** (Khatri–Rao product): row $r$ of $Z_{uv}$ is $z_{ru}^\top \otimes z_{rv}^\top$.

### 2.3 Identifiability for the Interaction

The interaction $f_{u,v}$ must be a **pure interaction** — it must not contain any component that looks like a main effect of $u$, a main effect of $v$, or a constant. This requires the ANOVA-type constraints on the coefficient grid $\{\theta_{ml}\}$:

**Constraint (i) — no $f_v$ component:**

$$
\sum_{m=1}^{M_u} \theta_{ml} = 0 \qquad \text{for each } l = 1, \ldots, M_v
$$

(Each column of the $\theta$ grid sums to zero.)

**Constraint (ii) — no $f_u$ component:**

$$
\sum_{l=1}^{M_v} \theta_{ml} = 0 \qquad \text{for each } m = 1, \ldots, M_u
$$

(Each row of the $\theta$ grid sums to zero.)

**Constraint (iii) — no intercept component:**

$$
\sum_{m=1}^{M_u} \sum_{l=1}^{M_v} \theta_{ml} = 0
$$

Note: constraint (iii) is **implied** by either (i) or (ii), so (i) + (ii) are the independent constraints.

### 2.4 Enforcing Identifiability via $T_u \otimes T_v$

**Claim:** Setting $\theta_{uv} = (T_u \otimes T_v)\, \beta_{uv}$ enforces all three constraints simultaneously.

**Proof of Constraint (i):**

Fix column index $l$. We need to show $\sum_{m=1}^{M_u} \theta_{ml} = 0$.

In the vectorized form, the entries $\theta_{1l}, \theta_{2l}, \ldots, \theta_{M_u, l}$ correspond to positions $(1-1)\cdot M_v + l,\; (2-1)\cdot M_v + l,\; \ldots,\; (M_u-1)\cdot M_v + l$ in $\theta_{uv}$.

Using the Kronecker product structure $\theta_{uv} = (T_u \otimes T_v)\, \beta_{uv}$, the $(m, l)$-th entry of the $\theta$-grid is:

$$
\theta_{ml} = \sum_{i=1}^{M_u-1} \sum_{j=1}^{M_v-1} [T_u]_{mi}\, [T_v]_{lj}\; [\beta_{uv}]_{ij}
= \left(\sum_{i} [T_u]_{mi}\right) \cdot (\ldots)
$$

Wait — let's be more precise. Summing over $m$:

$$
\sum_{m=1}^{M_u} \theta_{ml}
= \sum_{m=1}^{M_u} \sum_{i=1}^{M_u-1} \sum_{j=1}^{M_v-1} [T_u]_{mi}\, [T_v]_{lj}\, \beta_{ij}
= \sum_{i} \sum_{j} \underbrace{\left(\sum_{m} [T_u]_{mi}\right)}_{= [\mathbf{1}^\top T_u]_i = 0} [T_v]_{lj}\, \beta_{ij}
= 0
$$

Since $\mathbf{1}^\top T_u = \mathbf{0}^\top$ (proved in Section 1.5), the inner sum vanishes for every $i$. $\square$

**Proof of Constraint (ii):**

Fix row index $m$. Sum over $l$:

$$
\sum_{l=1}^{M_v} \theta_{ml}
= \sum_{l} \sum_{i} \sum_{j} [T_u]_{mi}\, [T_v]_{lj}\, \beta_{ij}
= \sum_{i} \sum_{j} [T_u]_{mi}\, \underbrace{\left(\sum_{l} [T_v]_{lj}\right)}_{= [\mathbf{1}^\top T_v]_j = 0} \beta_{ij}
= 0
$$

Since $\mathbf{1}^\top T_v = \mathbf{0}^\top$. $\square$

**Constraint (iii) follows** from either (i) or (ii). $\square$

### 2.5 The Identified Interaction Design Matrix

Define the **identified interaction coefficient vector**:

$$
\beta_{uv} \in \mathbb{R}^{(M_u - 1)(M_v - 1)}
$$

with the relation $\theta_{uv} = (T_u \otimes T_v)\, \beta_{uv}$.

The **identified interaction design matrix** is:

$$
W_{uv} = Z_{uv}\, (T_u \otimes T_v) \in \mathbb{R}^{n \times (M_u - 1)(M_v - 1)}
$$

so that:

$$
\mathbf{f}_{uv} = W_{uv}\, \beta_{uv}
$$

**Dimension count:** With $M_u = M_v = 20$ knots per covariate, each interaction block has $(20-1)(20-1) = 361$ identified parameters. With $p = 3$ covariates there are $\binom{3}{2} = 3$ interaction pairs, giving $3 \times 361 = 1{,}083$ interaction parameters total.

---

## 3. The 2D Penalty (Option C: Curvature of Surface Heights)

### 3.1 Motivation

Since $\theta_{ml} = f_{u,v}(\tau_{mu}, \tau_{lv})$ are surface heights at grid points, penalizing the **curvature** of the surface corresponds to penalizing **second differences** of $\theta$ along each axis of the grid.

This is the direct 2D generalization of the 1D RW2 penalty, which penalizes $\sum_m (\theta_{m+1} - 2\theta_m + \theta_{m-1})^2$ — the squared second differences of knot ordinates.

### 3.2 Second Differences on the Grid

**Along the $u$-axis** (varying $m$, fixed $l$):

$$
\Delta^2_m\, \theta_{ml} = \theta_{m+1,l} - 2\,\theta_{ml} + \theta_{m-1,l}
$$

**Along the $v$-axis** (varying $l$, fixed $m$):

$$
\Delta^2_l\, \theta_{ml} = \theta_{m,l+1} - 2\,\theta_{ml} + \theta_{m,l-1}
$$

### 3.3 The Total Penalty

$$
\mathcal{P}(\theta_{uv}) = \underbrace{\sum_{l=1}^{M_v} \sum_{m=2}^{M_u-1} \big(\Delta^2_m\, \theta_{ml}\big)^2}_{\text{curvature along } u}
\;+\;
\underbrace{\sum_{m=1}^{M_u} \sum_{l=2}^{M_v-1} \big(\Delta^2_l\, \theta_{ml}\big)^2}_{\text{curvature along } v}
$$

### 3.4 Matrix Form on $\theta_{uv}$

Let $D_2^{(u)} \in \mathbb{R}^{(M_u - 2) \times M_u}$ and $D_2^{(v)} \in \mathbb{R}^{(M_v - 2) \times M_v}$ be the 1D second-difference matrices.

Define the 1D penalty matrices:

$$
K_u^{(\theta)} = \big(D_2^{(u)}\big)^\top D_2^{(u)} \in \mathbb{R}^{M_u \times M_u}
$$

$$
K_v^{(\theta)} = \big(D_2^{(v)}\big)^\top D_2^{(v)} \in \mathbb{R}^{M_v \times M_v}
$$

Then the 2D penalty matrix on the vectorized $\theta_{uv}$ is:

$$
\boxed{
K_{uv}^{(\theta)} = K_u^{(\theta)} \otimes I_{M_v} + I_{M_u} \otimes K_v^{(\theta)}
}
$$

**Why this works:** $K_u^{(\theta)} \otimes I_{M_v}$ applies the $u$-direction second differences to each "column" $l$ independently. $I_{M_u} \otimes K_v^{(\theta)}$ applies the $v$-direction second differences to each "row" $m$ independently. The sum penalizes curvature in both directions.

### 3.5 Transforming to the Identified Coefficients

Using $\theta_{uv} = (T_u \otimes T_v)\, \beta_{uv}$:

$$
\mathcal{P} = \theta_{uv}^\top\, K_{uv}^{(\theta)}\, \theta_{uv}
= \beta_{uv}^\top\, \underbrace{(T_u \otimes T_v)^\top\, K_{uv}^{(\theta)}\, (T_u \otimes T_v)}_{K_{uv}^{(\beta)}}\, \beta_{uv}
$$

Expanding using the mixed-product property $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$:

$$
K_{uv}^{(\beta)}
= (T_u^\top K_u^{(\theta)} T_u) \otimes (T_v^\top T_v)
+ (T_u^\top T_u) \otimes (T_v^\top K_v^{(\theta)} T_v)
$$

If we define the 1D identified penalties and Gram matrices:

$$
K_u^{(\beta)} = T_u^\top K_u^{(\theta)} T_u, \qquad G_u = T_u^\top T_u
$$

$$
K_v^{(\beta)} = T_v^\top K_v^{(\theta)} T_v, \qquad G_v = T_v^\top T_v
$$

then:

$$
\boxed{
K_{uv}^{(\beta)} = K_u^{(\beta)} \otimes G_v + G_u \otimes K_v^{(\beta)}
}
$$

**Note:** $G_u = T_u^\top T_u$ and $G_v = T_v^\top T_v$ are not identity matrices (they arise from the sum-to-zero contrast). If $T$ has the specific structure in Section 1.5, then $G = T^\top T$ has $G_{11} = M - 1$ (first row/col), $G_{1k} = G_{k1} = -1$ for $k \geq 2$, and $G_{kk} = 2$, $G_{kl} = 1$ for $k \neq l$, $k,l \geq 2$. This is a known, fixed matrix.

**Special case check:** If $T_u$ and $T_v$ were orthogonal matrices (which they are not in general), then $G_u = G_v = I$ and we would recover $K_{uv}^{(\beta)} = K_u^{(\beta)} \otimes I + I \otimes K_v^{(\beta)}$, the simple Kronecker sum. With our specific $T$ matrices, the $G$ factors introduce cross-coupling that reflects the sum-to-zero constraint structure.

---

## 4. The Full Model with Interactions

### 4.1 Model Equation

$$
Y(s_r) = \mu + \sum_{j=1}^{p} f_j\big(X_j(s_r)\big)
+ \sum_{u < v} f_{u,v}\big(X_u(s_r), X_v(s_r)\big)
+ b(s_r) + \varepsilon_r
$$

for $r = 1, \ldots, n$.

### 4.2 Matrix Form

$$
\mathbf{y} = \mu\, \mathbf{1}_n + \sum_{j=1}^{p} W_j\, \beta_j
+ \sum_{u < v} W_{uv}\, \beta_{uv}
+ \mathbf{b} + \boldsymbol{\varepsilon}
$$

where:
- $W_j = Z_j\, T_j \in \mathbb{R}^{n \times (M_j - 1)}$ — main effect design matrices (existing)
- $W_{uv} = Z_{uv}\, (T_u \otimes T_v) \in \mathbb{R}^{n \times (M_u-1)(M_v-1)}$ — interaction design matrices (new)
- $\mathbf{b} \sim N(\mathbf{0}, \sigma^2 R(\rho, \nu))$ — spatial random effect (existing)
- $\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \tau^2 I_n)$ — nugget (existing)

### 4.3 Combined Design Matrix

Define the full parameter vector and combined design matrix:

$$
\boldsymbol{\eta} = (\mu, \beta_1^\top, \ldots, \beta_p^\top, \beta_{12}^\top, \beta_{13}^\top, \ldots, \beta_{p-1,p}^\top)^\top
$$

$$
H = \big[\mathbf{1}_n \;\big|\; W_1 \;\big|\; \cdots \;\big|\; W_p \;\big|\; W_{12} \;\big|\; W_{13} \;\big|\; \cdots \;\big|\; W_{p-1,p}\big]
$$

so the model is:

$$
\mathbf{y} = H\, \boldsymbol{\eta} + \mathbf{b} + \boldsymbol{\varepsilon}
$$

This has the same structure as your existing model — $H$ just has additional columns for the interaction blocks.


---

## 5. Summary: What Is Inherited vs. What Is New

| Component | Status |
|-----------|--------|
| LS basis functions $\Phi, \Psi$ | Existing |
| $A^{-1}C$ reduction to $z_j(x)$ | Existing |
| Contrast $T_j$ with $\mathbf{1}^\top T_j = 0$ | Existing |
| 1D design matrix $W_j = Z_j T_j$ | Existing |
| RW2 penalty $K_j^{(\beta)}$ | Existing |
| Tensor product $z_u \otimes z_v$ giving 2D LS basis | **New** |
| 2D knot-ordinate property $\theta_{ml} = f(\tau_m, \tau_l)$ | **New** |
| ANOVA identifiability via $T_u \otimes T_v$ | **New** |
| 2D penalty as Kronecker sum of 1D penalties | **New** |
| Spatial random effect $b$, collapsed Gibbs sampler | Existing (unchanged) |


---

## Appendix: Notation Reference

| Symbol | Meaning | Dimension |
|--------|---------|-----------|
| $M_j$ | Number of knots for covariate $j$ | scalar |
| $d_j = M_j - 1$ | Number of identified coefficients per main effect | scalar |
| $\tau_{mj}$ | $m$-th knot for covariate $j$ | scalar |
| $z_j(x)$ | Reduced LS basis vector for covariate $j$ | $M_j \times 1$ |
| $Z_j$ | Design matrix for covariate $j$ (unreduced) | $n \times M_j$ |
| $T_j$ | Sum-to-zero contrast matrix | $M_j \times (M_j - 1)$ |
| $W_j = Z_j T_j$ | Identified design matrix (main effect) | $n \times (M_j - 1)$ |
| $\beta_j$ | Identified main effect coefficients | $(M_j - 1) \times 1$ |
| $K_j^{(\beta)}$ | RW2 penalty on identified coefficients | $(M_j-1) \times (M_j-1)$ |
| $Z_{uv}$ | Interaction design (Khatri–Rao of $Z_u, Z_v$) | $n \times M_u M_v$ |
| $T_u \otimes T_v$ | Interaction contrast | $M_u M_v \times (M_u-1)(M_v-1)$ |
| $W_{uv}$ | Identified interaction design matrix | $n \times (M_u-1)(M_v-1)$ |
| $\beta_{uv}$ | Identified interaction coefficients | $(M_u-1)(M_v-1) \times 1$ |
| $K_{uv}^{(\beta)}$ | 2D penalty on identified interaction coefficients | $(M_u-1)(M_v-1) \times (M_u-1)(M_v-1)$ |
