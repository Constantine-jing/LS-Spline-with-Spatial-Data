# MOC — Dissertation

Your home base. Everything starts here.

**Title:** Bayesian Spatial Additive Regression with Lancaster–Šalkauskas Splines and Matérn Gaussian Process Random Effects
**Advisor:** Sounak Chakraborty
**Department:** Statistics and Data Science, University of Missouri

---

## Current focus (April 2026)

Extending Chapter 1 with **pairwise interaction terms** $f_{u,v}(X_u, X_v)$ — the key novelty for completing Project 1.

- See [[Tensor-product LS basis]]
- Latest experiment: [[Interaction 2x2 v2]]
- Open question: dissertation notation alignment for the interaction construction

---

## The novelty stack

These are what makes this dissertation new. See [[Novelty stack]] for the full argument.

1. [[LS basis]] + [[RW2 prior]] — never previously combined
2. [[Matérn GP random effect]] on continuous locations within a fully Bayesian LS-spline framework
3. [[Collapsed Gibbs sampler]] marginalizing out the spatial random effect
4. [[Tensor-product LS basis]] pairwise interactions in a spatial setting

---

## Chapter 1 — Spatial additive LS-spline ✅ + extension in progress

**Core model**
$$
y_i = \beta_0 + \sum_{j=1}^{p} f_j(X_{ij}) + s(\mathbf{loc}_i) + \varepsilon_i
$$
plus pairwise interactions $f_{u,v}(X_{iu}, X_{iv})$ in the extension.

- Basis: [[LS basis]] with [[RW2 prior]]
- Spatial: [[Matérn GP random effect]]
- Estimation: [[Collapsed Gibbs sampler]]
- Extension: [[Tensor-product LS basis]] (current focus)

**Validation**
- [[Sim1 n1000]] · [[Sim2 n1000]] · [[Sim3 n1000]] · [[Sim4 plus n1000]]
- [[Interaction null test]] · [[Interaction 2x2 v2]]
- Robustness: [[Prior sensitivity B1-B5]]
- Comparison: [[REML vs Bayes]]

**Application**
- [[Munich rent data]]

**Outputs**
- Poster at Mizzou | Revvity Innovation Summit (March 2026)

## Chapter 2 — Variable selection (queued)

Spike-and-slab on $\tau^2_{s,j}$ for component selection.

## Chapter 3 — Scalability / ML extensions (paused April 2026)

Poisson extension and MCMC acceleration. SEER lung cancer data is queued — see [[SEER lung cancer queued]].

---

## Comparison papers

| Paper | Role |
|---|---|
| [[Lang and Brezger 2004]] | Primary blueprint for interactions; field benchmark |
| [[Chib and Greenberg 2010]] | LS basis in Bayesian regression (1D, no spatial) |
| [[Nandy Lim Maiti 2017]] | Frequentist additive + Matérn — comparison target |
| [[Eilers and Marx 2003]] | Kronecker-sum penalty source |
| [[Vehtari et al 2017]] | WAIC / LOO-CV (replaces DIC) |

---

## Key decisions and insights

The "why" behind the dissertation. These are the notes most worth re-reading before defenses or talks.

- [[Why LS over P-splines]]
- [[Why collapsed Gibbs]]
- [[Centering for ANOVA identifiability]]
- [[X1 undercoverage is basis resolution]]
- [[Prior sensitivity is a strength]]
- [[INLA as computational benchmark]]

---

## Quick links

- [[MOC - Methods]]
- [[MOC - Simulations]]
- [[MOC - Literature]]
- [[R file index]]
- [[Hellbender HPC notes]]
- [[Bug log]]
