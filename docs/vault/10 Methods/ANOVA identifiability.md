# ANOVA identifiability

The constraint structure that makes the additive decomposition identifiable.

## The problem

The model
$$
y = \mu + \sum_j f_j(X_j) + \sum_{u<v} f_{u,v}(X_u, X_v) + b + \varepsilon
$$
is **not identifiable** as written:

- Adding a constant to $f_j$ and subtracting it from $\mu$ leaves $y$ unchanged
- Adding $g(X_u)$ to $f_{u,v}$ and subtracting it from $f_u$ leaves $y$ unchanged

Without constraints, the posterior is improper or pathologically multimodal.

## The standard fix (ANOVA-type decomposition)

Each component is constrained to integrate (or sum) to zero:

- $\int f_j(x) \, dx = 0$ for all $j$
- $\int f_{u,v}(x_u, x_v) \, dx_u = 0$ for all $x_v$
- $\int f_{u,v}(x_u, x_v) \, dx_v = 0$ for all $x_u$

In practice, with finite samples, we use **sum-to-zero over the design points** (see [[Wood 2006]]):

- $\sum_i f_j(X_{ij}) = 0$
- For interactions: rows and columns of $f_{u,v}$ at design points sum to zero

## How it's implemented

Via contrast matrices applied to the [[LS basis]] design matrices.

For an interaction $f_{u,v}$:
$$
\tilde{Z}_{uv} = (Z_u T_u) \odot (Z_v T_v)
$$
where $T_u$ projects out the column space corresponding to $f_u$ (and similarly for $T_v$). This removes the main-effect directions from the interaction space.

## Why centering matters in practice

A subtle but important consequence:

- If the **true** $f_j$ (data-generating function) is **not** mean-zero on $[0,1]$, you must center it before generating $y$
- Otherwise the simulation gives the constant to the intercept, and your fitted $f_j$ matches a centered version while the true function is not centered → spurious bias and inflated RMSE

This is exactly the bug that [[Centering for ANOVA identifiability]] documents — the $2\sin(\pi x)$ → $2(\sin(\pi x) - 2/\pi)$ fix in the simulation DGPs.

## Connects to

- [[Tensor-product LS basis]] — the contrast matrices live here
- [[Centering for ANOVA identifiability]] — the bug and the fix
- [[Wood 2006]] — sum-to-zero constraint reference
