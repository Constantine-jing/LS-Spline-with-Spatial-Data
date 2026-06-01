# MOC — Methods

All methodological building blocks of the dissertation. Each is a standalone atomic note.

## Basis & smoothing

- [[LS basis]] — the core spline basis
- [[RW2 prior]] — penalty / prior on spline coefficients
- [[Tensor-product LS basis]] — interaction extension via Khatri-Rao
- [[Kronecker-sum RW2 penalty]] — 2D penalty for tensor product
- [[ANOVA identifiability]] — centering and contrast constraints

## Spatial

- [[Matérn GP random effect]] — spatial dependence on continuous locations
- [[Matérn covariance function]] — kernel parameterization

## Estimation

- [[Collapsed Gibbs sampler]] — marginalize out spatial random effect $b$
- [[MCMC convergence diagnostics]] — what we monitor

## Model selection & evaluation

- [[WAIC and LOO-CV]] — replaces DIC, see [[Vehtari et al 2017]]
- [[RMSE and coverage as metrics]]

## Comparison frameworks

- [[REML estimation]] — frequentist comparison via [[fit_spatial_reml.R]]
- [[INLA as computational benchmark]] — see also [[Why LS over P-splines]]
