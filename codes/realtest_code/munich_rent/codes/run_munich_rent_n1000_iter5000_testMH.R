# ============================================================
# run_munich_rent.R
# Real-data application: Munich rent index 1999
#
# Direct comparison with Lang & Brezger (2004, JCGS):
#   They used Bayesian P-splines + MRF spatial prior (BayesX).
#   We use LS-spline basis + Matern GP spatial prior (our model).
#
# Model:
#   rentsqm_i = mu + f1(area_i) + f2(yearc_i) + v_i'gamma + b_i + eps_i
#
# where:
#   rentsqm : net rent per square meter (Euro/m^2)
#   area    : floor space in m^2          -> nonparametric LS-spline
#   yearc   : year of construction        -> nonparametric LS-spline
#   v_i     : binary covariates           -> linear (bath, kitchen, cheating,
#                                            location dummies)
#   b_i     : spatial random effect        -> Matern GP on district centroids
#   eps_i   : iid noise
#
# Data: gamlss.data::rent99  (n=3082, from Infratest Sozialforschung 1998)
#       gamlss.data::rent99.polys (district boundary polygons -> centroids)
#
# Spatial handling:
#   L&B (2004) used a Markov random field on ~380 districts.
#   We instead compute district centroids from boundary polygons and
#   use a continuous Matern(nu=1.5) GP. Multiple observations in the
#   same district share identical coordinates -> duplicated rows in D.
#   This is fine: R(d=0) = 1 gives perfect within-district correlation
#   for the spatial effect, which mirrors the MRF "same region = same
#   spatial effect" logic.
#
# Implementation:
#   We fit the 2 continuous covariates (area, yearc) nonparametrically
#   with LS-spline bases (M=20 knots). Binary covariates enter as extra
#   linear columns in the design matrix H, placed after the intercept
#   and before the spline blocks. The sampler's col_map points ONLY to
#   the spline columns, so the binary columns receive a vague normal
#   prior (same as the intercept, controlled by kappa2).
#
#   IMPORTANT: this requires a small extension to build_block_prior_precision
#   to give non-intercept, non-col_map columns a vague prior instead of
#   zero precision (which would make Q0 singular). We handle this by
#   adding eps_ridge to all diagonal entries of Q0 at the end, or more
#   cleanly, by specifying n_linear extra columns that get 1/kappa2 prior.
#   We implement the latter below.
#
# Run: Rscript run_munich_rent.R
# ============================================================

source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
library(mvtnorm)

# ============================================================
# Settings
# ============================================================
M_C    <- 20         # knots for Bayes
M_reml <- 6          # knots for REML
n_iter <- 5000
n_burn <- 1500
n_thin <- 1
nu     <- 1.5        # Matern smoothness (fixed)

# MH proposal SDs
# Acceptance from first run: sigma2=0.608 tau2=0.191 rho=0.809
# Target: 0.20-0.40. Increase sigma2 and rho proposals.
mh_sd_sigma2 <- 1.0    # was 0.5, accept=0.608 -> increase
mh_sd_tau2   <- 0.3    # was 0.3, accept=0.191 -> keep (borderline low)
mh_sd_rho    <- 0.8    # was 0.3, accept=0.809 -> increase substantially

dir.create("munich_rent", showWarnings = FALSE)
stopifnot(ls_tests())

# ============================================================
# Load and prepare data
# ============================================================
# Pre-saved from gamlss.data (see save_munich_data.R)
load("rent99.rda")        # -> rent99
load("rent99_polys.rda")  # -> rent99.polys

# ---- Subsample for computational feasibility ----
# Full n=3082 requires ~3082^3 Cholesky per MH step -> days.
# Subsample to n_sub observations (stratified would be ideal,
# but random is fine for a first pass).
n_sub <- 1000
set.seed(123)
sub_idx <- sort(sample(nrow(rent99), n_sub))
rent99  <- rent99[sub_idx, ]

cat("\n============================================================\n")
cat("# Munich Rent Application\n")
cat("============================================================\n")
cat(sprintf("  n = %d observations (subsampled from 3082)\n", nrow(rent99)))
cat(sprintf("  n distinct districts = %d\n", length(unique(rent99$district))))

# ---- Response ----
y <- rent99$rentsqm
cat(sprintf("  rentsqm: mean=%.3f  sd=%.3f  range=[%.2f, %.2f]\n",
            mean(y), sd(y), min(y), max(y)))

# ---- Continuous covariates (nonparametric) ----
# Rescale to [0,1] for LS-spline basis
area_raw  <- rent99$area
yearc_raw <- rent99$yearc

area_min  <- min(area_raw);  area_max  <- max(area_raw)
yearc_min <- min(yearc_raw); yearc_max <- max(yearc_raw)

area_01  <- (area_raw  - area_min)  / (area_max  - area_min)
yearc_01 <- (yearc_raw - yearc_min) / (yearc_max - yearc_min)

X_smooth <- cbind(area = area_01, yearc = yearc_01)
p_smooth <- ncol(X_smooth)
smooth_names <- c("area", "yearc")

cat(sprintf("  area:  [%.0f, %.0f] m^2  -> scaled [0,1]\n", area_min, area_max))
cat(sprintf("  yearc: [%.0f, %.0f]      -> scaled [0,1]\n", yearc_min, yearc_max))

# ---- Binary covariates (linear) ----
# bath:     0=standard, 1=premium
# kitchen:  0=standard, 1=premium
# cheating: 0=no central heating, 1=yes
# location: factor with 3 levels (average, good, top) -> 2 dummies
bath_bin    <- as.integer(rent99$bath == "premium")
kitchen_bin <- as.integer(rent99$kitchen == "premium")
cheat_bin   <- as.integer(rent99$cheating == "yes")
loc_good    <- as.integer(rent99$location == "good")
loc_top     <- as.integer(rent99$location == "top")

V_linear <- cbind(bath = bath_bin, kitchen = kitchen_bin,
                  cheating = cheat_bin, loc_good = loc_good, loc_top = loc_top)
p_linear <- ncol(V_linear)
linear_names <- colnames(V_linear)

cat(sprintf("  %d binary covariates: %s\n", p_linear,
            paste(linear_names, collapse=", ")))

# ---- Spatial coordinates from district centroids ----
# rent99.polys is a list of polygons, one per district.
# Each element is a 2-column matrix (x, y).
# Compute centroid = column means of each polygon.
poly_names <- names(rent99.polys)
centroids  <- t(sapply(rent99.polys, function(p) colMeans(p)))
colnames(centroids) <- c("cx", "cy")
centroid_df <- data.frame(district = as.integer(poly_names),
                          cx = centroids[,1], cy = centroids[,2])

# Match each observation to its district centroid
idx <- match(rent99$district, centroid_df$district)
if (any(is.na(idx))) {
  # Some districts in data may not have polygons — drop those obs
  cat(sprintf("  WARNING: %d obs have no polygon match, dropping\n", sum(is.na(idx))))
  keep <- !is.na(idx)
  rent99   <- rent99[keep, ]
  y        <- y[keep]
  X_smooth <- X_smooth[keep, , drop = FALSE]
  V_linear <- V_linear[keep, , drop = FALSE]
  area_01  <- area_01[keep]
  yearc_01 <- yearc_01[keep]
  area_raw <- area_raw[keep]
  yearc_raw <- yearc_raw[keep]
  idx      <- idx[keep]
}

coords <- as.matrix(centroid_df[idx, c("cx", "cy")])
rownames(coords) <- NULL
n <- nrow(coords)
cat(sprintf("  Final n = %d (after polygon matching)\n", n))

# Rescale coordinates to [0,1] x [0,1] for numerical stability
coord_min <- apply(coords, 2, min)
coord_max <- apply(coords, 2, max)
coord_range <- coord_max - coord_min
coords_01 <- sweep(coords, 2, coord_min) / matrix(coord_range, n, 2, byrow = TRUE)

# Add small jitter so no two observations share exact coordinates.
# Without this, the n x n Matern correlation matrix has identical
# rows/columns for same-district obs (R=1), making sigma2*R + tau2*I
# nearly singular -> Cholesky failure.
# Jitter magnitude ~0.001 on [0,1] scale — keeps same-district obs
# very close (cor ~ 0.999+) while breaking exact singularity.
set.seed(42)
jitter_sd <- 0.001
coords_01 <- coords_01 + matrix(rnorm(n * 2, sd = jitter_sd), n, 2)
# Clamp back to [0,1]
coords_01 <- pmin(pmax(coords_01, 0), 1)
cat(sprintf("  Coords rescaled to [0,1]^2 + jitter (sd=%.4f)\n", jitter_sd))

D <- pairdist(coords_01)
cat(sprintf("  D: max=%.4f  median=%.4f\n", max(D), median(D[D > 0])))

# How many unique coordinate locations?
n_unique <- nrow(unique(round(coords_01, 6)))
cat(sprintf("  Unique locations: %d (of %d obs)\n", n_unique, n))


# ============================================================
# Extend build_block_prior_precision to handle extra linear cols
# ============================================================
# The original function assumes columns are:
#   col 1 = intercept
#   cols in col_map = spline basis
# We now have:
#   col 1           = intercept
#   cols 2 : (1+p_linear) = binary linear covariates
#   cols after that  = spline basis (referenced by col_map)
#
# Strategy: intercept gets vague 1/kappa2 prior,
# linear binary columns get a moderately informative prior (sd=10),
# spline blocks get RW2 penalty.
# sd=10 for binary coefs: plausible range of ±20 Euro/m^2 at 95% level.
# This is still weakly informative but prevents the sampler from
# wandering in an absurdly wide prior (sd=1000 was the old kappa2=1e6).

build_block_prior_precision_extended <- function(col_map, tau2_s, kappa2 = 1e6,
                                                  eps_ridge = 1e-6,
                                                  n_linear = 0,
                                                  kappa2_linear = 100) {
  # col_map here is ALREADY shifted (indices account for linear cols).
  # p_total = 1 (intercept) + n_linear + n_spline_cols
  p_covs   <- length(col_map)
  p_total  <- 1 + max(unlist(col_map))
  Q0       <- matrix(0, p_total, p_total)
  
  # Intercept: vague
  Q0[1, 1] <- 1 / kappa2
  
  # Linear columns: moderately informative (sd = sqrt(kappa2_linear))
  if (n_linear > 0) {
    for (k in 2:(1 + n_linear)) {
      Q0[k, k] <- 1 / kappa2_linear
    }
  }
  
  # Spline blocks: RW2 penalty
  # col_map already has the right indices (shifted), so just +1 for intercept
  for (j in 1:p_covs) {
    idx <- 1 + col_map[[j]]
    d_j <- length(idx)
    K_j <- build_rw2_penalty(d_j)
    Q0[idx, idx] <- Q0[idx, idx] + (1 / tau2_s[j]) * K_j + eps_ridge * diag(d_j)
  }
  Q0
}

# Monkey-patch: override the function used by gibbs_full_sampler
# We need to replace build_block_prior_precision inside the sampler.
# Instead, we create a wrapper sampler that handles the extended design.

# ============================================================
# Extended Gibbs sampler for mixed linear + spline design
# ============================================================
# Wraps gibbs_full_sampler by:
# 1) Adjusting col_map to account for n_linear offset
# 2) Overriding build_block_prior_precision

gibbs_full_sampler_extended <- function(y, H, D, nu = 1.5, col_map, n_linear = 0,
                                         ...) {
  # Shift col_map indices to account for linear columns
  col_map_shifted <- lapply(col_map, function(idx) idx + n_linear)
  
  # Temporarily override the global function
  orig_fn <- build_block_prior_precision
  
  build_block_prior_precision <<- function(col_map, tau2_s, kappa2 = 1e6,
                                            eps_ridge = 1e-6) {
    build_block_prior_precision_extended(col_map, tau2_s, kappa2 = kappa2,
                                          eps_ridge = eps_ridge,
                                          n_linear = n_linear)
  }
  
  on.exit(build_block_prior_precision <<- orig_fn)
  
  gibbs_full_sampler(y = y, H = H, D = D, nu = nu,
                     col_map = col_map_shifted, ...)
}


# ============================================================
# REML (M=6)
# ============================================================
cat("\n####################################################\n")
cat("# REML: M=6  (continuous covariates only, no linear terms)\n")
cat("####################################################\n")

# Use fit_ls_spatial which builds cbind(1, W) internally.
# Binary covariates excluded from REML — same as the simulation
# scripts where REML only handles the nonparametric + spatial part.
# The Bayes model adds them; this keeps the REML baseline clean.
t0_reml <- proc.time()
reml_obj <- fit_ls_spatial(y = y, X_raw = X_smooth, coords = coords_01,
                           M_vec = rep(M_reml, p_smooth), nu = nu,
                           rho_init = 0.1, lambda_init = 1.0,
                           jitter = 1e-6, verbose = TRUE)
t_reml <- (proc.time() - t0_reml)[3]
reml_obj$X_raw_for_marginal <- X_smooth

cat(sprintf("  REML done in %.1fs\n", t_reml))
cat(sprintf("  sigma2=%.4f  tau2=%.4f  rho=%.4f  lambda=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho,
            reml_obj$fit$tau2 / reml_obj$fit$sigma2))

# Extract REML coefficients
beta_reml <- as.numeric(reml_obj$fit$beta)
cat(sprintf("  Intercept (mu): %.4f\n", beta_reml[1]))

# Marginal curves for REML
x_grid <- seq(0, 1, length.out = 101)
reml_curves <- marginal_curves(reml_obj, x_grid = x_grid,
                                truth_f_list = NULL, clip = TRUE)
# NOTE: REML marginal plot deferred until after Bayes so we can share y-axes.


# ============================================================
# Full Bayes (M=20)
# ============================================================
cat(sprintf("\n####################################################\n"))
cat(sprintf("# Full Bayes: M=%d  n_iter=%d  n_burn=%d\n", M_C, n_iter, n_burn))
cat(sprintf("# MH sd: sigma2=%.1f  tau2=%.1f  rho=%.1f\n",
            mh_sd_sigma2, mh_sd_tau2, mh_sd_rho))
cat("####################################################\n")

des_bayes <- ls_additive_build(X_smooth, M_vec = rep(M_C, p_smooth))
H_bayes   <- cbind(1, V_linear, des_bayes$W)
cat(sprintf("  H: %d x %d  (1 intercept + %d linear + %d spline)\n",
            nrow(H_bayes), ncol(H_bayes), p_linear,
            ncol(des_bayes$W)))

# Initialise from REML
R_init     <- matern_cor(D, rho = reml_obj$fit$rho, nu = nu)
Sigma_init <- reml_obj$fit$sigma2 * R_init + reml_obj$fit$tau2 * diag(n)
L_init     <- chol(Sigma_init + diag(1e-6, n))
y_w        <- forwardsolve(t(L_init), y)
H_w        <- forwardsolve(t(L_init), H_bayes)
eta_init   <- as.vector(solve(crossprod(H_w) + diag(1e-4, ncol(H_bayes)),
                               crossprod(H_w, y_w)))
init <- list(eta    = eta_init,
             sigma2 = max(reml_obj$fit$sigma2, 0.01),
             tau2   = max(reml_obj$fit$tau2,   0.01),
             rho    = reml_obj$fit$rho,
             tau2_s = rep(1.0, p_smooth))
cat(sprintf("  Init from REML: rho=%.4f  sigma2=%.4f  tau2=%.4f\n",
            init$rho, init$sigma2, init$tau2))

# Prior on rho: diffuse log-normal.
# Median distance between distinct points gives a scale reference.
med_dist <- median(D[D > 0])
cat(sprintf("  Median pairwise dist = %.4f -> log_rho_mu = %.2f\n",
            med_dist, log(med_dist)))

t0_bayes <- proc.time()
gs <- gibbs_full_sampler_extended(
  y = y, H = H_bayes, D = D, nu = nu,
  col_map  = des_bayes$col_map,
  n_linear = p_linear,
  n_iter = n_iter, n_burn = n_burn, n_thin = n_thin,
  kappa2 = 1e6,
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 0.3,
  a_smooth = 1, b_smooth = 0.005,
  log_rho_mu = log(med_dist),
  log_rho_sd = 1.0,
  mh_sd_log_sigma2 = mh_sd_sigma2,
  mh_sd_log_tau2   = mh_sd_tau2,
  mh_sd_log_rho    = mh_sd_rho,
  init = init, verbose = TRUE
)
t_bayes <- (proc.time() - t0_bayes)[3]
cat(sprintf("  Full Bayes done in %.1fs\n", t_bayes))

# ---- Summaries ----
gsm <- gibbs_summary(gs)
rho_summary <- c(mean = mean(gs$rho_samples), sd = sd(gs$rho_samples),
                 quantile(gs$rho_samples, c(0.025, 0.5, 0.975)))
tau2s_means <- colMeans(gs$tau2_s_samples)
names(tau2s_means) <- smooth_names

cat(sprintf("  sigma2: mean=%.4f  95%%CI=[%.4f,%.4f]\n",
            gsm$sigma2["mean"], gsm$sigma2["2.5%"], gsm$sigma2["97.5%"]))
cat(sprintf("  tau2:   mean=%.4f  95%%CI=[%.4f,%.4f]\n",
            gsm$tau2["mean"],   gsm$tau2["2.5%"],   gsm$tau2["97.5%"]))
cat(sprintf("  rho:    mean=%.4f  95%%CI=[%.4f,%.4f]\n",
            rho_summary["mean"], rho_summary["2.5%"], rho_summary["97.5%"]))
cat(sprintf("  MH accept: sigma2=%.3f  tau2=%.3f  rho=%.3f\n",
            gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))
cat("  tau2_s (smoothing variances):\n"); print(round(tau2s_means, 5))

# Linear coefficients
eta_post_mean <- colMeans(gs$eta_samples)
cat("\n  Linear coefs (posterior mean +/- sd):\n")
for (k in 1:p_linear) {
  pm <- mean(gs$eta_samples[, 1 + k])
  psd <- sd(gs$eta_samples[, 1 + k])
  ci <- quantile(gs$eta_samples[, 1 + k], c(0.025, 0.975))
  cat(sprintf("    %10s: %.4f (%.4f)  95%%CI=[%.4f,%.4f]\n",
              linear_names[k], pm, psd, ci[1], ci[2]))
}

# ---- Spatial effect b ----
b_post <- colMeans(gs$b_samples)
cat(sprintf("\n  Spatial b: range=[%.3f, %.3f]  sd=%.3f\n",
            min(b_post), max(b_post), sd(b_post)))


# ============================================================
# Bayes marginal curves
# ============================================================
# Build synthetic obj for marginal plotting (spline part only)
bayes_spline_beta <- c(eta_post_mean[1],
                       eta_post_mean[(1 + p_linear + 1):length(eta_post_mean)])
obj_bayes <- list(
  des = des_bayes,
  X_fix = H_bayes,
  fit = list(beta = bayes_spline_beta),
  X_raw_for_marginal = X_smooth
)

# For marginal band plots, we need to adjust gs$eta_samples to
# remove the linear columns so the indexing matches col_map.
# Create a modified gs for plotting purposes:
gs_for_plot <- gs
gs_for_plot$eta_samples <- cbind(
  gs$eta_samples[, 1, drop = FALSE],               # intercept
  gs$eta_samples[, (1 + p_linear + 1):ncol(gs$eta_samples), drop = FALSE]  # spline
)

# Pre-compute Bayes bands for shared ylim
bayes_bands <- vector("list", p_smooth)
for (j in 1:p_smooth) {
  bayes_bands[[j]] <- gibbs_marginal_band(gs_for_plot, obj_bayes, j, x_grid,
                                           X_raw = X_smooth, truth_f = NULL, level = 0.95)
}

# Compute shared y-axis PER COVARIATE (across REML curve + Bayes band)
shared_ylim <- vector("list", p_smooth)
for (j in 1:p_smooth) {
  all_y <- c(reml_curves$fhat_grid[[j]],
             bayes_bands[[j]]$f_hat,
             bayes_bands[[j]]$lower,
             bayes_bands[[j]]$upper)
  rng <- range(all_y, finite = TRUE)
  pad <- 0.05 * diff(rng)
  shared_ylim[[j]] <- rng + c(-pad, pad)
}

# ---- Combined plot: REML (top row) vs Bayes (bottom row), shared y-axes ----
pdf("munich_rent/munich_marginals_comparison.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

# Row 1: REML
for (j in 1:p_smooth) {
  if (j == 1) {
    x_orig <- x_grid * (area_max - area_min) + area_min
    xlab <- "Floor space (m^2)"
  } else {
    x_orig <- x_grid * (yearc_max - yearc_min) + yearc_min
    xlab <- "Year of construction"
  }
  plot(x_orig, reml_curves$fhat_grid[[j]], type = "l", lwd = 2,
       main = paste0("REML (M=6): f(", smooth_names[j], ")"),
       xlab = xlab, ylab = paste0("f_", j, "(x)"),
       ylim = shared_ylim[[j]], col = "blue")
  abline(h = 0, lty = 3, col = "gray50")
}

# Row 2: Bayes
for (j in 1:p_smooth) {
  band <- bayes_bands[[j]]
  if (j == 1) {
    x_orig <- x_grid * (area_max - area_min) + area_min
    xlab <- "Floor space (m^2)"
  } else {
    x_orig <- x_grid * (yearc_max - yearc_min) + yearc_min
    xlab <- "Year of construction"
  }
  plot(x_orig, band$f_hat, type = "n",
       main = paste0("Bayes (M=20): f(", smooth_names[j], ")"),
       xlab = xlab, ylab = paste0("f_", j, "(x)"),
       ylim = shared_ylim[[j]])
  polygon(c(x_orig, rev(x_orig)), c(band$lower, rev(band$upper)),
          col = rgb(0.2, 0.4, 0.8, 0.2), border = NA)
  lines(x_orig, band$f_hat, lwd = 2, col = "blue")
  abline(h = 0, lty = 3, col = "gray50")
  legend("topleft", legend = c("posterior mean", "95% band"),
         col = c("blue", rgb(0.2, 0.4, 0.8, 0.4)),
         lty = c(1, NA), lwd = c(2, 10), bty = "n", cex = 0.85)
}
dev.off()
cat("  Plot: munich_rent/munich_marginals_comparison.pdf\n")


# ---- Trace plots ----
pdf("munich_rent/munich_bayes_trace.pdf", width = 10, height = 10)
plot_gibbs_trace_full(gs, true_sigma2 = NULL, true_tau2 = NULL, true_rho = NULL)
dev.off()
cat("  Plot: munich_rent/munich_bayes_trace.pdf\n")


# ---- Smoothing variance trace ----
pdf("munich_rent/munich_bayes_tau2s_trace.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
for (j in 1:p_smooth) {
  plot(gs$tau2_s_samples[, j], type = "l",
       main = bquote(paste("Trace: ", tau[s]^2, " (", .(smooth_names[j]), ")")),
       ylab = expression(tau[s]^2), xlab = "iteration",
       col = "purple")
}
dev.off()
cat("  Plot: munich_rent/munich_bayes_tau2s_trace.pdf\n")


# ============================================================
# Spatial effect map (by district)
# ============================================================
# Average b posterior mean by district
b_by_district <- tapply(b_post, rent99$district, mean)

pdf("munich_rent/munich_spatial_effect.pdf", width = 8, height = 8)
# Color scale
b_vals <- as.numeric(b_by_district)
b_range <- range(b_vals, na.rm = TRUE)
n_colors <- 100
pal <- colorRampPalette(c("blue", "white", "red"))(n_colors)
b_to_col <- function(b) {
  idx <- round((b - b_range[1]) / diff(b_range) * (n_colors - 1)) + 1
  idx <- pmin(pmax(idx, 1), n_colors)
  pal[idx]
}

plot(NULL, xlim = range(coords_01[,1]), ylim = range(coords_01[,2]),
     asp = 1, main = "Spatial effect b (posterior mean)",
     xlab = "x", ylab = "y")
for (dname in names(rent99.polys)) {
  d_int <- as.integer(dname)
  poly_xy <- rent99.polys[[dname]]
  # Rescale polygon coordinates to [0,1]^2
  poly_01 <- sweep(poly_xy, 2, coord_min) / matrix(coord_range, nrow(poly_xy), 2, byrow = TRUE)
  b_val <- b_by_district[dname]
  if (!is.na(b_val)) {
    polygon(poly_01[,1], poly_01[,2], col = b_to_col(b_val), border = "gray30", lwd = 0.3)
  }
}
# Legend
legend_vals <- seq(b_range[1], b_range[2], length.out = 5)
legend("bottomleft",
       legend = sprintf("%.2f", legend_vals),
       fill = b_to_col(legend_vals),
       title = "b(s)", bty = "n", cex = 0.8)
dev.off()
cat("  Plot: munich_rent/munich_spatial_effect.pdf\n")


# ============================================================
# Summary comparison
# ============================================================
cat("\n\n####################################################\n")
cat("# REML vs Bayes: Munich Rent 1999\n")
cat("####################################################\n\n")

cat("--- Variance estimates ---\n")
cat(sprintf("  REML:  sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            reml_obj$fit$sigma2, reml_obj$fit$tau2, reml_obj$fit$rho))
cat(sprintf("  Bayes: sigma2=%.4f  tau2=%.4f  rho=%.4f\n",
            gsm$sigma2["mean"], gsm$tau2["mean"], rho_summary["mean"]))

cat("\n--- Linear coefs (Bayes only) ---\n")
for (k in 1:p_linear) {
  pm <- mean(gs$eta_samples[, 1 + k])
  cat(sprintf("  %-10s: %.4f\n", linear_names[k], pm))
}

cat("\n--- Smoothing variances (Bayes tau2_s) ---\n")
for (j in 1:p_smooth) {
  cat(sprintf("  %s: %.5f\n", smooth_names[j], tau2s_means[j]))
}

cat(sprintf("\n--- Spatial b: range=[%.3f, %.3f]  sd=%.3f ---\n",
            min(b_post), max(b_post), sd(b_post)))

cat(sprintf("\n--- Timing: REML=%.1fs  Bayes=%.1fs ---\n", t_reml, t_bayes))


# ============================================================
# Save outputs
# ============================================================
comp_table <- data.frame(
  method   = c("REML", "Bayes"),
  sigma2   = c(reml_obj$fit$sigma2, gsm$sigma2["mean"]),
  tau2     = c(reml_obj$fit$tau2,   gsm$tau2["mean"]),
  rho      = c(reml_obj$fit$rho,    rho_summary["mean"]),
  time_sec = c(t_reml, t_bayes)
)
write.csv(comp_table, "munich_rent/munich_reml_vs_bayes.csv", row.names = FALSE)

# Bayes posterior summaries
bayes_summary <- data.frame(
  parameter = c("sigma2", "tau2", "rho", "intercept",
                linear_names, paste0("tau2_s_", smooth_names)),
  post_mean = c(gsm$sigma2["mean"], gsm$tau2["mean"], rho_summary["mean"],
                mean(gs$eta_samples[,1]),
                sapply(2:(1+p_linear), function(k) mean(gs$eta_samples[,k])),
                tau2s_means),
  post_sd   = c(gsm$sigma2["sd"], gsm$tau2["sd"], rho_summary["sd"],
                sd(gs$eta_samples[,1]),
                sapply(2:(1+p_linear), function(k) sd(gs$eta_samples[,k])),
                apply(gs$tau2_s_samples, 2, sd)),
  ci_025    = c(gsm$sigma2["2.5%"], gsm$tau2["2.5%"], rho_summary["2.5%"],
                quantile(gs$eta_samples[,1], 0.025),
                sapply(2:(1+p_linear), function(k) quantile(gs$eta_samples[,k], 0.025)),
                apply(gs$tau2_s_samples, 2, function(x) quantile(x, 0.025))),
  ci_975    = c(gsm$sigma2["97.5%"], gsm$tau2["97.5%"], rho_summary["97.5%"],
                quantile(gs$eta_samples[,1], 0.975),
                sapply(2:(1+p_linear), function(k) quantile(gs$eta_samples[,k], 0.975)),
                apply(gs$tau2_s_samples, 2, function(x) quantile(x, 0.975)))
)
write.csv(bayes_summary, "munich_rent/munich_bayes_summary.csv", row.names = FALSE)

cat("\nOutputs in munich_rent/:\n")
cat("  munich_marginals_comparison.pdf\n")
cat("  munich_bayes_trace.pdf\n")
cat("  munich_bayes_tau2s_trace.pdf\n")
cat("  munich_spatial_effect.pdf\n")
cat("  munich_reml_vs_bayes.csv\n")
cat("  munich_bayes_summary.csv\n")
cat("\nDone.\n")
