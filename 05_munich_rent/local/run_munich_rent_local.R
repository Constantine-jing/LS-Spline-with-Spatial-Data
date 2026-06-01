# ============================================================
# run_munich_rent_local.R
# Munich rent — FULL DATA (n=3082) with tensor-product interaction
#
# Model (matches Lang & Brezger 2004):
#   rentsqm = mu + f1(area) + f2(yearc) + f12(area,yearc)
#             + v'gamma + b(s) + eps
#
# Runtime: ~4-5 days on i7-12700F  |  RAM: ~2 GB peak
# ============================================================

cat("Loading source files...\n")
source("ls_basis.R")
source("spatial_utils.R")
source("fit_spatial_reml.R")
source("marginal_utils.R")
source("baby_bayes.R")
source("gibbs_bayes.R")
source("gibbs_stage_c_full.R")
source("ls_interaction.R")
source("gibbs_interaction.R")

M_C <- 20; n_iter <- 5000; n_burn <- 1500; n_thin <- 1; nu <- 1.5
mh_sd_sigma2 <- 1.0; mh_sd_tau2 <- 0.3; mh_sd_rho <- 1.2
outdir <- "munich_rent_full"
dir.create(outdir, showWarnings = FALSE)
stopifnot(ls_tests())

# ---- Load data ----
load("rent99.rda"); load("rent99_polys.rda")
cat(sprintf("\n# Munich Rent FULL (n=%d) with interaction\n", nrow(rent99)))

y <- rent99$rentsqm
area_raw <- rent99$area; yearc_raw <- rent99$yearc
area_min <- min(area_raw); area_max <- max(area_raw)
yearc_min <- min(yearc_raw); yearc_max <- max(yearc_raw)
area_01  <- (area_raw - area_min)/(area_max - area_min)
yearc_01 <- (yearc_raw - yearc_min)/(yearc_max - yearc_min)
smooth_names <- c("area","yearc")

bath_bin <- as.integer(rent99$bath=="premium")
kitchen_bin <- as.integer(rent99$kitchen=="premium")
cheat_bin <- as.integer(rent99$cheating=="yes")
loc_good <- as.integer(rent99$location=="good")
loc_top <- as.integer(rent99$location=="top")
V_lin_noexp <- cbind(bath=bath_bin, kitchen=kitchen_bin, cheating=cheat_bin)
V_lin_exp <- cbind(bath=bath_bin, kitchen=kitchen_bin, cheating=cheat_bin,
                   loc_good=loc_good, loc_top=loc_top)

# ---- Spatial coords ----
centroid_df <- data.frame(
  district=as.integer(names(rent99.polys)),
  cx=sapply(rent99.polys, function(p) mean(p[,1])),
  cy=sapply(rent99.polys, function(p) mean(p[,2])))
idx <- match(rent99$district, centroid_df$district)
if (any(is.na(idx))) {
  keep <- !is.na(idx); rent99 <- rent99[keep,]; y <- y[keep]
  area_01 <- area_01[keep]; yearc_01 <- yearc_01[keep]
  area_raw <- area_raw[keep]; yearc_raw <- yearc_raw[keep]
  V_lin_noexp <- V_lin_noexp[keep,,drop=FALSE]
  V_lin_exp <- V_lin_exp[keep,,drop=FALSE]; idx <- idx[keep]
}
coords <- as.matrix(centroid_df[idx, c("cx","cy")]); rownames(coords) <- NULL
n <- nrow(coords)
coord_min <- apply(coords,2,min); coord_max <- apply(coords,2,max)
coord_range <- coord_max - coord_min
coords_01 <- sweep(coords,2,coord_min)/matrix(coord_range,n,2,byrow=TRUE)
set.seed(42)
coords_01 <- coords_01 + matrix(rnorm(n*2,sd=0.001),n,2)
coords_01 <- pmin(pmax(coords_01,0),1)
D <- pairdist(coords_01)
cat(sprintf("  n=%d  D: max=%.4f  median=%.4f\n", n, max(D), median(D[D>0])))

# ---- Build LS basis ----
cat("Building LS basis (M=20)...\n")
obj_area  <- ls_build_one_full(area_01, M=M_C)
obj_yearc <- ls_build_one_full(yearc_01, M=M_C)
W_area <- obj_area$W; W_yearc <- obj_yearc$W
K_main_list <- list(build_rw2_penalty(M_C-1), build_rw2_penalty(M_C-1))
int_12 <- ls_build_interaction(obj_area, obj_yearc, orthogonalize=TRUE)
W_int <- int_12$W_uv; K_int <- int_12$K_uv
cat(sprintf("  W_area:%dx%d  W_yearc:%dx%d  W_int:%dx%d\n",
            nrow(W_area),ncol(W_area),nrow(W_yearc),ncol(W_yearc),nrow(W_int),ncol(W_int)))

# ---- Extended prior precision with linear columns ----
build_interaction_prior_precision_ext <- function(
    col_map_full, K_main_list, K_int_list,
    tau2_s_main, tau2_s_int, kappa2=1e6, eps_ridge=1e-6,
    n_linear=0, kappa2_linear=4) {
  # p_total: intercept (1) + number of non-intercept positions (max_0idx + 1)
  all_idx <- c(unlist(col_map_full$main), unlist(col_map_full$interaction))
  p_total <- 1L + max(all_idx) + 1L
  Q0 <- matrix(0, p_total, p_total)
  Q0[1,1] <- 1/kappa2
  # Linear columns sit at positions 2..(1+n_linear)
  if (n_linear > 0) for (k in 2:(1+n_linear)) Q0[k,k] <- 1/kappa2_linear
  # Main effect blocks (indices already correct)
  for (j in seq_along(K_main_list)) {
    idx <- 1 + col_map_full$main[[j]]
    Kj <- K_main_list[[j]]
    Q0[idx,idx] <- Q0[idx,idx] + (1/tau2_s_main[j])*Kj + eps_ridge*diag(length(idx))
  }
  # Interaction blocks (indices already correct)
  for (k in seq_along(K_int_list)) {
    idx <- 1 + col_map_full$interaction[[k]]
    Kk <- K_int_list[[k]]
    Q0[idx,idx] <- Q0[idx,idx] + (1/tau2_s_int[k])*Kk + eps_ridge*diag(length(idx))
  }
  Q0
}

# ---- Run one model ----
run_model <- function(V_linear, linear_names, label,
                      init_rho=0.03, init_sigma2=0.30, init_tau2=3.85) {
  p_linear <- ncol(V_linear)
  cat(sprintf("\n### %s (n=%d, %d linear, %d total cols)\n",
              label, n, p_linear, 1+p_linear+ncol(W_area)+ncol(W_yearc)+ncol(W_int)))
  H <- cbind(1, V_linear, W_area, W_yearc, W_int)
  p_H <- ncol(H)

  # col_map: 0-indexed positions in eta AFTER the intercept.
  # eta layout: [intercept | p_linear linear | 19 area | 19 yearc | 361 int]
  # eta[1] = intercept (not in col_map)
  # eta[2..(1+p_linear)] = linear (not in col_map, handled by prior builder)
  # eta[(2+p_linear)..(1+p_linear+19)] = area main
  # etc.
  # 0-indexed: position k means eta[k+1].
  # So area starts at 0-indexed position p_linear.
  off <- p_linear
  cm_main <- list(off + seq_len(ncol(W_area)) - 1L,
                  off + ncol(W_area) + seq_len(ncol(W_yearc)) - 1L)
  off <- off + ncol(W_area) + ncol(W_yearc)
  cm_int <- list("1_2" = off + seq_len(ncol(W_int)) - 1L)

  # Verify: max 0-indexed position + 1 (for intercept) should equal p_H
  stopifnot(1 + max(unlist(cm_int)) == p_H - 1)  # last 0-indexed pos + 1 for intercept = p_H
  # Actually: p_H = 1 + max(0-indexed) + 1, so max(0-indexed) = p_H - 2
  cat(sprintf("  H: %dx%d  max_col_idx=%d (expect %d)\n",
              nrow(H), p_H, max(unlist(cm_int)), p_H - 2))

  orig_fn <- build_interaction_prior_precision
  build_interaction_prior_precision <<- function(col_map_full, K_main_list, K_int_list,
                                                  tau2_s_main, tau2_s_int, kappa2=1e6, eps_ridge=1e-6) {
    build_interaction_prior_precision_ext(col_map_full, K_main_list, K_int_list,
      tau2_s_main, tau2_s_int, kappa2=kappa2, eps_ridge=eps_ridge,
      n_linear=p_linear, kappa2_linear=4)
  }
  on.exit(build_interaction_prior_precision <<- orig_fn)

  R_init <- matern_cor(D, rho=init_rho, nu=nu)
  Sigma_init <- init_sigma2*R_init + init_tau2*diag(n)
  L_init <- chol(Sigma_init + diag(1e-6,n))
  y_w <- forwardsolve(t(L_init), y)
  H_w <- forwardsolve(t(L_init), H)
  eta_init <- as.vector(solve(crossprod(H_w)+diag(1e-4,ncol(H)), crossprod(H_w,y_w)))
  init <- list(eta=eta_init, sigma2=init_sigma2, tau2=init_tau2, rho=init_rho,
               tau2_s_main=rep(1.0,2), tau2_s_int=rep(0.1,1))

  t0 <- proc.time()
  gs <- gibbs_interaction_sampler(
    y=y, H=H, D=D, nu=nu,
    col_map_main=cm_main, K_main_list=K_main_list,
    col_map_int=cm_int, K_int_list=list(K_int),
    n_iter=n_iter, n_burn=n_burn, n_thin=n_thin, kappa2=1e6,
    a_sigma=2, b_sigma=1, a_tau=2, b_tau=0.3,
    a_smooth=1, b_smooth=0.005,
    log_rho_mu=log(median(D[D>0])), log_rho_sd=1.0,
    mh_sd_log_sigma2=mh_sd_sigma2, mh_sd_log_tau2=mh_sd_tau2,
    mh_sd_log_rho=mh_sd_rho, init=init, verbose=TRUE)
  t_el <- (proc.time()-t0)[3]

  eta_post <- colMeans(gs$eta_samples); b_post <- colMeans(gs$b_samples)
  cat(sprintf("  Done in %.1fh  sigma2=%.4f tau2=%.4f rho=%.4f\n",
              t_el/3600, mean(gs$sigma2_samples), mean(gs$tau2_samples), mean(gs$rho_samples)))
  cat(sprintf("  MH accept: sigma2=%.3f tau2=%.3f rho=%.3f\n",
              gs$accept_rate["sigma2"], gs$accept_rate["tau2"], gs$accept_rate["rho"]))
  for (k in 1:p_linear) {
    pm <- mean(gs$eta_samples[,1+k]); psd <- sd(gs$eta_samples[,1+k])
    cat(sprintf("    %s: %.3f (%.3f)\n", linear_names[k], pm, psd))
  }
  cat(sprintf("  Spatial b: [%.3f, %.3f] sd=%.3f\n", min(b_post), max(b_post), sd(b_post)))
  save(gs, eta_post, b_post, t_el, p_linear, linear_names, cm_main, cm_int,
       file=file.path(outdir, paste0(gsub(" ","_",label),".rda")))
  list(gs=gs, eta_post=eta_post, b_post=b_post, t_el=t_el,
       p_linear=p_linear, linear_names=linear_names,
       cm_main=cm_main, cm_int=cm_int)
}

# ---- Run ----
cat(sprintf("\nStarted: %s\n", Sys.time()))
resA <- run_model(V_lin_noexp, colnames(V_lin_noexp), "Model A experts excluded")
resB <- run_model(V_lin_exp, colnames(V_lin_exp), "Model B experts included")
cat(sprintf("Finished: %s\n", Sys.time()))

# ---- Plots ----
x_grid <- seq(0,1,length.out=101)
pl_B <- resB$p_linear; eta_samp <- resB$gs$eta_samples

# Main effects
pdf(file.path(outdir,"marginals_modelB.pdf"), width=10, height=5)
par(mfrow=c(1,2))
for (j in 1:2) {
  idx_j <- 1 + pl_B + resB$cm_main[[j]]
  if (j==1) { W_g <- obj_area$design_new(x_grid,"W",TRUE); x_o <- x_grid*(area_max-area_min)+area_min; xl <- "Floor space (m^2)"
  } else    { W_g <- obj_yearc$design_new(x_grid,"W",TRUE); x_o <- x_grid*(yearc_max-yearc_min)+yearc_min; xl <- "Year of construction" }
  cm <- t(apply(eta_samp[,idx_j,drop=FALSE], 1, function(b) as.vector(W_g %*% b)))
  fh <- colMeans(cm); fl <- apply(cm,2,quantile,0.025); fu <- apply(cm,2,quantile,0.975)
  plot(x_o, fh, type="n", ylim=range(c(fl,fu)), main=paste0("f(",smooth_names[j],")"), xlab=xl, ylab="")
  polygon(c(x_o,rev(x_o)), c(fl,rev(fu)), col=rgb(.2,.4,.8,.2), border=NA)
  lines(x_o, fh, lwd=2, col="blue"); abline(h=0, lty=3)
}
dev.off()

# Interaction surface (Figure 6c)
grid2d <- expand.grid(area=x_grid, yearc=x_grid)
Wgu <- obj_area$design_new(grid2d$area,"W",TRUE)
Wgv <- obj_yearc$design_new(grid2d$yearc,"W",TRUE)
Wgi <- khatri_rao_rowwise_R(Wgu, Wgv)
if (int_12$recipe$orthogonalize && !is.null(int_12$recipe$A_uv))
  Wgi <- Wgi - cbind(1, Wgu, Wgv) %*% int_12$recipe$A_uv
idx_int <- 1 + pl_B + resB$cm_int[[1]]
f12_grid <- as.vector(Wgi %*% resB$eta_post[idx_int])
f12_mat <- matrix(f12_grid, length(x_grid), length(x_grid))

pdf(file.path(outdir,"interaction_surface.pdf"), width=8, height=7)
persp(x_grid*(area_max-area_min)+area_min,
      x_grid*(yearc_max-yearc_min)+yearc_min,
      f12_mat, theta=-40, phi=25, expand=0.6,
      xlab="Floor space", ylab="Year of construction", zlab="f_12",
      main="Interaction: floor space x year of construction",
      col="lightblue", shade=0.3, border=NA, ticktype="detailed")
dev.off()

# Traces
pdf(file.path(outdir,"trace_modelB.pdf"), width=10, height=10)
par(mfrow=c(3,2))
plot(resB$gs$sigma2_samples,type="l",main=expression(sigma^2))
plot(density(resB$gs$sigma2_samples),main=expression(sigma^2))
plot(resB$gs$tau2_samples,type="l",main=expression(tau^2))
plot(density(resB$gs$tau2_samples),main=expression(tau^2))
plot(resB$gs$rho_samples,type="l",main=expression(rho))
plot(density(resB$gs$rho_samples),main=expression(rho))
dev.off()

pdf(file.path(outdir,"tau2s_trace.pdf"), width=12, height=4)
par(mfrow=c(1,3))
plot(resB$gs$tau2_s_main_samp[,1],type="l",main="tau2_s area",col="purple")
plot(resB$gs$tau2_s_main_samp[,2],type="l",main="tau2_s yearc",col="purple")
plot(resB$gs$tau2_s_int_samp[,1],type="l",main="tau2_s interaction",col="darkred")
dev.off()

# Spatial maps
plot_smap <- function(bp, title, sr) {
  bd <- tapply(bp, rent99$district, mean); nc <- 100
  pal <- colorRampPalette(c("blue","white","red"))(nc)
  bc <- function(b) { i <- round((b-sr[1])/diff(sr)*(nc-1))+1; pal[pmin(pmax(i,1),nc)] }
  plot(NULL, xlim=range(coords_01[,1]), ylim=range(coords_01[,2]), asp=1, main=title, xlab="x", ylab="y")
  for (d in names(rent99.polys)) {
    p <- rent99.polys[[d]]
    p01 <- sweep(p,2,coord_min)/matrix(coord_range,nrow(p),2,byrow=TRUE)
    v <- bd[d]; if (!is.na(v)) polygon(p01[,1],p01[,2],col=bc(v),border="gray30",lwd=0.3)
  }
  lv <- seq(sr[1],sr[2],length.out=5)
  legend("bottomleft",legend=sprintf("%.2f",lv),fill=bc(lv),title="b(s)",bty="n",cex=0.8)
}
bA <- tapply(resA$b_post,rent99$district,mean)
bB <- tapply(resB$b_post,rent99$district,mean)
sr <- range(c(bA,bB),na.rm=TRUE)
pdf(file.path(outdir,"spatial_experts_comparison.pdf"), width=14, height=7)
par(mfrow=c(1,2)); plot_smap(resA$b_post,"a) Experts excluded",sr); plot_smap(resB$b_post,"b) Experts included",sr)
dev.off()

# Summary
cat("\n--- Summary ---\n")
cat(sprintf("  ModelA: sigma2=%.4f tau2=%.4f rho=%.4f b_sd=%.4f\n",
            mean(resA$gs$sigma2_samples),mean(resA$gs$tau2_samples),mean(resA$gs$rho_samples),sd(resA$b_post)))
cat(sprintf("  ModelB: sigma2=%.4f tau2=%.4f rho=%.4f b_sd=%.4f\n",
            mean(resB$gs$sigma2_samples),mean(resB$gs$tau2_samples),mean(resB$gs$rho_samples),sd(resB$b_post)))
cat(sprintf("  Spatial sd ratio (incl/excl): %.3f\n", sd(resB$b_post)/sd(resA$b_post)))
cat(sprintf("  tau2_s_int (ModelB): %.5f\n", mean(resB$gs$tau2_s_int_samp[,1])))
comp <- data.frame(model=c("excl","incl"),
  sigma2=c(mean(resA$gs$sigma2_samples),mean(resB$gs$sigma2_samples)),
  tau2=c(mean(resA$gs$tau2_samples),mean(resB$gs$tau2_samples)),
  rho=c(mean(resA$gs$rho_samples),mean(resB$gs$rho_samples)),
  b_sd=c(sd(resA$b_post),sd(resB$b_post)),
  time_hr=c(resA$t_el/3600,resB$t_el/3600))
write.csv(comp, file.path(outdir,"comparison.csv"), row.names=FALSE)
cat(sprintf("Total time: %.1f hours\n", (resA$t_el+resB$t_el)/3600))
cat("Done.\n")
