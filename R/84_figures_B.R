## =====================================================================
## 84_figures_B.R -- the two figures of PaperB, from the stored fits
## =====================================================================
## Produces, with seeded posterior sampling (the earlier figures were
## drawn from unseeded draws, so they could not be reproduced exactly):
##
##   figures/figure_field.pdf    the shared field in 2018 at the 373
##                               survey locations: (a) posterior mean on
##                               the logit scale, (b) posterior SD.
##   figures/figure_species.pdf  predicted occupancy in 2018 for the
##                               twelve species: (a) posterior mean,
##                               (b) posterior SD.
##
## MODEL selects the fit the species panels are drawn from:
##   "M1" the reported pooled-detection model, read from fits.rds;
##   "M3" the grouped-detection model of Section 5.8, which this script
##        refits (about 7 minutes) rather than storing, since the fitted
##        object is several GB. The field figure always comes from M1,
##        whose hyperparameters Section 5.2 reports.
## Output names carry the model: figure_species_M3.pdf, and so on.
##
## Points are plotted at the survey locations only. The field is not
## interpolated between them, because between the locations it is
## informed largely by the prior (PaperB, Section 5.6).
##
## The numbers quoted in Section 5.6 -- the range of the field's posterior
## mean in 2018, its median posterior SD, and the median ratio of |mean|
## to SD -- are recomputed and printed, so the text can be checked
## against the figure it describes.
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR  <- "results/real"
FIGDIR  <- "figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
SEED    <- 20260201L
N_DRAWS <- 1000
YEAR    <- 2018L
if (!exists("MODEL")) MODEL <- "M3"
stopifnot(MODEL %in% c("M1", "M3"))
message("species panels from ", MODEL)

fits   <- readRDS(file.path(OUTDIR, "fits.rds"))
m_msom <- fits$msom          # M1: also the source of the field figure

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP   <- length(SP_NAMES)
N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
N_ROW  <- N_SITE * N_YEAR
YEARS  <- 2010:2018
t_idx  <- match(YEAR, YEARS)
stopifnot(!is.na(t_idx))

coords   <- hbefTrends$coords / 1000
elev_s   <- as.numeric(scale(hbefTrends$occ.covs$elev))
mesh <- fmesher::fm_mesh_2d(
  boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
  max.edge = c(0.35, 1.2), min.angle = 21,
  offset = c(0.05, 1.5), cutoff = 0.35)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

## ---------------------------------------------------------------------
## Draws: field in YEAR, and occupancy in YEAR for each species
## ---------------------------------------------------------------------
s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = m_msom,
                                  seed = SEED, num.threads = "1:1")
nm <- rownames(s[[1]]$latent)
i1 <- function(pat) { i <- grep(pat, nm); stopifnot(length(i) == 1L); i }
i_field <- grep("^spatialfield:", nm); n_spde <- length(i_field) / N_YEAR
stopifnot(n_spde == mesh$n)
sl_year <- i_field[((t_idx - 1) * n_spde + 1):(t_idx * n_spde)]
i_int <- i1("^Int_occ:"); i_e1 <- i1("^elev_s:"); i_e2 <- i1("^elev_s2:")
i_sp  <- vapply(seq_len(N_SP), \(i) i1(sprintf("^species_id:%d$", i)), 1L)

fld <- matrix(NA_real_, N_SITE, N_DRAWS)
psi <- array(NA_real_, c(N_SITE, N_SP, N_DRAWS))
for (d in seq_len(N_DRAWS)) {
  v <- s[[d]]$latent
  w <- as.numeric(A_site %*% v[sl_year])
  fld[, d] <- w
  base <- v[i_int] + v[i_e1] * elev_s + v[i_e2] * elev_s^2 + w
  for (i in seq_len(N_SP)) psi[, i, d] <- plogis(base + v[i_sp[i]])
}
f_mean <- rowMeans(fld);        f_sd <- apply(fld, 1, sd)
p_mean <- apply(psi, c(1, 2), mean); p_sd <- apply(psi, c(1, 2), sd)

## ---------------------------------------------------------------------
## M3: refit with species-specific elevation slopes and grouped detection,
## and replace the species panels. The field figure stays with M1.
## ---------------------------------------------------------------------
if (MODEL == "M3") {
  ph  <- read.csv(file.path(OUTDIR, "detection_heterogeneity.csv"))
  ph  <- ph[match(SP_NAMES, ph$species), ]
  grp <- cut(ph$p_hat, c(-Inf, 0.30, 0.50, Inf), labels = FALSE)
  message("detection groups: ", paste(table(grp), collapse = "/"))

  stack_species <- function(sp) {
    y <- hbefTrends$y[sp, , , ]
    do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
  }
  stack_det_cov <- function(nmc) {
    z <- hbefTrends$det.covs[[nmc]]
    M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
    M <- (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE); M[is.na(M)] <- 0; M
  }
  K_VIS   <- dim(hbefTrends$y)[4]
  site_id <- rep(seq_len(N_SITE), times = N_YEAR)
  time_id <- rep(seq_len(N_YEAR), each = N_SITE)
  Xday <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")

  det_group_cols <- function(g) {
    lapply(sort(unique(g))[-1], function(l)
      matrix(rep(rep(as.numeric(g == l), each = N_ROW), K_VIS), ncol = K_VIS))
  }
  Xd <- build_det_X(c(list(do.call(rbind, rep(list(Xday), N_SP)),
                           do.call(rbind, rep(list(Xtod), N_SP))),
                      det_group_cols(grp)), K = K_VIS)
  Y_h <- do.call(rbind, lapply(SP_NAMES, stack_species))

  spde <- INLA::inla.spde2.pcmatern(mesh,
            prior.range = PRIORS$spde_real$range,
            prior.sigma = PRIORS$spde_real$sigma, constr = TRUE)
  Xocc <- data.frame(
    Int_occ = 1,
    elev_s  = rep(elev_s[site_id], times = N_SP),
    elev_s2 = rep(elev_s[site_id]^2, times = N_SP),
    species_id = rep(seq_len(N_SP), each = N_ROW),
    x = rep(coords[site_id, 1], times = N_SP),
    y = rep(coords[site_id, 2], times = N_SP),
    time = rep(time_id, times = N_SP))
  Xocc$sp_e1 <- Xocc$species_id; Xocc$sp_e2 <- Xocc$species_id
  stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                       time = Xocc$time, n_time = N_YEAR)
  f3 <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
    f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(sp_e1, elev_s,  model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(sp_e2, elev_s2, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))
  t0 <- Sys.time()
  m3 <- fit_occupancy(f3, c(as.list(Xocc), stf$idx, list(Y = Y_h, X = Xd)),
                      Xd, compute_cv = FALSE)
  message(sprintf("M3 fitted in %.1f min",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  s3 <- INLA::inla.posterior.sample(n = N_DRAWS, result = m3,
                                    seed = SEED, num.threads = "1:1")
  nm3 <- rownames(s3[[1]]$latent)
  j1  <- function(p) { i <- grep(p, nm3); stopifnot(length(i) == 1L); i }
  jf  <- grep("^spatialfield:", nm3); n3 <- length(jf) / N_YEAR
  sl3 <- jf[((t_idx - 1) * n3 + 1):(t_idx * n3)]
  j_int <- j1("^Int_occ:"); j_e1 <- j1("^elev_s:"); j_e2 <- j1("^elev_s2:")
  j_sp  <- vapply(seq_len(N_SP), \(i) j1(sprintf("^species_id:%d$", i)), 1L)
  j_s1  <- vapply(seq_len(N_SP), \(i) j1(sprintf("^sp_e1:%d$", i)), 1L)
  j_s2  <- vapply(seq_len(N_SP), \(i) j1(sprintf("^sp_e2:%d$", i)), 1L)

  psi3 <- array(NA_real_, c(N_SITE, N_SP, N_DRAWS))
  for (d in seq_len(N_DRAWS)) {
    v <- s3[[d]]$latent
    w <- as.numeric(A_site %*% v[sl3])
    for (i in seq_len(N_SP))
      psi3[, i, d] <- plogis(v[j_int] + v[j_sp[i]] + w +
                             (v[j_e1] + v[j_s1[i]]) * elev_s +
                             (v[j_e2] + v[j_s2[i]]) * elev_s^2)
  }
  p_mean_M1 <- p_mean
  p_mean <- apply(psi3, c(1, 2), mean); p_sd <- apply(psi3, c(1, 2), sd)
  cat(sprintf("\nmean occupancy in %d, M1 -> M3, by species:\n", YEAR))
  print(data.frame(species = SP_NAMES,
                   M1 = round(colMeans(p_mean_M1), 3),
                   M3 = round(colMeans(p_mean), 3),
                   ratio = round(colMeans(p_mean) / colMeans(p_mean_M1), 2)),
        row.names = FALSE)
}

cat(sprintf("\nfield in %d: posterior mean from %.2f to %.2f (span %.2f)\n",
            YEAR, min(f_mean), max(f_mean), diff(range(f_mean))))
cat(sprintf("median posterior SD %.2f; median |mean|/SD %.2f\n",
            median(f_sd), median(abs(f_mean) / f_sd)))
cat(sprintf("occupancy SD across species and sites: %.3f to %.3f\n",
            min(p_sd), max(p_sd)))

## ---------------------------------------------------------------------
## Plot helpers: points coloured on a shared scale, with a colour bar
## ---------------------------------------------------------------------
pal <- function(n) hcl.colors(n, "Viridis")
## panel sizes follow the shape of the survey area, so that asp = 1 does
## not leave a band of white space on either side
ar <- diff(range(coords[, 1])) / diff(range(coords[, 2]))
pad <- 0.02 * diff(range(coords[, 1]))
xlim <- range(coords[, 1]) + c(-pad, pad)
ylim <- range(coords[, 2]) + c(-pad, pad) * ar

panel <- function(z, zlim, main, cex = 0.5) {
  col <- pal(101)[1 + round(100 * (pmin(pmax(z, zlim[1]), zlim[2]) - zlim[1]) /
                              diff(zlim))]
  plot(coords, pch = 19, cex = cex, col = col, asp = 1,
       xlim = xlim, ylim = ylim,
       xlab = "", ylab = "", axes = FALSE, main = main, cex.main = 0.9)
  box(col = "grey70")
}
colbar <- function(zlim, lab) {
  op <- par(mar = c(2.5, 0.6, 1.2, 2.6)); on.exit(par(op))
  z <- seq(zlim[1], zlim[2], length.out = 101)
  image(1, z, matrix(z, 1), col = pal(101), axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.7); mtext(lab, 3, 0.3, cex = 0.65); box(col = "grey70")
}

## --- figure 1: the shared field -------------------------------------
pw <- 3.3                                   # panel width, inches
pdf(file.path(FIGDIR, "figure_field.pdf"),
    width = 2 * (pw + 0.62), height = pw / ar + 0.75)
layout(matrix(c(1, 2, 3, 4), 1), widths = c(1, 0.22, 1, 0.22))
par(mar = c(1.2, 1.2, 2.0, 0.4))
zl <- range(f_mean); panel(f_mean, zl, sprintf("(a) posterior mean, %d", YEAR), 0.8)
colbar(zl, "logit")
par(mar = c(1.2, 1.2, 2.0, 0.4))
zs <- range(f_sd); panel(f_sd, zs, "(b) posterior standard deviation", 0.8)
colbar(zs, "logit")
dev.off()

## --- figure 2: occupancy by species ---------------------------------
ord <- order(colMeans(p_mean))          # rarest first, as in the tables
lay <- rbind(cbind(matrix(1:12, 3, 4, byrow = TRUE), 13),
             cbind(matrix(14:25, 3, 4, byrow = TRUE), 26))
pw2 <- 1.85                                 # panel width, inches
ph2 <- pw2 / ar + 0.30                      # plus room for the species name
pdf(file.path(FIGDIR, sprintf("figure_species_%s.pdf", MODEL)),
    width = 4 * pw2 + 0.75, height = 6 * ph2 + 0.9)
layout(lay, widths = c(1, 1, 1, 1, 0.3))
zl <- c(0, 1)
par(mar = c(0.6, 0.6, 1.6, 0.4), oma = c(0, 0, 2.4, 0))
for (k in ord) panel(p_mean[, k], zl, SP_NAMES[k])
colbar(zl, expression(hat(psi)))
mtext("(a) posterior mean occupancy probability", 3, 0.4, outer = TRUE, cex = 0.95)
zs <- c(0, max(p_sd))
par(mar = c(0.6, 0.6, 1.6, 0.4))
for (j in seq_along(ord)) {
  panel(p_sd[, ord[j]], zs, SP_NAMES[ord[j]])
  ## the second block needs its own heading; oma covers only the first
  ## centred over the four columns, to match the heading of block (a):
  ## `at` is in the user coordinates of the first panel, so two panel widths
  ## to its right is the middle of the row
  if (j == 1L) mtext("(b) posterior standard deviation", side = 3, line = 2.4,
                     at = xlim[1] + 2.05 * diff(xlim), adj = 0.5,
                     cex = 0.95, xpd = NA)
}
colbar(zs, "SD")
dev.off()

write.csv(data.frame(year = YEAR, min = min(f_mean), max = max(f_mean),
                     span = diff(range(f_mean)), median_sd = median(f_sd),
                     median_ratio = median(abs(f_mean) / f_sd)),
          file.path(OUTDIR, "field_summary_2018.csv"), row.names = FALSE)

cat("\nwritten:\n  ", file.path(FIGDIR, "figure_field.pdf"),
    "\n  ", file.path(FIGDIR, sprintf("figure_species_%s.pdf", MODEL)), "\n", sep = "")
cat("Check the printed field summaries against Section 5.6 of the manuscript.\n")
cat("NOTE: inla.posterior.sample(seed=, num.threads='1:1') is NOT exactly\n")
cat("repeatable for this model: two calls with the same seed in one session\n")
cat("gave latent values differing by up to 0.0069. Summaries over 1000 draws\n")
cat("are stable to about 0.001 (means) and 0.01 (extremes), so the extremes\n")
cat("quoted in Section 5.6 are stated as approximate.\n")
