## =====================================================================
## 63_avg_occupancy.R -- average occupancy and P(recorded), corrected
## =====================================================================
## Replaces the earlier version. Three corrections and one speed-up:
##
##  1. REPRODUCIBILITY. inla.posterior.sample() is now called with
##     seed = SEED and num.threads = "1:1". Without both, repeated calls
##     do not reproduce (set.seed() alone does not control it), so the
##     published psi_bar values could not be regenerated exactly.
##
##  2. VISITS ACTUALLY MADE. P(recorded) multiplied over all K visits,
##     including visits that did not take place, which overstates it for
##     the 18.6% of site-years with a missing visit. It now multiplies
##     over the visits made, using the NA pattern of the species' data.
##
##  3. P(recorded) IS REPORTED with its ratio to the observed frequency.
##     This is the quantity comparable with naive occupancy; psi_bar is
##     not (naive occupancy is a lower bound on occupancy).
##
##  4. SPEED. The hierarchical model is sampled ONCE and the draws are
##     reused for all twelve species, instead of being resampled per
##     species. With a fixed seed this also makes the species share a
##     single set of draws, which is what a joint posterior means.
##
## Because of 1 and 4 the numbers move slightly from the published table
## (third decimal of psi_bar). If the old file is present, the script
## prints the differences.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR  <- "results/real"
SEED    <- 20260201L
N_DRAWS <- 1000

fits   <- readRDS(file.path(OUTDIR, "fits.rds"))
m_msom <- fits$msom
ssom   <- fits$ssom

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP   <- length(SP_NAMES)
N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
K      <- dim(hbefTrends$y)[4]
N_ROW  <- N_SITE * N_YEAR

stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
}

coords   <- hbefTrends$coords / 1000
elev_s   <- as.numeric(scale(hbefTrends$occ.covs$elev))
elev_rep <- rep(elev_s, times = N_YEAR)
time_rep <- rep(seq_len(N_YEAR), each = N_SITE)

Y_all <- lapply(SP_NAMES, stack_species); names(Y_all) <- SP_NAMES
made  <- lapply(Y_all, \(Y) !is.na(Y))          # visits actually made
surv  <- lapply(made,  \(M) rowSums(M) > 0)

observed <- vapply(SP_NAMES, function(sp) {
  Y <- Y_all[[sp]]; ok <- surv[[sp]]
  mean(rowSums(Y[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0)
}, numeric(1))

## ---------------------------------------------------------------------
## Projector, as in the fit
## ---------------------------------------------------------------------
bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

## ---------------------------------------------------------------------
## Posterior mean occupancy, all species from one set of draws.
## `species` = NULL for a single-species fit.
## ---------------------------------------------------------------------
psi_mean <- function(fit, species = NULL) {
  s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = fit,
                                    seed = SEED, num.threads = "1:1")
  nm <- rownames(s[[1]]$latent)
  idx1 <- function(pat) {
    i <- grep(pat, nm)
    if (length(i) != 1L) stop("pattern '", pat, "' matched ", length(i))
    i
  }
  i_int <- idx1("^Int_occ:"); i_e1 <- idx1("^elev_s:"); i_e2 <- idx1("^elev_s2:")
  i_sp  <- if (is.null(species)) NULL else
             vapply(species, \(i) idx1(sprintf("^species_id:%d$", i)), 1L)
  i_field <- grep("^spatialfield:", nm)
  n_spde  <- length(i_field) / N_YEAR
  stopifnot(n_spde == mesh$n)

  ncol_out <- if (is.null(i_sp)) 1L else length(i_sp)
  acc <- matrix(0, N_ROW, ncol_out)
  for (d in seq_len(N_DRAWS)) {
    v   <- s[[d]]$latent
    fld <- numeric(N_ROW)
    for (tt in seq_len(N_YEAR)) {
      sl <- i_field[((tt - 1) * n_spde + 1):(tt * n_spde)]
      fld[time_rep == tt] <- as.numeric(A_site %*% v[sl])
    }
    base <- v[i_int] + v[i_e1] * elev_rep + v[i_e2] * elev_rep^2 + fld
    if (is.null(i_sp)) acc[, 1] <- acc[, 1] + plogis(base)
    else for (j in seq_along(i_sp)) acc[, j] <- acc[, j] + plogis(base + v[i_sp[j]])
  }
  acc / N_DRAWS
}

## Detection probability per row and visit, averaged over draws.
p_detect_mean <- function(fit, Xday, Xtod) {
  dr <- det_draws(fit, n = N_DRAWS)     # columns: intercept, day, tod
  pm <- matrix(NA_real_, N_ROW, K)
  for (j in seq_len(K)) {
    eta <- outer(rep(1, N_ROW), dr[, 1]) +
           outer(Xday[, j], dr[, 2]) +
           outer(Xtod[, j], dr[, 3])
    pm[, j] <- rowMeans(plogis(eta))
  }
  pm
}

## P(recorded at least once), over the visits ACTUALLY MADE
det_any <- function(pm, M) {
  pm[!M] <- 0
  1 - exp(rowSums(log1p(-pm)))
}

Xday <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")
Xday[is.na(Xday)] <- 0;       Xtod[is.na(Xtod)] <- 0

## ---------------------------------------------------------------------
## Compute
## ---------------------------------------------------------------------
message("hierarchical model: one set of ", N_DRAWS, " draws for all species")
psi_h_all <- psi_mean(m_msom, species = seq_len(N_SP))
p_h       <- p_detect_mean(m_msom, Xday, Xtod)

res <- list()
for (i in seq_len(N_SP)) {
  sp <- SP_NAMES[i]
  message("species ", sp)
  M  <- made[[sp]]; ok <- surv[[sp]]

  psi_h <- psi_h_all[, i]
  psi_s <- psi_mean(ssom[[sp]])[, 1]
  p_s   <- p_detect_mean(ssom[[sp]], Xday, Xtod)

  rec_h <- psi_h * det_any(p_h, M)
  rec_s <- psi_s * det_any(p_s, M)

  res[[i]] <- data.frame(
    species       = sp,
    observed      = observed[[sp]],
    psi_bar_msom  = mean(psi_h[ok]),
    psi_bar_ssom  = mean(psi_s[ok]),
    psi_mean_msom = plogis(m_msom$summary.fixed["Int_occ", "mean"] +
                           m_msom$summary.random$species_id$mean[i]),
    psi_mean_ssom = plogis(ssom[[sp]]$summary.fixed["Int_occ", "mean"]),
    prec_msom     = mean(rec_h[ok]),
    prec_ssom     = mean(rec_s[ok]))
  res[[i]]$ratio_msom <- res[[i]]$psi_bar_msom / res[[i]]$observed
  res[[i]]$ratio_ssom <- res[[i]]$psi_bar_ssom / res[[i]]$observed
  res[[i]]$prec_ratio_msom <- res[[i]]$prec_msom / res[[i]]$observed
  res[[i]]$prec_ratio_ssom <- res[[i]]$prec_ssom / res[[i]]$observed
}
out <- do.call(rbind, res)
out <- out[order(out$observed), ]

cat("\n", strrep("=", 92), "\n", sep = "")
cat("psi_mean_* : occupancy at the covariate and field means (not comparable with 'observed')\n")
cat("psi_bar_*  : average posterior occupancy over surveyed site-years\n")
cat("prec_*     : posterior predictive P(recorded at least once in the visits MADE)\n\n")
print(out[, c("species", "observed", "psi_bar_msom", "psi_bar_ssom",
              "ratio_msom", "ratio_ssom")], row.names = FALSE, digits = 3)

cat("\n--- P(recorded) against the observed frequency ---\n")
print(out[, c("species", "observed", "prec_msom", "prec_ssom",
              "prec_ratio_msom", "prec_ratio_ssom")], row.names = FALSE, digits = 3)
cat(sprintf("\nmean |P(recorded)/observed - 1|: hierarchical %.3f, independent %.3f\n",
            mean(abs(out$prec_ratio_msom - 1)), mean(abs(out$prec_ratio_ssom - 1))))
cat("Both classes are fitted to reproduce the observed frequency, so neither\n")
cat("is expected to depart from it; this comparison does not discriminate\n")
cat("between them (PaperB, Section 2.4). The discriminating check is the\n")
cat("number of detections per recorded site-year (68/70/52 scripts).\n")

cat("\n--- does either class fall below the observed frequency? ---\n")
for (nmn in c("psi_mean", "psi_bar", "prec"))
  cat(sprintf("  %-9s: hierarchical %d, independent %d (of %d species)\n", nmn,
              sum(out[[paste0(nmn, "_msom")]] < out$observed),
              sum(out[[paste0(nmn, "_ssom")]] < out$observed), nrow(out)))

## ---------------------------------------------------------------------
## What changed relative to the stored (unseeded, all-K) version?
## ---------------------------------------------------------------------
f_old <- file.path(OUTDIR, "avg_occupancy.csv")
if (file.exists(f_old)) {
  old <- read.csv(f_old)
  cmp <- merge(old[, c("species", "psi_bar_msom", "psi_bar_ssom",
                       "prec_msom", "prec_ssom")],
               out[, c("species", "psi_bar_msom", "psi_bar_ssom",
                       "prec_msom", "prec_ssom")],
               by = "species", suffixes = c("_old", "_new"))
  cat("\n--- change from the stored version (new - old) ---\n")
  d <- data.frame(species = cmp$species,
                  d_psi_bar_msom = cmp$psi_bar_msom_new - cmp$psi_bar_msom_old,
                  d_psi_bar_ssom = cmp$psi_bar_ssom_new - cmp$psi_bar_ssom_old,
                  d_prec_msom    = cmp$prec_msom_new    - cmp$prec_msom_old,
                  d_prec_ssom    = cmp$prec_ssom_new    - cmp$prec_ssom_old)
  print(d, row.names = FALSE, digits = 2)
  cat(sprintf("largest absolute change: %.4f\n",
              max(abs(as.matrix(d[, -1])), na.rm = TRUE)))
  file.rename(f_old, file.path(OUTDIR, "avg_occupancy_unseeded.csv"))
  cat("old file kept as avg_occupancy_unseeded.csv\n")
}

write.csv(out, file.path(OUTDIR, "avg_occupancy.csv"), row.names = FALSE)
cat("\nwritten: ", file.path(OUTDIR, "avg_occupancy.csv"), "\n", sep = "")
