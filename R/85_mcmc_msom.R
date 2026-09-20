## =====================================================================
## 85_mcmc_msom.R -- MCMC benchmark for the hierarchical model
## =====================================================================
## PaperB, Section 3.3. The INLA fits are checked against an independent
## MCMC implementation: stMsPGOcc() in spOccupancy, a multi-species
## multi-season spatial occupancy model fitted with Polya-Gamma data
## augmentation and a nearest-neighbour Gaussian process.
##
## WHAT IS AND IS NOT COMPARABLE. The two models are not the same model,
## and the comparison is only as good as this list:
##   * detection: stMsPGOcc estimates SPECIES-SPECIFIC coefficients, drawn
##     from community distributions. That is what the INLA occupancy
##     family cannot do (Sec. 3.2), so the MCMC fit is the reference for
##     detection, not a replication of M1.
##   * occupancy: species-specific intercepts and elevation terms, drawn
##     from community distributions, so it is closest to M2/M3.
##   * the field: a spatial factor model. With n.factors = 1 there is one
##     latent field and species have loadings on it, the first fixed to 1
##     for identifiability. The INLA model fixes ALL loadings to 1, so
##     n.factors = 1 is the closest analogue, not an identical one.
##   * time: ar1 = TRUE gives each species its own AR(1); INLA uses one
##     AR(1) shared by the community.
##   * covariance: exponential NNGP here, Matern (nu = 1) SPDE there.
##
## The comparison therefore answers: do the species-level quantities the
## paper reports -- average occupancy, and detection -- agree with an
## independent implementation that is free of the detection ceiling? It
## does not answer whether INLA and MCMC agree on identical models.
##
## TWO REFERENCE MODELS (ENGINE):
##   "stMsPGOcc" -- spatial factor version, closest to the INLA model but,
##      on these data, not convergent at the species level: the chains
##      disagree on the signs of the factor loadings (83% of species), and
##      R-hat reaches 3 or more. Fixing the first loading at 1 rules out
##      the usual sign symmetry, so this is genuine multimodality, and
##      reordering the species to anchor on the most common one does not
##      cure it.
##   "tMsPGOcc" -- the same multi-species multi-season model WITHOUT the
##      spatial field. It drops the component that fails to identify while
##      keeping what the comparison is for: species-specific detection and
##      species-level occupancy. Its psi is not directly comparable with a
##      spatial model's, and that must be said wherever it is used.
##
## STAGES. Run the pilot first to time it and to see the chains move:
##     PILOT <- TRUE; source("R/85_mcmc_msom.R")
## then the full run (hours; leave it overnight):
##     rm(PILOT); source("R/85_mcmc_msom.R")
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT. Argument names were taken from the
## stMsPGOcc documentation; check them against the installed version.
## =====================================================================

suppressPackageStartupMessages({
  library(spOccupancy)
  library(coda)
})
data(hbefTrends)

OUTDIR <- "results/mcmc"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
if (!exists("PILOT"))  PILOT  <- FALSE
if (!exists("ENGINE")) ENGINE <- "tMsPGOcc"
stopifnot(ENGINE %in% c("tMsPGOcc", "stMsPGOcc"))
message("engine: ", ENGINE)
SEED <- 20260201L
set.seed(SEED)

## --- MCMC settings ----------------------------------------------------
if (PILOT) {
  N_BATCH <- 50L; BATCH <- 25L; N_BURN <- 500L; N_THIN <- 1L; N_CHAINS <- 1L
} else {
  N_BATCH <- 1200L; BATCH <- 25L      # 30,000 iterations per chain
  N_BURN  <- 15000L; N_THIN <- 15L; N_CHAINS <- 3L   # 1000 kept per chain
}
N_FACTORS <- 1L          # closest analogue to the single shared field
N_THREADS <- 4L

## ANCHOR SPECIES. stMsPGOcc fixes the first species' loading on each
## latent factor at 1 for identifiability. In hbefTrends the first species
## is AMRE, the second rarest, whose data barely determine the sign of the
## field: chains then settle on opposite signs for every other species'
## loading (loading sign-switching), which inflates R-hat for the
## species-level coefficients without any of them being wrong. Ordering the
## species so that the anchor is the most frequently recorded one pins the
## sign. Set to FALSE to reproduce the alphabetical-order run.
ANCHOR_MOST_COMMON <- TRUE

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP     <- length(SP_NAMES)
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K_VIS    <- dim(hbefTrends$y)[4]

## --- data in spOccupancy form ----------------------------------------
## hbefTrends is already y[species, site, year, visit] with occ.covs,
## det.covs and coords; only the coordinate scale is changed, to km, so
## that the range is on the same scale as the INLA fits.
naive_all <- vapply(SP_NAMES, function(sp) {
  a <- hbefTrends$y[sp, , , ]
  rec <- apply(a, c(1, 2), function(v) any(v == 1, na.rm = TRUE))
  srv <- apply(!is.na(a), c(1, 2), any)
  mean(rec[srv])
}, numeric(1))
ord <- if (ANCHOR_MOST_COMMON) order(naive_all, decreasing = TRUE) else seq_len(N_SP)
message("anchor species (loading fixed at 1): ", SP_NAMES[ord[1]],
        " (naive ", round(naive_all[ord[1]], 3), ")")

dat <- hbefTrends
dat$y <- hbefTrends$y[ord, , , , drop = FALSE]
SP_NAMES <- SP_NAMES[ord]
naive_all <- naive_all[ord]
dat$coords <- hbefTrends$coords / 1000
dat$occ.covs$elev_s  <- as.numeric(scale(hbefTrends$occ.covs$elev))
dat$occ.covs$elev_s2 <- dat$occ.covs$elev_s^2
dat$det.covs$day_s <- (hbefTrends$det.covs$day - mean(hbefTrends$det.covs$day, na.rm = TRUE)) /
                       sd(hbefTrends$det.covs$day, na.rm = TRUE)
dat$det.covs$tod_s <- (hbefTrends$det.covs$tod - mean(hbefTrends$det.covs$tod, na.rm = TRUE)) /
                       sd(hbefTrends$det.covs$tod, na.rm = TRUE)

## Distances, to set the range prior on the same footing as the PC prior
d <- dist(dat$coords)
message(sprintf("distances (km): min %.2f, median %.2f, max %.2f",
                min(d), median(d), max(d)))

## --- fit --------------------------------------------------------------
common <- list(
  occ.formula = ~ elev_s + elev_s2,
  det.formula = ~ day_s + tod_s,
  data        = dat,
  ar1         = TRUE,
  n.batch     = N_BATCH, batch.length = BATCH, accept.rate = 0.43,
  n.burn      = N_BURN, n.thin = N_THIN, n.chains = N_CHAINS,
  n.omp.threads = N_THREADS, verbose = TRUE, n.report = 100)

t0 <- Sys.time()
out <- if (ENGINE == "stMsPGOcc") {
  do.call(stMsPGOcc, c(common, list(n.factors = N_FACTORS,
                                    cov.model = "exponential",
                                    NNGP = TRUE, n.neighbors = 15)))
} else {
  dat$coords <- NULL                      # not used by tMsPGOcc
  common$data <- dat
  do.call(tMsPGOcc, common)
}
run_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("\n%s finished in %.1f min", ENGINE, run_min))
tag_file <- paste0("mcmc_", ENGINE, if (PILOT) "_pilot" else "_full", ".rds")
saveRDS(out, file.path(OUTDIR, tag_file))

## loading signs by chain: the diagnostic for the problem described above
if (N_CHAINS > 1 && !is.null(out$lambda.samples)) {
  nper <- nrow(out$lambda.samples) / N_CHAINS
  bych <- split(seq_len(nrow(out$lambda.samples)), rep(seq_len(N_CHAINS), each = nper))
  lam  <- round(sapply(bych, function(i) colMeans(out$lambda.samples[i, , drop = FALSE])), 3)
  cat("\n=== factor loadings, posterior mean by chain ===\n"); print(lam)
  flip <- mean(apply(sign(lam), 1, function(v) length(unique(v)) > 1))
  cat(sprintf("species whose loading sign differs between chains: %.0f%%\n", 100 * flip))
}

## --- convergence ------------------------------------------------------
cat("\n=== convergence ===\n")
conv <- function(nm) {
  if (is.null(out[[nm]])) return(NULL)
  x <- as.mcmc(out[[nm]])
  rh <- if (N_CHAINS > 1 && !is.null(out$rhat[[sub("\\.samples$", "", nm)]]))
          max(out$rhat[[sub("\\.samples$", "", nm)]], na.rm = TRUE) else NA_real_
  es <- min(effectiveSize(x), na.rm = TRUE)
  data.frame(block = nm, max_rhat = rh, min_ess = round(es))
}
cv <- do.call(rbind, lapply(c("beta.comm.samples", "alpha.comm.samples",
                              "beta.samples", "alpha.samples",
                              "theta.samples"), conv))
print(cv, row.names = FALSE)
cat("Rule of thumb: R-hat below 1.05 and effective sample size above 400\n")
if (PILOT) cat("(the pilot will fail both; it is only a timing and sanity run)\n")

## --- species-level detection and occupancy ---------------------------
## detection probability at the covariate means: plogis(alpha_i intercept)
al <- out$alpha.samples          # columns are species x detection terms
int_cols <- grep("Intercept", colnames(al))
stopifnot(length(int_cols) == N_SP)
p_mcmc <- apply(plogis(al[, int_cols, drop = FALSE]), 2,
                function(v) c(mean = mean(v), lo = quantile(v, .025),
                              hi = quantile(v, .975)))

## average occupancy over SURVEYED site-years, per species
surveyed <- apply(!is.na(dat$y[1, , , , drop = TRUE]), c(1, 2), any)
psi <- out$psi.samples           # [sample, species, site, year]
stopifnot(length(dim(psi)) == 4L, dim(psi)[2] == N_SP)
psi_bar <- vapply(seq_len(N_SP), function(i) {
  m <- apply(psi[, i, , , drop = FALSE], c(3, 4), mean)   # site x year
  mean(m[surveyed])
}, numeric(1))

naive <- naive_all

tab <- data.frame(species = SP_NAMES, naive = naive,
                  p_mcmc = p_mcmc["mean", ], p_lo = p_mcmc[2, ], p_hi = p_mcmc[3, ],
                  psi_bar_mcmc = psi_bar, ratio_mcmc = psi_bar / naive)

## bring in the conditional estimates and the INLA fits, if present
f_ph <- "results/real/detection_heterogeneity.csv"
if (file.exists(f_ph)) {
  ph <- read.csv(f_ph)
  tab <- merge(tab, ph[, c("species", "p_hat", "se_p")], by = "species", all.x = TRUE)
}
f_lad <- "results/ladder/ladder_full_summary.csv"
if (file.exists(f_lad)) {
  lad <- read.csv(f_lad)
  for (m in c("M1", "M3")) {
    sub <- lad[lad$model == m, c("species", "ratio")]
    names(sub)[2] <- paste0("ratio_", m)
    tab <- merge(tab, sub, by = "species", all.x = TRUE)
  }
}
tab <- tab[order(tab$naive), ]
cat("\n=== species-level comparison ===\n")
print(tab, row.names = FALSE, digits = 3)

if ("p_hat" %in% names(tab)) {
  d <- tab$p_mcmc - tab$p_hat
  cat(sprintf("\ndetection, MCMC minus conditional estimate: between %+.3f and %+.3f, mean %+.3f\n",
              min(d), max(d), mean(d)))
  cat(sprintf("correlation of the two across species: %.3f\n", cor(tab$p_mcmc, tab$p_hat)))
}
if ("ratio_M3" %in% names(tab))
  cat(sprintf("psibar/naive: MCMC %.2f-%.2f, M3 %.2f-%.2f, M1 %.2f-%.2f\n",
              min(tab$ratio_mcmc), max(tab$ratio_mcmc),
              min(tab$ratio_M3), max(tab$ratio_M3),
              min(tab$ratio_M1), max(tab$ratio_M1)))

## --- community-level and field parameters ----------------------------
cat("\n=== community-level occupancy (beta.comm) ===\n")
print(summary(out$beta.comm.samples)$quantiles[, c("2.5%", "50%", "97.5%")], digits = 3)
cat("\n=== community-level detection (alpha.comm) ===\n")
print(summary(out$alpha.comm.samples)$quantiles[, c("2.5%", "50%", "97.5%")], digits = 3)
cat("\n=== temporal (and, for stMsPGOcc, spatial) parameters (theta) ===\n")
th <- summary(out$theta.samples)$quantiles[, c("2.5%", "50%", "97.5%")]
print(th, digits = 3)
if (ENGINE == "stMsPGOcc")
  cat("phi is the decay of the exponential covariance. Note that 3/phi is the\n",
      "distance at which the correlation falls to 0.05, whereas the INLA\n",
      "practical range of 3.3 km is where it falls to 0.13, so the two are\n",
      "not on the same definition.\n", sep = "")
if (any(grepl("^phi", rownames(th))))
  print(round(3 / th[grep("^phi", rownames(th)), c("97.5%", "50%", "2.5%")], 2))

## if the other engine has already been run, put the two side by side
f_other <- file.path(OUTDIR, paste0("mcmc_species_",
             setdiff(c("tMsPGOcc", "stMsPGOcc"), ENGINE), ".csv"))
if (file.exists(f_other)) {
  oth <- read.csv(f_other)[, c("species", "p_mcmc", "ratio_mcmc")]
  names(oth)[-1] <- paste0(names(oth)[-1], "_other")
  cmp <- merge(tab[, c("species", "p_mcmc", "ratio_mcmc")], oth, by = "species")
  cat("\n=== the two MCMC references side by side ===\n")
  print(cmp, row.names = FALSE, digits = 3)
  cat(sprintf("detection: max absolute difference %.3f, correlation %.3f\n",
              max(abs(cmp$p_mcmc - cmp$p_mcmc_other)),
              cor(cmp$p_mcmc, cmp$p_mcmc_other)))
}

write.csv(tab, file.path(OUTDIR, paste0("mcmc_species_", ENGINE, ".csv")),
          row.names = FALSE)
write.csv(cv,  file.path(OUTDIR, paste0("mcmc_convergence_", ENGINE, ".csv")),
          row.names = FALSE)
cat("\nwritten to ", OUTDIR, " (run time ", round(run_min, 1), " min)\n", sep = "")
