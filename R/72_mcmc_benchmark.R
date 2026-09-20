## =====================================================================
## 72_mcmc_benchmark.R -- INLA-SPDE against MCMC (spOccupancy)
## =====================================================================
## The Limitations section states that no comparison against MCMC or
## against spOccupancy was made. That gap is closeable: M1 and M2 are
## ordinary multi-season single-species spatial occupancy models, and M2
## differs only in carrying one additional detection covariate, which
## happens to be another species' detection indicator. stPGOcc() fits
## exactly this class, with an arbitrary detection formula and an
## optional AR(1) term, on the same hbefTrends object.
##
## WHAT THIS ESTABLISHES, AND WHAT IT DOES NOT.
##
## It establishes whether the Laplace approximation and the MCMC sampler
## agree on this model and these data. It does not establish that the
## approximation is accurate in general: the occupancy log-likelihood is
## not log-concave, and one empirical comparison does not change that.
## Running the 800 simulation replicates under MCMC is not feasible and
## is not attempted.
##
## THE FIRST RUN REVEALED A STRUCTURAL DIFFERENCE worth stating plainly.
## stPGOcc with ar1 = TRUE fits a SEPARATE temporal random effect, with its
## own variance sigma.sq.t, alongside the spatial field. The INLA
## specification has one separable space-time field whose temporal
## structure is the AR(1) grouping. So spOccupancy carries one variance
## component more. In the first run it estimated sigma.sq.t = 1.60, the
## field's AR(1) rho collapsed to 0.16 against 0.92 under INLA, and the
## elevation coefficients moved by up to three posterior standard
## deviations, while every detection coefficient agreed to within half of
## one. The occupancy discrepancy is therefore a difference of model, not
## of inference method -- the same difference already noted for the
## archived code of Belmont et al.
##
## THE TWO MODELS ARE NOT IDENTICAL in other ways too. The
## INLA fit represents the field by an SPDE approximation to a Matern
## GRF; stPGOcc uses a Nearest Neighbor Gaussian Process, and for
## multi-season single-species models only NNGP is supported. Disagreement
## can therefore be a difference of model as much as of inference method.
## The priors also differ in parameterisation: PC priors on range and
## marginal standard deviation here, uniform priors on phi and sigma.sq
## there. Choose the spOccupancy priors to be as close as the two
## parameterisations allow, and report both.
##
## RUN TIME IS NOT A CLEAN COMPARISON either. Thread counts, neighbour
## counts, chain length and convergence criteria all enter. Report the
## settings alongside the timings and claim no advantage, consistent with
## the position already taken in the Strengths section.
## =====================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(spOccupancy)
  library(coda)
})

data(hbefTrends)
OUTDIR <- "results/mcmc_benchmark"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "mcmc_benchmark")

## ---------------------------------------------------------------------
## MCMC settings. n.batch * batch.length is the chain length per chain.
## These are starting values: check the diagnostics at the end and raise
## them if R-hat or ESS say to.
## ---------------------------------------------------------------------
## First run: 20,000 iterations per chain. The detection parameters were
## fully converged there (ESS 2752-3000, R-hat 1.00), but the occupancy
## intercept was not (ESS 107, R-hat 1.12). Raised threefold so that the
## occupancy comparison is defensible too; without it a reviewer can
## attribute the occupancy discrepancy to the sampler rather than to the
## structural difference between the two models, and we could not answer.
N_BATCH   <- 2400
BATCH_LEN <- 25            # 60,000 iterations per chain
N_BURN    <- 30000
N_THIN    <- 20            # keeps 1500 post-burn-in draws per chain
N_CHAINS  <- 3
N_NEIGHB  <- 15            # Datta et al. suggest 15 is usually enough
N_THREADS <- 1             # keep at 1: the timing comparison is against a
                           # single-threaded INLA fit (OPTS$num_threads)

## ---------------------------------------------------------------------
## Data, in spOccupancy's layout: y is sites x seasons x visits.
## ---------------------------------------------------------------------
SP_FOCAL <- "REVI"
SP_BIOTIC <- "OVEN"
sp_names <- dimnames(hbefTrends$y)[[1]]
stopifnot(SP_FOCAL %in% sp_names, SP_BIOTIC %in% sp_names)

y_revi <- hbefTrends$y[SP_FOCAL, , , ]
y_oven <- hbefTrends$y[SP_BIOTIC, , , ]
coords <- hbefTrends$coords / 1000            # km, as in 60_real_data.R

zs <- function(a) (a - mean(a, na.rm = TRUE)) / sd(a, na.rm = TRUE)
elev_s <- as.numeric(scale(hbefTrends$occ.covs$elev))

## The biotic covariate enters the detection component with the same
## convention as build_det_X(..., na_action = "zero"): a missing record
## is entered as a non-detection. Section 2.3 of the manuscript explains
## why that is immaterial here, and 62_missing_sensitivity.R tests it.
oven_z <- y_oven
oven_z[is.na(oven_z)] <- 0

occ_covs <- list(elev_s = elev_s, elev_s2 = elev_s^2)
det_covs <- list(day = zs(hbefTrends$det.covs$day),
                 tod = zs(hbefTrends$det.covs$tod),
                 oven = oven_z)

dl <- list(y = y_revi, occ.covs = occ_covs, det.covs = det_covs,
           coords = coords)

cat(sprintf("sites %d | seasons %d | visits %d | observed records %d\n",
            dim(y_revi)[1], dim(y_revi)[2], dim(y_revi)[3],
            sum(!is.na(y_revi))))

## ---------------------------------------------------------------------
## Priors. spOccupancy parameterises the spatial range through phi, the
## decay, on 3/range for an exponential covariance. The bounds below span
## practical ranges from roughly a tenth of the domain to several times
## it, which is deliberately wide: the point is not to reproduce the PC
## prior but to avoid a prior that drives the answer.
## ---------------------------------------------------------------------
d <- dist(coords)
lo <- min(d[d > 0]); hi <- max(d)
priors <- list(beta.normal   = list(mean = 0, var = 1),     # occupancy
               alpha.normal  = list(mean = 0, var = 3),     # detection, prec 1/3
               sigma.sq.ig   = c(2, 1),
               phi.unif      = c(3 / hi, 3 / (lo / 2)),
               rho.unif      = c(-1, 1))
cat(sprintf("phi prior spans practical range %.2f to %.2f km\n",
            3 / priors$phi.unif[2], 3 / priors$phi.unif[1]))

tuning <- list(phi = 0.5, rho = 0.5)

fit_one <- function(det_formula, tag) {
  message("\n=== ", tag, " : ", deparse(det_formula), " ===")
  t0 <- proc.time()[3]
  out <- stPGOcc(occ.formula = ~ elev_s + elev_s2,
                 det.formula = det_formula,
                 data = dl, priors = priors, tuning = tuning,
                 cov.model = "exponential", NNGP = TRUE,
                 n.neighbors = N_NEIGHB, ar1 = TRUE,
                 n.batch = N_BATCH, batch.length = BATCH_LEN,
                 n.burn = N_BURN, n.thin = N_THIN, n.chains = N_CHAINS,
                 n.omp.threads = N_THREADS, verbose = TRUE, n.report = 200)
  el <- proc.time()[3] - t0
  saveRDS(out, file.path(OUTDIR, paste0("stPGOcc_", tag, ".rds")))
  message(sprintf("%s finished in %.0f s", tag, el))
  list(fit = out, seconds = el)
}

M1 <- fit_one(~ day + tod,        "M1")
M2 <- fit_one(~ day + tod + oven, "M2")

## ---------------------------------------------------------------------
## Convergence first. A comparison against an unconverged chain says
## nothing, so this is checked before anything is tabulated.
## ---------------------------------------------------------------------
conv <- function(res, tag) {
  f <- res$fit
  ## Split occupancy from detection. The first run showed these converge
  ## at very different rates, and a single worst-case number hides that
  ## the detection block -- which is what this paper is about -- was
  ## already fine.
  blk <- function(rh, es, part)
    data.frame(model = tag, part = part,
               max_rhat = max(rh, na.rm = TRUE),
               min_ess  = min(es, na.rm = TRUE))
  rbind(blk(f$rhat$beta,  f$ESS$beta,  "occupancy"),
        blk(f$rhat$alpha, f$ESS$alpha, "detection"),
        blk(f$rhat$theta, f$ESS$theta, "covariance"))
}
cv <- rbind(conv(M1, "M1"), conv(M2, "M2"))
cv$seconds <- rep(c(M1$seconds, M2$seconds), each = 3)
print(cv, row.names = FALSE, digits = 4)
write.csv(cv, file.path(OUTDIR, "convergence.csv"), row.names = FALSE)

bad <- cv[cv$max_rhat > 1.05 | cv$min_ess < 400, ]
if (nrow(bad)) {
  print(bad, row.names = FALSE, digits = 4)
  warning("the blocks above have not converged well enough to compare. ",
          "Raise N_BATCH or N_BURN. A block that HAS converged can still ",
          "be reported; say which.")
} else {
  message("all blocks converged: R-hat < 1.05 and ESS > 400 throughout")
}

## ---------------------------------------------------------------------
## Side by side with the INLA fits
## ---------------------------------------------------------------------
pull_mcmc <- function(res, tag) {
  b <- summary(res$fit$beta.samples)$statistics
  a <- summary(res$fit$alpha.samples)$statistics
  rbind(
    data.frame(model = tag, part = "occupancy", parameter = rownames(b),
               mean = b[, "Mean"], sd = b[, "SD"]),
    data.frame(model = tag, part = "detection", parameter = rownames(a),
               mean = a[, "Mean"], sd = a[, "SD"]))
}
mcmc_tab <- rbind(pull_mcmc(M1, "M1"), pull_mcmc(M2, "M2"))
write.csv(mcmc_tab, file.path(OUTDIR, "mcmc_estimates.csv"), row.names = FALSE)
print(mcmc_tab, row.names = FALSE, digits = 4)

## The INLA side, if it is on disk. 82_slim_outputs.R writes these; the
## full fits.rds is not needed.
f_inla <- "results/real/m1_m2_summaries.rds"
if (file.exists(f_inla)) {
  inla_sum <- readRDS(f_inla)
  cat("\n--- INLA-SPDE, for comparison ---\n")
  for (nm in names(inla_sum)) {
    cat("\n", toupper(nm), " occupancy:\n", sep = "")
    print(inla_sum[[nm]]$fixed[, c("mean", "sd")], digits = 4)
    cat(toupper(nm), " detection (hyperpar):\n", sep = "")
    print(inla_sum[[nm]]$hyper[, c("mean", "sd")], digits = 4)
  }
  cat("\nMatch the rows by name, not by position: the two packages order\n")
  cat("and label parameters differently, and the detection coefficients\n")
  cat("are hyperparameters in INLA but fixed effects in spOccupancy.\n")
} else {
  message("results/real/m1_m2_summaries.rds not found; run 82_slim_outputs.R")
}

cat("\n--- reporting this in the manuscript ---\n")
cat("Report posterior means and standard deviations side by side, and the\n")
cat("difference in units of the INLA posterior SD, which is the scale on\n")
cat("which a discrepancy matters. Report run times with the settings\n")
cat("above and claim no advantage. State that the spatial models differ\n")
cat("(SPDE-Matern against NNGP-exponential), so a discrepancy is not\n")
cat("necessarily approximation error.\n")
