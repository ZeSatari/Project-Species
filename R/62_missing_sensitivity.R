## =====================================================================
## 62_missing_sensitivity.R -- how much does the NA -> 0 convention matter?
## =====================================================================
## Records of the co-occurring species are missing for 12.9% of surveys and
## are entered as non-detections. The manuscript states that this attenuates
## alpha_2 toward zero, but until now offered no measurement of how much.
##
## Three treatments of the same data:
##
##   zero      missing -> non-detection. The convention used throughout.
##             Introduces measurement error in the covariate; classical
##             attenuation pulls alpha_2 toward zero.
##
##   complete  survey OCCASIONS with any missing y_B are dropped from the
##             detection design entirely, by setting the focal species'
##             response to NA for those visits. No covariate error, but a
##             smaller sample, and valid only if missingness is unrelated
##             to detectability.
##
##   indicator missing -> 0, plus a separate indicator column flagging
##             which records were imputed. The coefficient on y_B is then
##             estimated from observed records alone, with the imputed ones
##             absorbed by the indicator.
##
## Agreement across the three means the convention is not driving the
## result. Disagreement means it must be reported as a caveat on alpha_2.
##
## Reads results/real/fits.rds only for comparison; refits three models.
## Runtime roughly 10 minutes with the coarsened mesh.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/real"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K        <- dim(hbefTrends$y)[4]

stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
}

coords <- hbefTrends$coords / 1000
elev_s <- as.numeric(scale(hbefTrends$occ.covs$elev))

Xocc <- data.frame(
  time = rep(seq_len(N_YEAR), each = N_SITE),
  x = rep(coords[, 1], times = N_YEAR),
  y = rep(coords[, 2], times = N_YEAR),
  elev_s = rep(elev_s, times = N_YEAR),
  elev_s2 = rep(elev_s^2, times = N_YEAR),
  Int_occ = 1
)
N_ROW <- nrow(Xocc)

Y_revi   <- stack_species("REVI")
Y_oven   <- stack_species("OVEN")
Xdet_day <- stack_det_cov("day")
Xdet_tod <- stack_det_cov("tod")

miss <- is.na(Y_oven)
cat(sprintf("Missing OVEN records: %d of %d visit-level records (%.1f%%)\n",
            sum(miss), length(miss), 100 * mean(miss)))
cat(sprintf("Site-years with at least one missing OVEN record: %d of %d (%.1f%%)\n",
            sum(rowSums(miss) > 0), nrow(miss),
            100 * mean(rowSums(miss) > 0)))

## Is missingness related to whether the focal species was detected? If it
## is, dropping incomplete occasions is not innocuous either.
det_focal <- !is.na(Y_revi) & Y_revi == 1
tab <- table(oven_missing = as.vector(miss),
             revi_detected = as.vector(det_focal))
print(tab)
cat(sprintf("REVI detection rate where OVEN observed: %.3f | where missing: %.3f\n",
            mean(det_focal[!miss]), mean(det_focal[miss])))

bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = PRIORS$spde_real$range,
                                  prior.sigma = PRIORS$spde_real$sigma,
                                  constr = TRUE)
stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                     time = Xocc$time, n_time = N_YEAR)
base <- c(as.list(Xocc), stf$idx)

f_one <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
  f(spatialfield, model = spde, A.local = stf$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

## ---------------------------------------------------------------------
## Three treatments.
## ---------------------------------------------------------------------
fit_one <- function(label, Yfocal, det_list) {
  Xd <- build_det_X(det_list, K = K, na_action = "zero")
  message("fitting: ", label)
  fit <- fit_occupancy(f_one, c(base, list(Y = Yfocal, X = Xd)), Xd,
                       compute_cv = FALSE)
  dr <- det_draws(fit)
  ## the biotic term is the last covariate supplied
  b <- dr[, ncol(dr)]
  data.frame(
    treatment = label,
    n_used    = sum(!is.na(Yfocal)),
    alpha_biotic = mean(b),
    sd        = sd(b),
    lwr       = quantile(b, 0.025),
    upr       = quantile(b, 0.975),
    waic      = fit$waic$waic,
    row.names = NULL
  )
}

## (a) the convention used in the manuscript
res_zero <- fit_one("missing -> 0", Y_revi, list(Xdet_day, Xdet_tod, Y_oven))

## (b) drop visits whose OVEN record is missing, by making the focal
##     response missing there too
Y_revi_cc <- Y_revi
Y_revi_cc[miss] <- NA
res_cc <- fit_one("incomplete visits dropped", Y_revi_cc,
                  list(Xdet_day, Xdet_tod, Y_oven))

## (c) missing -> 0, with an indicator for imputed records
IND <- matrix(as.numeric(miss), nrow = nrow(miss))
res_ind <- fit_one("missing -> 0 + indicator", Y_revi,
                   list(Xdet_day, Xdet_tod, IND, Y_oven))

out <- rbind(res_zero, res_cc, res_ind)
print(out, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------
## Read the spread in units of the primary analysis's posterior SD. A
## shift of 0.02 means nothing if sd = 0.3 and a great deal if sd = 0.01.
## ---------------------------------------------------------------------
ref <- out$alpha_biotic[1]; ref_sd <- out$sd[1]
out$shift_in_sd <- abs(out$alpha_biotic - ref) / ref_sd
cat("\nshift from the primary analysis, in its posterior SD units:\n")
print(out[, c("treatment", "alpha_biotic", "shift_in_sd")],
      row.names = FALSE, digits = 3)

worst <- max(out$shift_in_sd)
cat(sprintf("\nLargest shift: %.2f posterior SD -> %s\n", worst,
            if (worst < 0.25) "the convention does not drive the estimate."
            else if (worst < 1) "modest sensitivity; report it."
            else "the estimate depends on the convention; report prominently."))

write.csv(out, file.path(OUTDIR, "missing_sensitivity.csv"), row.names = FALSE)
