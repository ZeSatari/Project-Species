## =====================================================================
## 67_sigma_gamma_check.R -- is sigma_gamma reproducible?
## =====================================================================
## Two fits of what appears to be the same model returned different
## values for the species-effect standard deviation:
##
##   60_real_data.R          sigma_gamma = 1.608  (precision 0.4185)
##   66_prior_sensitivity.R  sigma_gamma = 1.171
##
## Everything else agreed. WAIC differed by 0.6 in 77,923; beta_0 by
## 0.004; the range by 0.02; and the species-level intercepts
## beta_0 + gamma_i by less than 0.01 for every species. So the fits are
## the same in every respect that the paper reports, except this one
## hyperparameter.
##
## Two explanations, with different consequences:
##
##   (a) sigma_gamma is weakly identified. With twelve species the
##       likelihood is flat in that direction and the optimiser settles
##       wherever it starts. The reported interval (1.138, 2.129) already
##       covers both values, so nothing is wrong, but the point estimate
##       should not be leaned on.
##
##   (b) the two scripts differ in some way not yet found, in which case
##       one of them is not the model the paper describes.
##
## Refitting the identical specification several times distinguishes
## them: under (a) the spread across refits should be comparable to the
## 1.17-1.61 gap; under (b) the refits should agree with each other and
## the discrepancy lies between scripts.
##
## Runtime: N_REFIT fits, roughly 5 minutes each.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/prior"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_REFIT <- 3

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP   <- length(SP_NAMES)
N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
K_VIS  <- dim(hbefTrends$y)[4]
N_ROW  <- N_SITE * N_YEAR

stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  M <- (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
  M[is.na(M)] <- 0
  M
}

coords  <- hbefTrends$coords / 1000
elev_s  <- as.numeric(scale(hbefTrends$occ.covs$elev))
site_id <- rep(seq_len(N_SITE), times = N_YEAR)
time_id <- rep(seq_len(N_YEAR), each = N_SITE)

Y_all <- do.call(rbind, lapply(SP_NAMES, stack_species))
Xd    <- build_det_X(list(do.call(rbind, rep(list(stack_det_cov("day")), N_SP)),
                          do.call(rbind, rep(list(stack_det_cov("tod")), N_SP))),
                     K = K_VIS)
Xocc <- data.frame(
  Int_occ = 1,
  elev_s  = rep(elev_s[site_id], times = N_SP),
  elev_s2 = rep(elev_s[site_id]^2, times = N_SP),
  species_id = rep(seq_len(N_SP), each = N_ROW),
  x = rep(coords[site_id, 1], times = N_SP),
  y = rep(coords[site_id, 2], times = N_SP),
  time = rep(time_id, times = N_SP))

mesh <- fmesher::fm_mesh_2d(
  boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
  max.edge = c(0.35, 1.2), min.angle = 21,
  offset = c(0.05, 1.5), cutoff = 0.35)
spde <- INLA::inla.spde2.pcmatern(mesh,
          prior.range = PRIORS$spde_real$range,
          prior.sigma = PRIORS$spde_real$sigma, constr = TRUE)
stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                     time = Xocc$time, n_time = N_YEAR)

f <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
  f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
  f(spatialfield, model = spde, A.local = stf$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

rows <- list()
for (k in seq_len(N_REFIT)) {
  message("refit ", k, " of ", N_REFIT)
  fit <- fit_occupancy(f, c(as.list(Xocc), stf$idx, list(Y = Y_all, X = Xd)),
                       Xd, compute_cv = FALSE)
  hp  <- fit$summary.hyperpar
  tau <- hp["Precision for species_id", "mean"]
  pm  <- fit$marginals.hyperpar[["Precision for species_id"]]
  sg  <- INLA::inla.zmarginal(
           INLA::inla.tmarginal(function(x) 1 / sqrt(x), pm), silent = TRUE)
  gr  <- function(p) { i <- grep(p, rownames(hp))
                       if (length(i) == 1) hp[i, "mean"] else NA_real_ }
  rows[[k]] <- data.frame(
    refit = k,
    precision = tau,
    one_over_sqrt_mean = 1 / sqrt(tau),   # the wrong transform, for contrast
    sigma_gamma = sg$mean,                # E[1/sqrt(tau)], the right one
    sg_lwr = sg$quant0.025, sg_upr = sg$quant0.975,
    beta0 = fit$summary.fixed["Int_occ", "mean"],
    range = gr("^Range for"), field_sd = gr("^Stdev for"),
    rho = gr("^GroupRho for"), waic = fit$waic$waic)
}

out <- do.call(rbind, rows)
cat("\n", strrep("=", 80), "\n", sep = "")
print(out, row.names = FALSE, digits = 5)

cat("\nspread across refits (max - min):\n")
for (v in c("precision", "sigma_gamma", "beta0", "range", "rho", "waic")) {
  cat(sprintf("  %-12s %.5f\n", v, diff(range(out[[v]]))))
}

cat("\nvalues to compare against:\n")
cat("  60_real_data.R         sigma_gamma 1.608, precision 0.4185\n")
cat("  66_prior_sensitivity.R sigma_gamma 1.171\n\n")
sp_range <- diff(range(out$sigma_gamma))
if (sp_range > 0.2) {
  cat("The refits disagree among themselves by ", signif(sp_range, 3),
      ", comparable to the\ngap between scripts: sigma_gamma is weakly ",
      "identified and its point estimate\nshould not be relied on.\n", sep = "")
} else {
  cat("The refits agree with each other to within ", signif(sp_range, 3),
      ". The gap between the\ntwo scripts is therefore NOT run-to-run ",
      "variation, and something differs\nbetween them that has not been ",
      "found. Do not report either value until it is.\n", sep = "")
}

write.csv(out, file.path(OUTDIR, "sigma_gamma_refits.csv"), row.names = FALSE)
