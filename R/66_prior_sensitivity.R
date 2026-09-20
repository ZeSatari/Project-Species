## =====================================================================
## 66_prior_sensitivity.R -- does the species-effect prior drive the fit?
## =====================================================================
## The empirical hierarchical model estimates sigma_gamma at 1.608, with
## a 95% credible interval of (1.138, 2.129). The prior states
## P(sigma_gamma > 1) = 0.01. The whole interval therefore sits in the
## upper 1% tail of the prior, which is a prior-data conflict: the data
## are pulling hard against a prior that considers this region nearly
## impossible.
##
## That is not automatically a problem -- with 40,284 observation rows
## the data may simply dominate -- but it is untested, and the
## heterogeneous-prevalence simulation showed that conflict of this kind
## can produce estimates that behave strangely. This script tests it.
##
## Four priors on sigma_gamma, everything else identical:
##   P(sigma > 1) = 0.01   as used throughout (the strongest shrinkage)
##   P(sigma > 2) = 0.01
##   P(sigma > 3) = 0.01
##   P(sigma > 5) = 0.01   the weakest
##
## If sigma_gamma and the species-level intercepts are stable across
## these, the prior is not driving the result and the paper can say so.
## If they move materially, the reported values are a compromise between
## prior and data and must be described that way.
##
## Runtime: four fits of the full multi-species model, roughly 20 minutes.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/prior"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR)

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

PRIOR_SET <- list(
  "P(sigma>1)=0.01" = c(1, 0.01),
  "P(sigma>2)=0.01" = c(2, 0.01),
  "P(sigma>3)=0.01" = c(3, 0.01),
  "P(sigma>5)=0.01" = c(5, 0.01)
)

fit_one <- function(param, tag) {
  hyper <- list(prec = list(prior = "pc.prec", param = param))
  f <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
    f(species_id, model = "iid", constr = TRUE, hyper = hyper) +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))
  message("fitting ", tag)
  fit <- tryCatch(
    fit_occupancy(f, c(as.list(Xocc), stf$idx, list(Y = Y_all, X = Xd)),
                  Xd, compute_cv = FALSE),
    error = function(e) e)
  if (inherits(fit, "error")) {
    message("  FAILED: ", conditionMessage(fit)); return(NULL)
  }

  ## sigma_gamma from the transformed marginal, not 1/sqrt(mean precision):
  ## the mean of a function is not the function of the mean.
  pm <- fit$marginals.hyperpar[["Precision for species_id"]]
  sg <- INLA::inla.zmarginal(
          INLA::inla.tmarginal(function(x) 1 / sqrt(x), pm), silent = TRUE)

  hp <- fit$summary.hyperpar
  gr <- function(p) { i <- grep(p, rownames(hp))
                      if (length(i) == 1) hp[i, "mean"] else NA_real_ }
  ## species_intercepts() takes the NUMBER of species, not their names,
  ## and returns them indexed 1..N in the order of the species_id factor.
  si <- species_intercepts(fit, N_SP)

  list(summary = data.frame(
         prior = tag,
         sigma_gamma = sg$mean, sg_lwr = sg$quant0.025, sg_upr = sg$quant0.975,
         beta0 = fit$summary.fixed["Int_occ", "mean"],
         elev  = fit$summary.fixed["elev_s", "mean"],
         range = gr("^Range for"), field_sd = gr("^Stdev for"),
         rho   = gr("^GroupRho for"), waic = fit$waic$waic),
       intercepts = setNames(si$mean, SP_NAMES[si$species]))
}

res <- list(); ints <- list()
for (tag in names(PRIOR_SET)) {
  r <- fit_one(PRIOR_SET[[tag]], tag)
  if (!is.null(r)) { res[[tag]] <- r$summary; ints[[tag]] <- r$intercepts }
}

out <- do.call(rbind, res); rownames(out) <- NULL
cat("\n", strrep("=", 88), "\n", sep = "")
print(out, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------
## Species-level intercepts are the quantity the paper reports, so they
## matter more than the hyperparameter itself.
## ---------------------------------------------------------------------
if (length(ints) >= 2) {
  M <- do.call(cbind, ints)
  cat("\nspecies-level intercepts (beta0 + gamma_i) under each prior:\n")
  print(round(M, 3))
  cat("\nmaximum shift across priors, by species:\n")
  print(round(apply(M, 1, function(z) diff(range(z))), 3))
  cat(sprintf("\nlargest shift for any species: %.3f on the logit scale\n",
              max(apply(M, 1, function(z) diff(range(z))))))
  base <- M[, 1]
  cat("correlation with the reported fit:\n")
  print(round(apply(M, 2, function(z) cor(z, base)), 4))
}

cat("\nReading this. If sigma_gamma rises steadily as the prior is weakened,\n")
cat("the reported value is a compromise and the paper must say so. If it\n")
cat("settles, the data are determining it and the prior is not binding.\n")
cat("The species-level intercepts matter more than sigma_gamma itself,\n")
cat("since those are what the paper tabulates.\n")

write.csv(out, file.path(OUTDIR, "prior_sensitivity.csv"), row.names = FALSE)
if (length(ints)) saveRDS(ints, file.path(OUTDIR, "prior_intercepts.rds"))
