## =====================================================================
## 75_field_prior_msom.R -- does the field prior drive the fit?
## =====================================================================
## PaperB, Section 3.1. The companion study varied the field priors for
## the single-species model; this does the same for the hierarchical one
## (the reported model M1: pooled detection, common elevation response).
##
## Seven specifications: five for the practical range, holding the sigma
## prior at its baseline, and two for the marginal standard deviation,
## holding the range prior at its baseline. The baseline is the reported
## fit, PRIORS$spde_real.
##
## Reported per specification: range, marginal SD and autoregressive
## correlation of the field; the community intercept and the elevation
## coefficients; and the largest change in any species-level intercept
## beta0 + gamma_i relative to the baseline fit. Species intercepts come
## from the marginal summaries, so no posterior sampling is needed and
## the script costs one fit per specification (about 5 min each here).
##
## The question is not whether the hyperparameters move -- a PC prior on
## the range is informative by construction when the data are of limited
## extent -- but whether the quantities the paper reports move with them.
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/prior_sens"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "field_prior")

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP     <- length(SP_NAMES)
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K_VIS    <- dim(hbefTrends$y)[4]
N_ROW    <- N_SITE * N_YEAR

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
Xday    <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")

mesh <- fmesher::fm_mesh_2d(
  boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
  max.edge = c(0.35, 1.2), min.angle = 21,
  offset = c(0.05, 1.5), cutoff = 0.35)
message("mesh: ", mesh$n, " nodes")

Y_h <- do.call(rbind, lapply(SP_NAMES, stack_species))
Xd  <- build_det_X(list(do.call(rbind, rep(list(Xday), N_SP)),
                        do.call(rbind, rep(list(Xtod), N_SP))), K = K_VIS)
Xocc <- data.frame(
  Int_occ = 1,
  elev_s  = rep(elev_s[site_id], times = N_SP),
  elev_s2 = rep(elev_s[site_id]^2, times = N_SP),
  species_id = rep(seq_len(N_SP), each = N_ROW),
  x = rep(coords[site_id, 1], times = N_SP),
  y = rep(coords[site_id, 2], times = N_SP),
  time = rep(time_id, times = N_SP))

## ---------------------------------------------------------------------
## Specifications. prior.range = c(r0, p): P(range < r0) = p, in km.
## prior.sigma = c(s0, p): P(sigma > s0) = p.
## ---------------------------------------------------------------------
base <- PRIORS$spde_real
cat("baseline field prior:\n"); print(base)

SPECS <- list(
  `range<1km`   = list(range = c(1,  0.7), sigma = base$sigma),
  `range<2km`   = list(range = c(2,  0.7), sigma = base$sigma),
  baseline      = list(range = base$range, sigma = base$sigma),
  `range<10km`  = list(range = c(10, 0.7), sigma = base$sigma),
  `range<20km`  = list(range = c(20, 0.7), sigma = base$sigma),
  `sigma>0.5`   = list(range = base$range, sigma = c(0.5, 0.5)),
  `sigma>2`     = list(range = base$range, sigma = c(2.0, 0.5)))

fit_one <- function(spec) {
  spde <- INLA::inla.spde2.pcmatern(mesh,
            prior.range = spec$range, prior.sigma = spec$sigma, constr = TRUE)
  stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                       time = Xocc$time, n_time = N_YEAR)
  f <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
    f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))
  fit_occupancy(f, c(as.list(Xocc), stf$idx, list(Y = Y_h, X = Xd)),
                Xd, compute_cv = FALSE)
}

grab <- function(fit, tag) {
  hp <- fit$summary.hyperpar
  g  <- function(pat, col = "mean") hp[grep(pat, rownames(hp))[1], col]
  fx <- fit$summary.fixed
  sp_int <- fx["Int_occ", "mean"] + fit$summary.random$species_id$mean
  list(row = data.frame(
         spec      = tag,
         range     = g("^Range"),
         range_lo  = g("^Range", "0.025quant"), range_hi = g("^Range", "0.975quant"),
         sigma     = g("^Stdev"),
         rho       = g("^GroupRho"),
         beta0     = fx["Int_occ", "mean"],
         beta0_sd  = fx["Int_occ", "sd"],
         elev      = fx["elev_s", "mean"],
         elev2     = fx["elev_s2", "mean"],
         det_int   = g("beta\\[0\\]")),
       sp_int = sp_int)
}

res <- list(); ints <- list()
for (tag in names(SPECS)) {
  message("\n=== ", tag, " ===")
  t0 <- Sys.time()
  fit <- fit_one(SPECS[[tag]])
  message(sprintf("  %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  gg <- grab(fit, tag)
  res[[tag]] <- gg$row; ints[[tag]] <- gg$sp_int
}
tab <- do.call(rbind, res)

## largest change in any species-level intercept, against the baseline
b <- ints[["baseline"]]
tab$max_d_species <- vapply(ints, \(v) max(abs(v - b)), numeric(1))
tab$max_d_beta0_sd <- abs(tab$beta0 - tab$beta0[tab$spec == "baseline"]) /
                      tab$beta0_sd[tab$spec == "baseline"]

cat("\n", strrep("=", 100), "\n", sep = "")
print(tab, row.names = FALSE, digits = 3)
cat("\nmax_d_species: largest absolute change in beta0 + gamma_i vs the baseline fit\n")
cat("max_d_beta0_sd: change in the community intercept, in baseline posterior SDs\n")
cat(sprintf("\nrange over all specifications: %.2f to %.2f km\n",
            min(tab$range), max(tab$range)))
cat(sprintf("largest change in any species intercept: %.3f on the logit scale\n",
            max(tab$max_d_species)))

sp_tab <- data.frame(species = SP_NAMES, do.call(cbind, ints))
write.csv(tab,    file.path(OUTDIR, "field_prior_sensitivity.csv"), row.names = FALSE)
write.csv(sp_tab, file.path(OUTDIR, "field_prior_species_intercepts.csv"), row.names = FALSE)
cat("\nwritten to ", OUTDIR, "\n", sep = "")
