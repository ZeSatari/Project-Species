## =====================================================================
## 74_strategy_check_msom.R -- Gaussian vs simplified vs full Laplace
## =====================================================================
## PaperB, Section 3.3. The occupancy log-likelihood is not log-concave,
## so the accuracy of the Laplace approximation cannot be assumed. The
## companion study made this comparison for the single-species model;
## this does it for the hierarchical one (the reported model: pooled
## detection, common elevation response), refitting it under the three
## strategies with everything else held fixed.
##
## Reported: posterior means and standard deviations of the fixed
## effects, the hyperparameters and the twelve species-level intercepts
## beta0 + gamma_i, and the discrepancy of each strategy from the full
## Laplace fit measured in units of the full-Laplace posterior standard
## deviation -- the scale on which such differences matter, since a shift
## small against the posterior spread changes no conclusion.
##
## Three fits, roughly 5 min each for the simplified strategy; the full
## Laplace is slower, so allow half an hour in total.
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/strategy"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "strategy")

STRATEGIES <- c("gaussian", "simplified.laplace", "laplace")
REFERENCE  <- "laplace"

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
dat <- c(as.list(Xocc), stf$idx, list(Y = Y_h, X = Xd))

## ---------------------------------------------------------------------
## Fit under each strategy
## ---------------------------------------------------------------------
fits <- list(); times <- numeric(0)
for (st in STRATEGIES) {
  message("\n=== ", st, " ===")
  t0 <- Sys.time()
  fits[[st]] <- fit_occupancy(f, dat, Xd, strategy = st, compute_cv = FALSE)
  times[st] <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("  %.1f min", times[st]))
}
saveRDS(fits, file.path(OUTDIR, "strategy_fits.rds"))

## ---------------------------------------------------------------------
## Collect every reported quantity: mean and posterior SD
## ---------------------------------------------------------------------
collect <- function(fit) {
  fx <- fit$summary.fixed[, c("mean", "sd"), drop = FALSE]
  hp <- fit$summary.hyperpar[, c("mean", "sd"), drop = FALSE]
  sr <- fit$summary.random$species_id
  sp <- data.frame(
    mean = fit$summary.fixed["Int_occ", "mean"] + sr$mean,
    ## SD of beta0 + gamma_i is not the SD of gamma_i alone; the marginal
    ## summaries do not give their covariance, so we use the wider of the
    ## two as a conservative scale for the comparison below.
    sd   = pmax(sr$sd, fit$summary.fixed["Int_occ", "sd"]),
    row.names = paste0("psi_int_", SP_NAMES))
  rbind(fx, hp, sp)
}
tabs <- lapply(fits, collect)
stopifnot(all(vapply(tabs, \(t) identical(rownames(t), rownames(tabs[[1]])), TRUE)))

ref <- tabs[[REFERENCE]]
out <- data.frame(parameter = rownames(ref))
for (st in STRATEGIES) {
  out[[paste0("mean_", st)]] <- tabs[[st]]$mean
  out[[paste0("d_", st)]]    <- (tabs[[st]]$mean - ref$mean) / ref$sd
}
out$sd_ref <- ref$sd

cat("\n", strrep("=", 100), "\n", sep = "")
cat("Posterior means, and differences from the ", REFERENCE,
    " fit in units of its posterior SD\n\n", sep = "")
print(out, row.names = FALSE, digits = 3)

cat("\n--- largest discrepancy from the full Laplace fit, by strategy ---\n")
for (st in setdiff(STRATEGIES, REFERENCE)) {
  d <- abs(out[[paste0("d_", st)]])
  i <- which.max(d)
  cat(sprintf("  %-20s max %.3f posterior SD (at %s); median %.3f\n",
              st, d[i], out$parameter[i], median(d)))
}
cat("\n--- species-level intercepts only ---\n")
sp_rows <- grep("^psi_int_", out$parameter)
for (st in setdiff(STRATEGIES, REFERENCE))
  cat(sprintf("  %-20s max change %.4f on the logit scale, %.3f posterior SD\n", st,
              max(abs(out[[paste0("mean_", st)]][sp_rows] -
                      out[[paste0("mean_", REFERENCE)]][sp_rows])),
              max(abs(out[[paste0("d_", st)]][sp_rows]))))

cat("\n--- run time (minutes) ---\n"); print(round(times, 1))

write.csv(out, file.path(OUTDIR, "strategy_comparison.csv"), row.names = FALSE)
cat("\nwritten to ", OUTDIR, "\n", sep = "")
