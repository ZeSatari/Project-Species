## =====================================================================
## 73_prior_sensitivity_field.R -- PC priors on range and sigma
## =====================================================================
## The Limitations section states that prior sensitivity was assessed
## only for the species effects, and not for the field hyperparameters.
## That is the cheapest of the open gaps to close: the empirical fit
## takes a few minutes, so the whole grid is under an hour.
##
## WHAT WE EXPECT, AND WHY REPORTING IT IS NOT AN ADMISSION.
##
## Section 3.1.2 already reports that the practical range sits near the
## limit of what a 7.7 km domain can identify, and that under a fine mesh
## it exceeded the domain diameter and was informed by the prior rather
## than the data. Range and sigma should therefore MOVE with the prior.
## That is the finding already stated, not a new defect.
##
## The question this script answers is a different one: does the estimand
## the paper is about move too? If alpha_2 and the other detection
## coefficients are stable while the field hyperparameters shift, the
## conclusion is that the detection inference does not rest on the prior
## for the field. If alpha_2 moves, that is a real problem and the paper
## must say so.
##
## DESIGN. One factor at a time from the baseline, not a full grid: with
## two hyperparameters and an already-flagged identifiability problem, a
## grid would mostly show that the two trade off against each other,
## which is not in question.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/prior_field"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "prior_field")

## ---------------------------------------------------------------------
## Data, exactly as 60_real_data.R prepares it
## ---------------------------------------------------------------------
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
  time    = rep(seq_len(N_YEAR), each  = N_SITE),
  x       = rep(coords[, 1], times = N_YEAR),
  y       = rep(coords[, 2], times = N_YEAR),
  elev_s  = rep(elev_s,      times = N_YEAR),
  elev_s2 = rep(elev_s^2,    times = N_YEAR),
  Int_occ = 1)

Y_revi <- stack_species("REVI")
Y_oven <- stack_species("OVEN")
X_M1   <- build_det_X(list(stack_det_cov("day"), stack_det_cov("tod")), K = K)
X_M2   <- build_det_X(list(stack_det_cov("day"), stack_det_cov("tod"),
                           Y_oven), K = K)

bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
message("mesh: ", mesh$n, " nodes (expect 283)")

## ---------------------------------------------------------------------
## The grid. Baseline is PRIORS$spde_real; each row changes one thing.
## The range quantiles span a tenth of the domain to its full diameter,
## which brackets the identifiability problem rather than avoiding it.
## ---------------------------------------------------------------------
base_range <- PRIORS$spde_real$range
base_sigma <- PRIORS$spde_real$sigma
cat("baseline: Pr(range <", base_range[1], "km) =", base_range[2],
    "| Pr(sigma >", base_sigma[1], ") =", base_sigma[2], "\n")

GRID <- list(
  list(tag = "baseline",     range = base_range,  sigma = base_sigma),
  list(tag = "range_short",  range = c(1,  0.7),  sigma = base_sigma),
  list(tag = "range_long",   range = c(10, 0.7),  sigma = base_sigma),
  list(tag = "range_vague",  range = c(5,  0.3),  sigma = base_sigma),
  list(tag = "sigma_tight",  range = base_range,  sigma = c(1, 0.1)),
  list(tag = "sigma_loose",  range = base_range,  sigma = c(3, 0.5))
)

fit_pair <- function(cfg) {
  spde <- INLA::inla.spde2.pcmatern(mesh,
                                    prior.range = cfg$range,
                                    prior.sigma = cfg$sigma,
                                    constr = TRUE)
  stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                       time = Xocc$time, n_time = N_YEAR)
  ar1 <- list(rho = list(prior = "pc.cor0", param = PRIORS$ar1$rho$param))
  base <- c(as.list(Xocc), stf$idx)

  f <- as.formula(paste(
    "inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +",
    "f(spatialfield, model = spde, A.local = stf$A,",
    "  group = spatialfield.group,",
    "  control.group = list(model = 'ar1', hyper = ar1))"))
  environment(f) <- environment()

  one <- function(Xdet, tag_m) {
    t0 <- proc.time()[3]
    fit <- tryCatch(fit_occupancy(f, c(base, list(Y = Y_revi, X = Xdet)),
                                  Xdet, compute_cv = FALSE),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      message("  FAILED ", cfg$tag, " / ", tag_m, ": ", conditionMessage(fit))
      return(NULL)
    }
    hp <- fit$summary.hyperpar
    ## Detection coefficients are hyperparameters, indexed by name so a
    ## change in their number cannot silently shift the columns.
    bet <- hp[grep("^beta\\[", rownames(hp)), , drop = FALSE]
    gr <- function(p) { i <- grep(p, rownames(hp))
                        if (length(i)) hp[i[1], "mean"] else NA_real_ }
    data.frame(
      prior = cfg$tag, model = tag_m,
      pr_range_q = cfg$range[1], pr_range_p = cfg$range[2],
      pr_sigma_q = cfg$sigma[1], pr_sigma_p = cfg$sigma[2],
      det_int  = bet[1, "mean"], det_int_sd = bet[1, "sd"],
      det_day  = bet[2, "mean"],
      det_tod  = bet[3, "mean"],
      alpha2   = if (nrow(bet) >= 4) bet[4, "mean"]  else NA_real_,
      alpha2_sd= if (nrow(bet) >= 4) bet[4, "sd"]    else NA_real_,
      beta0    = fit$summary.fixed["Int_occ", "mean"],
      elev     = fit$summary.fixed["elev_s",  "mean"],
      elev2    = fit$summary.fixed["elev_s2", "mean"],
      range    = gr("^Range for"), sigma = gr("^Stdev for"),
      rho      = gr("^GroupRho for"),
      waic     = fit$waic$waic,
      seconds  = round(proc.time()[3] - t0, 1))
  }
  rbind(one(X_M1, "M1"), one(X_M2, "M2"))
}

res <- list()
for (k in seq_along(GRID)) {
  message(sprintf("[%d/%d] %s", k, length(GRID), GRID[[k]]$tag))
  res[[k]] <- fit_pair(GRID[[k]])
  saveRDS(do.call(rbind, res), file.path(OUTDIR, "prior_field.rds"))
}
tab <- do.call(rbind, res)
write.csv(tab, file.path(OUTDIR, "prior_field.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## Report. The comparison that matters is the spread of alpha_2 relative
## to its own posterior standard deviation: a prior that moves a
## parameter by a fraction of its posterior SD has not driven the answer.
## ---------------------------------------------------------------------
print(tab[, c("prior", "model", "alpha2", "det_int", "beta0",
              "elev", "range", "sigma", "waic")],
      row.names = FALSE, digits = 4)

m2 <- tab[tab$model == "M2" & !is.na(tab$alpha2), ]
if (nrow(m2)) {
  sd_ref <- m2$alpha2_sd[m2$prior == "baseline"]
  cat(sprintf("\nalpha_2 across priors: %.4f to %.4f (spread %.4f)\n",
              min(m2$alpha2), max(m2$alpha2), diff(range(m2$alpha2))))
  cat(sprintf("  that spread is %.2f of one posterior SD (%.4f)\n",
              diff(range(m2$alpha2)) / sd_ref, sd_ref))
}
for (v in c("range", "sigma", "beta0", "elev")) {
  x <- tab[[v]][tab$model == "M2"]
  cat(sprintf("%-7s across priors: %.3f to %.3f\n", v, min(x), max(x)))
}

cat("\n--- how to read this ---\n")
cat("The field hyperparameters are expected to move: Section 3.1.2 already\n")
cat("reports that the range is at the limit of what this domain identifies.\n")
cat("What matters is whether alpha_2 moves with them. A spread well under\n")
cat("one posterior SD means the detection inference does not rest on the\n")
cat("prior for the field, and the Limitations paragraph can say so.\n")
