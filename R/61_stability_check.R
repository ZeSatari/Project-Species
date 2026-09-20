## =====================================================================
## 61_stability_check.R -- is the empirical single-species fit stable?
## =====================================================================
## The same model converged on one run and crashed on the next with no
## change to its inputs. Either the fit is genuinely unstable -- which
## would be a reportable property of this likelihood, not a nuisance --
## or something in the environment differed. This script settles it by
## refitting the identical model several times and recording outcomes.
##
## Self-contained: run with source(), not line by line.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

N_TRIES <- 3

## ---------------------------------------------------------------------
## Rebuild exactly the objects Part A of 60_real_data.R uses.
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

Xocc_one <- data.frame(
  site = rep(seq_len(N_SITE), times = N_YEAR),
  time = rep(seq_len(N_YEAR), each = N_SITE),
  x = rep(coords[, 1], times = N_YEAR),
  y = rep(coords[, 2], times = N_YEAR),
  elev_s  = rep(elev_s,   times = N_YEAR),
  elev_s2 = rep(elev_s^2, times = N_YEAR),
  Int_occ = 1
)

Y_revi   <- stack_species("REVI")
Y_oven   <- stack_species("OVEN")
Xdet_day <- stack_det_cov("day")
Xdet_tod <- stack_det_cov("tod")

## ---------------------------------------------------------------------
## Mesh resolution.
##
## The original specification -- max.edge = c(0.1, 0.7), cutoff = 0.12 --
## produced 1724 nodes, so the latent field held 1724 * 9 = 15,516
## elements against 3,357 observation rows: 4.6 field elements per
## observation. Repeated fits of the IDENTICAL model then returned
## community intercepts spanning 0.183 on the logit scale, run times
## varying by a factor of six (1900 s to 13,100 s), and outright failures
## in 20-80% of attempts depending on the model and strategy. The
## optimiser was landing in different local optima.
##
## Coarsening to max.edge = c(0.35, 1.2), cutoff = 0.35 gives 283 nodes,
## hence 2,547 field elements and a ratio of 0.76 -- comparable to the
## simulation setting, which was stable. This is not a loss of quality:
## a mesh finer than the data can support promises a resolution that
## 373 sites do not contain, and buys it with irreproducibility.
##
## Verify stability with R/61_stability_check.R after any change here.
## ---------------------------------------------------------------------
bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = PRIORS$spde_real$range,
                                  prior.sigma = PRIORS$spde_real$sigma,
                                  constr = TRUE)
stf_one <- make_st_field(mesh, spde, Xocc_one[, c("x", "y")],
                         time = Xocc_one$time, n_time = N_YEAR)

X_M1 <- build_det_X(list(Xdet_day, Xdet_tod), K = K)
X_M2 <- build_det_X(list(Xdet_day, Xdet_tod, Y_oven), K = K)
base_one <- c(as.list(Xocc_one), stf_one$idx)

f_one <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
  f(spatialfield, model = spde, A.local = stf_one$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

message("mesh nodes: ", mesh$n,
        " | field elements: ", stf_one$n_field,
        " | observation rows: ", nrow(Xocc_one))

## ---------------------------------------------------------------------
## Repeat the identical fit and record what happens.
## ---------------------------------------------------------------------
try_fit <- function(label, det_X, strategy = "simplified.laplace") {
  out <- vector("list", N_TRIES)
  for (k in seq_len(N_TRIES)) {
    r <- tryCatch(
      fit_occupancy(f_one, c(base_one, list(Y = Y_revi, X = det_X)),
                    det_X, strategy = strategy, compute_cv = FALSE),
      error = function(e) e
    )
    out[[k]] <- if (inherits(r, "error")) {
      data.frame(model = label, strategy = strategy, try = k, ok = FALSE,
                 beta0 = NA_real_, secs = NA_real_)
    } else {
      data.frame(model = label, strategy = strategy, try = k, ok = TRUE,
                 beta0 = r$summary.fixed["Int_occ", "mean"],
                 secs  = r$elapsed)
    }
    cat(sprintf("  %s / %s / try %d: %s\n", label, strategy, k,
                if (out[[k]]$ok) sprintf("ok, beta0 = %.4f", out[[k]]$beta0)
                else "FAILED"))
  }
  do.call(rbind, out)
}

res <- rbind(
  try_fit("M1", X_M1),
  try_fit("M2", X_M2)
)

## If the default strategy is unreliable, does a more forgiving one help?
if (any(!res$ok)) {
  message("\nFailures observed; retrying under strategy = 'gaussian'.")
  res <- rbind(res,
               try_fit("M1", X_M1, strategy = "gaussian"),
               try_fit("M2", X_M2, strategy = "gaussian"))
}

cat("\n", strrep("=", 60), "\n", sep = "")
print(res, row.names = FALSE, digits = 5)

summ <- aggregate(ok ~ model + strategy, res, mean)
names(summ)[3] <- "success_rate"
cat("\n")
print(summ, row.names = FALSE)

## Identical inputs should give identical output. Any spread across
## successful tries means the fit is not reproducible, which matters
## more than the failures themselves: a number that changes between runs
## cannot go in a table.
for (m in unique(res$model)) {
  b <- res$beta0[res$model == m & res$ok]
  if (length(b) > 1) {
    cat(sprintf("%s: beta0 range across successful tries = %.3e\n",
                m, diff(range(b))))
  }
}

dir.create("results/real", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/real/stability_check.csv", row.names = FALSE)
