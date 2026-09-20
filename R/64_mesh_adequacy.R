## =====================================================================
## 64_mesh_adequacy.R -- is the coarsened mesh adequate, or merely stable?
## =====================================================================
## The manuscript coarsened the mesh from 1724 nodes to 283 because the
## fine mesh gave irreproducible fits. That establishes stability. It does
## not establish that the coarse mesh resolves the spatial structure: a
## mesh can be perfectly stable and still be too crude to represent the
## field, in which case the range would be inflated, sigma attenuated,
## and the predicted surface oversmoothed.
##
## Stability and adequacy are different properties and are checked
## differently. Stability is repeated fits of one mesh agreeing with each
## other; adequacy is different meshes agreeing with one another.
##
## Three meshes, one model (M1 for REVI), everything else held fixed:
##   coarse       283 nodes   the specification used in the paper
##   intermediate ~600 nodes
##   fine        1724 nodes   the original, known to be unstable
##
## The fine mesh is fitted N_REPEAT times, since a single fit from it is
## not reproducible; the spread across those fits sets the scale against
## which any difference between meshes should be judged. A difference
## smaller than the fine mesh's own run-to-run variation is not evidence
## of inadequacy.
##
## Runtime: dominated by the fine mesh. Allow several hours.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/mesh"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR)

N_REPEAT_FINE <- 3      # the fine mesh is refitted; the others once

N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
K      <- dim(hbefTrends$y)[4]

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
Y_revi <- stack_species("REVI")
X_M1   <- build_det_X(list(stack_det_cov("day"), stack_det_cov("tod")), K = K)

bnd <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)

MESHES <- list(
  coarse       = list(max.edge = c(0.35, 1.2), cutoff = 0.35,
                      offset = c(0.05, 1.5), min.angle = 21),
  intermediate = list(max.edge = c(0.20, 0.9), cutoff = 0.20,
                      offset = c(0.03, 1.2), min.angle = 21),
  fine         = list(max.edge = c(0.10, 0.7), cutoff = 0.12,
                      offset = c(0.01, 1.0), min.angle = 20)
)

build_mesh <- function(cfg) {
  fmesher::fm_mesh_2d(boundary = bnd, max.edge = cfg$max.edge,
                      cutoff = cfg$cutoff, offset = cfg$offset,
                      min.angle = cfg$min.angle)
}
for (nm in names(MESHES)) {
  m <- build_mesh(MESHES[[nm]])
  message(sprintf("%-13s %4d nodes | field %5d | ratio %.2f",
                  nm, m$n, m$n * N_YEAR, m$n * N_YEAR / nrow(Xocc)))
}

## ---------------------------------------------------------------------
## Fit one mesh. Priors, data and formula are identical throughout; only
## the triangulation changes.
## ---------------------------------------------------------------------
fit_mesh <- function(cfg, tag) {
  mesh <- build_mesh(cfg)
  spde <- INLA::inla.spde2.pcmatern(mesh,
                                    prior.range = PRIORS$spde_real$range,
                                    prior.sigma = PRIORS$spde_real$sigma,
                                    constr = TRUE)
  stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                       time = Xocc$time, n_time = N_YEAR)
  dat <- c(as.list(Xocc), stf$idx, list(Y = Y_revi, X = X_M1))
  f <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))

  fit <- tryCatch(fit_occupancy(f, dat, X_M1, compute_cv = FALSE),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    message("  ", tag, ": FAILED -- ", conditionMessage(fit))
    return(NULL)
  }

  hp <- fit$summary.hyperpar
  gr <- function(p) { i <- grep(p, rownames(hp)); if (length(i) == 1) hp[i, "mean"] else NA_real_ }

  ## Predicted occupancy at each site in the final year, so that meshes
  ## are compared on the surface they produce and not only on parameters.
  A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))
  s <- INLA::inla.posterior.sample(n = 300, result = fit)
  nm_l <- rownames(s[[1]]$latent)
  i_int <- grep("^Int_occ:", nm_l); i_e1 <- grep("^elev_s:", nm_l)
  i_e2 <- grep("^elev_s2:", nm_l); i_f <- grep("^spatialfield:", nm_l)
  n_spde <- length(i_f) / N_YEAR
  sl <- i_f[((N_YEAR - 1) * n_spde + 1):(N_YEAR * n_spde)]
  psi <- rowMeans(vapply(s, function(d) {
    v <- d$latent
    plogis(v[i_int] + v[i_e1] * elev_s + v[i_e2] * elev_s^2 +
             as.numeric(A_site %*% v[sl]))
  }, numeric(N_SITE)))

  list(
    summary = data.frame(
      mesh = tag, nodes = mesh$n,
      beta0   = fit$summary.fixed["Int_occ", "mean"],
      elev    = fit$summary.fixed["elev_s", "mean"],
      elev2   = fit$summary.fixed["elev_s2", "mean"],
      range   = gr("^Range for"), sigma = gr("^Stdev for"),
      rho     = gr("^GroupRho for"),
      waic    = fit$waic$waic, secs = fit$elapsed
    ),
    psi = psi
  )
}

## ---------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------
res <- list(); psis <- list()

for (nm in c("coarse", "intermediate")) {
  message("fitting ", nm)
  r <- fit_mesh(MESHES[[nm]], nm)
  if (!is.null(r)) { res[[nm]] <- r$summary; psis[[nm]] <- r$psi }
}

## The fine mesh is refitted: one fit from it is not reproducible, and
## its own spread is the yardstick for the between-mesh comparison.
for (k in seq_len(N_REPEAT_FINE)) {
  tag <- paste0("fine_", k)
  message("fitting ", tag)
  r <- fit_mesh(MESHES$fine, tag)
  if (!is.null(r)) { res[[tag]] <- r$summary; psis[[tag]] <- r$psi }
}

out <- do.call(rbind, res); rownames(out) <- NULL
cat("\n", strrep("=", 84), "\n", sep = "")
print(out, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------
## Between-mesh differences against the fine mesh's own variability.
## ---------------------------------------------------------------------
fine <- out[grepl("^fine", out$mesh), ]
if (nrow(fine) >= 2) {
  cat("\nfine-mesh run-to-run spread (max - min):\n")
  for (v in c("beta0", "range", "sigma", "rho")) {
    cat(sprintf("  %-6s %.4f\n", v, diff(range(fine[[v]], na.rm = TRUE))))
  }
}
if (nrow(fine) >= 1 && "coarse" %in% out$mesh) {
  co <- out[out$mesh == "coarse", ]
  cat("\ncoarse minus fine (fine averaged over its repeats):\n")
  for (v in c("beta0", "range", "sigma", "rho")) {
    d <- co[[v]] - mean(fine[[v]], na.rm = TRUE)
    cat(sprintf("  %-6s %+.4f\n", v, d))
  }
}

## Predicted surfaces. If the coarse mesh is oversmoothing, its
## predictions will differ from the fine mesh's by more than the fine
## mesh differs from itself.
if (length(psis) >= 2) {
  cat("\npredicted occupancy across the 373 sites, pairwise:\n")
  nms <- names(psis)
  for (i in seq_along(nms)) for (j in seq_along(nms)) if (i < j) {
    a <- psis[[nms[i]]]; b <- psis[[nms[j]]]
    cat(sprintf("  %-14s vs %-14s  cor %.4f  max abs diff %.4f  mean abs diff %.4f\n",
                nms[i], nms[j], cor(a, b), max(abs(a - b)), mean(abs(a - b))))
  }
}

cat("\nReading this. The coarse mesh is adequate, not merely stable, if its\n")
cat("hyperparameters and predicted surface sit within the fine mesh's own\n")
cat("run-to-run variation. A systematically inflated range or attenuated sigma\n")
cat("at the coarse resolution would indicate oversmoothing, and the manuscript\n")
cat("would then need to say that the mesh was chosen for stability at some cost\n")
cat("in spatial resolution.\n")

write.csv(out, file.path(OUTDIR, "mesh_comparison.csv"), row.names = FALSE)
saveRDS(psis, file.path(OUTDIR, "mesh_predictions.rds"))
