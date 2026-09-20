## =====================================================================
## 99_diag_group.R -- which A.local / group convention is correct?
## =====================================================================
## Self-contained. Does not source anything. Runs in seconds: tiny mesh,
## few observations, Gaussian likelihood (the question is about the
## f() term, not the likelihood).
##
## Ground truth is a known AR(1) field with rho = 0.8. A convention is
## correct if it (a) runs, (b) emits no "unused groups" warning, and
## (c) recovers GroupRho near 0.8. A variant that runs cleanly but
## returns GroupRho near 0 has silently disabled the temporal structure
## -- which is the failure mode we are hunting.
## =====================================================================

suppressPackageStartupMessages({ library(INLA); library(fmesher) })
INLA::inla.setOption(num.threads = "1:1")
set.seed(42)

## --- small domain -----------------------------------------------------
nT    <- 4
n_obs <- 600
loc   <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 10))
time  <- sample(seq_len(nT), n_obs, replace = TRUE)

mesh <- fmesher::fm_mesh_2d(loc.domain = cbind(c(0,10,10,0), c(0,0,10,10)),
                            max.edge = c(1.5, 3), offset = c(0.5, 2))
spde <- inla.spde2.pcmatern(mesh, prior.range = c(3, 0.5),
                            prior.sigma = c(1, 0.5))
cat("mesh nodes:", mesh$n, " spde$n.spde:", spde$n.spde, "\n")

## --- simulate a genuine AR(1)-in-time field ---------------------------
RHO_TRUE <- 0.8
Q  <- inla.spde.precision(spde, theta = c(log(3), log(1)))
e  <- as.matrix(inla.qsample(n = nT, Q = Q, seed = 7))
u  <- e
for (t in 2:nT) u[, t] <- RHO_TRUE * u[, t-1] + sqrt(1-RHO_TRUE^2) * e[, t]

A_grouped <- inla.spde.make.A(mesh, loc = loc, group = time, n.group = nT)
A_spatial <- inla.spde.make.A(mesh, loc = loc)
cat("dim(A_grouped):", dim(A_grouped),
    "  expected cols:", spde$n.spde * nT, "\n")
cat("dim(A_spatial):", dim(A_spatial), "\n")

eta <- as.numeric(A_grouped %*% as.vector(u))
y   <- 1.5 + eta + rnorm(n_obs, 0, 0.3)

idx_stack <- inla.spde.make.index("sf", n.spde = spde$n.spde, n.group = nT)

## --- the four candidate conventions -----------------------------------
variants <- list(
  "A: A_grouped + group=time (obs-level)" = list(
    data = list(y = y, Int = rep(1, n_obs), sf = rep(NA, n_obs), grp = time),
    form = y ~ -1 + Int + f(sf, model = spde, A.local = A_grouped,
                            group = grp,
                            control.group = list(model = "ar1"))
  ),
  "B: A_grouped + ngroup=nT, no group=" = list(
    data = list(y = y, Int = rep(1, n_obs), sf = rep(NA, n_obs)),
    form = y ~ -1 + Int + f(sf, model = spde, A.local = A_grouped,
                            ngroup = nT,
                            control.group = list(model = "ar1"))
  ),
  "C: A_grouped + make.index (field-level)" = list(
    data = c(list(y = y, Int = rep(1, n_obs)), idx_stack),
    form = y ~ -1 + Int + f(sf, model = spde, A.local = A_grouped,
                            group = sf.group,
                            control.group = list(model = "ar1"))
  ),
  "D: A_spatial + group=time (obs-level)" = list(
    data = list(y = y, Int = rep(1, n_obs), sf = rep(NA, n_obs), grp = time),
    form = y ~ -1 + Int + f(sf, model = spde, A.local = A_spatial,
                            group = grp,
                            control.group = list(model = "ar1"))
  )
)

## --- run each, capturing warnings -------------------------------------
run_variant <- function(v) {
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      inla(v$form, data = v$data, family = "gaussian",
           control.compute = list(waic = TRUE),
           control.inla = list(int.strategy = "eb"),
           verbose = FALSE),
      error = function(e) structure(conditionMessage(e), class = "failed")
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warns = warns)
}

results <- lapply(variants, run_variant)

## --- report -----------------------------------------------------------
cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("TRUE GroupRho = %.2f\n\n", RHO_TRUE))

for (nm in names(results)) {
  r <- results[[nm]]
  cat(nm, "\n")
  if (inherits(r$fit, "failed")) {
    cat("   FAILED: ", as.character(r$fit), "\n\n", sep = "")
    next
  }
  hp  <- r$fit$summary.hyperpar
  rho <- if ("GroupRho for sf" %in% rownames(hp))
           hp["GroupRho for sf", "mean"] else NA_real_
  unused <- any(grepl("unused groups", r$warns))
  cat(sprintf("   GroupRho = %-7s  WAIC = %-9.1f  unused-groups warning: %s\n",
              formatC(rho, format = "f", digits = 3),
              r$fit$waic$waic, if (unused) "YES" else "no"))
  if (length(r$warns)) {
    for (w in r$warns) cat("     warn: ", substr(w, 1, 90), "\n", sep = "")
  }
  cat("\n")
}

cat(strrep("=", 74), "\n")
cat("Adopt the variant that runs, emits no unused-groups warning, and\n")
cat("recovers GroupRho near", RHO_TRUE, "with the lowest WAIC.\n")
cat("A variant that runs cleanly but returns GroupRho ~ 0 has silently\n")
cat("switched the temporal structure off -- reject it.\n")
