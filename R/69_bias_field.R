## =====================================================================
## 69_bias_field.R -- does the latent field produce the slope bias?
## =====================================================================
## 68_bias_source.R stripped the model to one species, a linear
## covariate, no field and no random effects, and found no bias: across
## sample size, detection probability, visits and approximation strategy
## the slope came back within one or two Monte Carlo standard errors of
## 0.700, with coverage near nominal. The occupancy likelihood is not the
## source.
##
## Four things present in the simulations were absent from that test:
## the spatio-temporal field, a quadratic covariate, species random
## effects, and detection covariates. This script adds the first two
## back, separately and together, holding everything else fixed. Species
## random effects are left for a later test, since they require the
## multi-species design.
##
## The competing explanations:
##
##   (a) The covariate and the field compete. Both are smooth spatial
##       surfaces; where the covariate is itself spatially structured,
##       the field can absorb part of its effect. This is the spatial
##       confounding already documented in the empirical elevation
##       coefficients, and it would bias the slope downward.
##
##   (b) The covariate's own geometry. A quadratic covariate is
##       chi-squared, and its right tail drives psi toward one where the
##       data carry little slope information. 68_ used a linear
##       covariate; the single-species simulation used x^2.
##
## Design: 2 covariate forms x 2 field settings x 2 covariate spatial
## structures, R replicates each. Field present or absent is the axis
## that matters; the rest separate (a) from (b).
##
## Runtime: the cells with a field are much slower. Roughly 1-2 hours.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

OUTDIR <- "results/bias"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

R_REP <- 100
SEED0 <- 20260601L
N_SITE <- 500
N_TIME <- 5
K_VIS  <- 3

TRUTH <- list(beta0 = -0.5, beta1 = 0.7, alpha0 = qlogis(0.30))

## Domain and mesh, matching 20_simulate_single.R.
DOM <- local({
  g <- expand.grid(x = seq(1.5, 298.5, by = 3), y = seq(1.5, 298.5, by = 3))
  bnd <- fmesher::fm_nonconvex_hull(as.matrix(g), convex = 30)
  mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(25, 50),
                              offset = c(-0.1, -0.2))
  list(grid = g, mesh = mesh)
})
message(sprintf("mesh: %d nodes, %d grid cells", DOM$mesh$n, nrow(DOM$grid)))

det_X_intercept_only <- function(n, K) {
  X <- matrix(1, nrow = n, ncol = K)
  attr(X, "n_det_coef") <- 1L
  attr(X, "na_action")  <- "zero"
  X
}

## ---------------------------------------------------------------------
## Generator. `cov_spatial` decides whether the covariate is a smooth
## spatial surface (as in the paper's simulations) or independent noise;
## `cov_form` whether it enters as x or x^2.
## ---------------------------------------------------------------------
sim_one <- function(seed, cov_spatial, cov_form, field_in_truth = TRUE) {
  set.seed(seed)
  idx <- sample(nrow(DOM$grid), N_SITE)
  loc <- as.matrix(DOM$grid[idx, ])

  ## A Matern surface over the sites, drawn directly from the SPDE
  ## precision. Generating parameters (range 100, sigma 0.5) match
  ## 20_simulate_single.R; the FITTED prior differs deliberately, so that
  ## recovery is not flattered by a prior centred on the truth.
  Q  <- INLA::inla.spde2.precision(
          INLA::inla.spde2.pcmatern(DOM$mesh, prior.range = c(100, 0.5),
                                    prior.sigma = c(0.5, 0.5)),
          theta = c(log(100), log(0.5)))
  w  <- as.numeric(INLA::inla.qsample(1, Q, seed = seed))
  A  <- INLA::inla.spde.make.A(DOM$mesh, loc = loc)
  fld_site <- as.numeric(A %*% w)
  rho <- 0.65
  fld <- numeric(N_SITE * N_TIME)
  prev <- fld_site
  for (tt in seq_len(N_TIME)) {
    prev <- if (tt == 1) fld_site
            else rho * prev + sqrt(1 - rho^2) * as.numeric(A %*%
                   as.numeric(INLA::inla.qsample(1, Q, seed = seed + 1000 * tt)))
    fld[((tt - 1) * N_SITE + 1):(tt * N_SITE)] <- prev
  }

  ## covariate: spatially structured or not
  if (cov_spatial) {
    w2 <- as.numeric(INLA::inla.qsample(1, Q, seed = seed + 99991))
    xs <- as.numeric(A %*% w2)
  } else {
    xs <- rnorm(N_SITE)
  }
  xs <- as.numeric(scale(xs))
  fx <- if (cov_form == "quadratic") xs^2 else xs
  fx_rep <- rep(fx, times = N_TIME)

  eta <- TRUTH$beta0 + TRUTH$beta1 * fx_rep + if (field_in_truth) fld else 0
  z   <- rbinom(length(eta), 1, plogis(eta))
  Y   <- matrix(rbinom(length(z) * K_VIS, 1, z * 0.30), length(z), K_VIS)

  list(Y = Y, fx = fx_rep, loc = loc,
       time = rep(seq_len(N_TIME), each = N_SITE))
}

fit_one <- function(d, fit_field) {
  n <- nrow(d$Y)
  Xd <- det_X_intercept_only(n, K_VIS)
  base <- list(Y = d$Y, X = Xd, Int_occ = rep(1, n), xx = d$fx)

  if (fit_field) {
    spde <- INLA::inla.spde2.pcmatern(DOM$mesh, prior.range = c(150, 0.5),
                                      prior.sigma = c(1.0, 0.5))
    stf <- make_st_field(DOM$mesh, spde,
                         as.data.frame(d$loc[rep(seq_len(N_SITE), N_TIME), ]),
                         time = d$time, n_time = N_TIME)
    f <- inla.mdata(Y, X) ~ -1 + Int_occ + xx +
      f(spatialfield, model = spde, A.local = stf$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))
    dat <- c(base, stf$idx)
  } else {
    f <- inla.mdata(Y, X) ~ -1 + Int_occ + xx
    dat <- base
  }

  fit <- tryCatch(fit_occupancy(f, dat, Xd, compute_cv = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(beta1 = NA_real_, covered = NA_real_))
  s <- fit$summary.fixed["xx", ]
  c(beta1 = s[["mean"]],
    covered = as.numeric(s[["0.025quant"]] <= TRUTH$beta1 &
                         TRUTH$beta1 <= s[["0.975quant"]]))
}

run_cell <- function(cov_spatial, cov_form, fit_field, R = R_REP) {
  lab <- sprintf("cov %s / %s | field %s",
                 if (cov_spatial) "spatial" else "iid   ",
                 substr(cov_form, 1, 6),
                 if (fit_field) "fitted " else "omitted")
  message("  ", lab)
  est <- vapply(seq_len(R), function(r) {
    d <- sim_one(SEED0 + r, cov_spatial, cov_form, field_in_truth = TRUE)
    fit_one(d, fit_field)
  }, numeric(2))
  b <- est["beta1", ]; ok <- !is.na(b)
  data.frame(cov_spatial = cov_spatial, cov_form = cov_form,
             field_fitted = fit_field, n_ok = sum(ok),
             bias = mean(b[ok]) - TRUTH$beta1,
             mc_se = sd(b[ok]) / sqrt(sum(ok)),
             coverage = mean(est["covered", ok]))
}

## One cell at R = 5 before committing to 800 fits, since the cells with
## a fitted field are the slow ones and a broken design would waste hours.
{
  message("smoke test: one cell at R = 5")
  s <- run_cell(TRUE, "linear", TRUE, R = 5)
  message(sprintf("  n_ok %d, bias %+.3f -- proceeding", s$n_ok, s$bias))
  if (s$n_ok == 0) stop("smoke test: every fit failed; check the field design.")
}

cells <- list()
for (cs in c(FALSE, TRUE))
  for (cf in c("linear", "quadratic"))
    for (ff in c(FALSE, TRUE))
      cells[[length(cells) + 1]] <- run_cell(cs, cf, ff)

out <- do.call(rbind, cells); rownames(out) <- NULL
cat("\n", strrep("=", 88), "\n", sep = "")
print(out, row.names = FALSE, digits = 4)
write.csv(out, file.path(OUTDIR, "bias_field.csv"), row.names = FALSE)

cat("\n", strrep("-", 88), "\n", sep = "")
cat("The field is generated in every cell; only whether it is FITTED varies.\n")
cat("Omitting a field that is present is a misspecification, so those cells\n")
cat("are not a clean control -- they show what the covariate alone recovers.\n\n")

for (cs in c(FALSE, TRUE)) for (cf in c("linear", "quadratic")) {
  d <- out[out$cov_spatial == cs & out$cov_form == cf, ]
  if (nrow(d) == 2) {
    cat(sprintf("cov %-8s %-9s : field omitted %+.4f | field fitted %+.4f\n",
                if (cs) "spatial" else "iid", cf,
                d$bias[!d$field_fitted], d$bias[d$field_fitted]))
  }
}

cat("\nFor reference: -0.024 in the single-species simulation (spatial\n")
cat("covariate, quadratic, field fitted) and -0.050 in the multi-species one\n")
cat("(spatial covariate, linear, field fitted, plus species effects).\n")
cat("\nIf the bias appears only where the covariate is spatially structured\n")
cat("AND the field is fitted, competition between them is the explanation.\n")
