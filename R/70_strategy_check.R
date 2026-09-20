## =====================================================================
## 70_strategy_check.R -- how good is the Laplace approximation here?
## =====================================================================
## Belmont et al. (2024) state that the occupancy log-likelihood is NOT
## log-concave. The accuracy of INLA's approximation therefore cannot be
## assumed for this model class -- it has to be checked.
##
## This refits ONE model under three latent-field strategies. If the
## posteriors agree, that is a positive check worth reporting. If they
## disagree, an MCMC benchmark is not optional.
##
## Cost: three fits of a single model. Run it early.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")
source("R/30_fit_occupancy.R")

OUTDIR <- "results/strategy"
record_session(OUTDIR)

## ---------------------------------------------------------------------
## Build ONE dataset. Use "independent": the regime in which the method
## is supposed to work, so any disagreement below is attributable to the
## approximation rather than to model misspecification.
## ---------------------------------------------------------------------
DOM <- make_domain()
sim <- simulate_single(DOM, scenario = "independent", seed = 20260101L)
d   <- sim$data
K   <- sim$truth$K

spde <- INLA::inla.spde2.pcmatern(DOM$mesh,
                                  prior.range = PRIORS$spde_sim$range,
                                  prior.sigma = PRIORS$spde_sim$sigma)
stf <- make_st_field(DOM$mesh, spde, d[, c("x.loc", "y.loc")],
                     time = d$time, n_time = sim$truth$nT)

Y  <- sim$Y_A
gA <- as.matrix(d[, paste0("gA", seq_len(K))])
yB <- as.matrix(d[, paste0("yB", seq_len(K))])
X2 <- build_det_X(list(gA, yB), K = K, na_action = "zero")

## A list, not a data.frame: the field index is longer than the
## observation vectors.
dat <- c(list(Y = Y, X = X2,
              Int_occ = rep(1, nrow(d)),
              x_sq    = d$x_A^2),
         stf$idx)

f <- inla.mdata(Y, X) ~ -1 + Int_occ + x_sq +
  f(spatialfield, model = spde, A.local = stf$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

## ---------------------------------------------------------------------
## Refit under each strategy. Everything else -- data, formula, priors,
## int.strategy -- is held fixed, so the only thing varying is the
## approximation to the latent marginals.
##
## Routed through fit_occupancy(): calling inla() directly here would
## bypass the bookkeeping that det_draws() depends on, which is exactly
## the kind of divergence this pipeline is meant to prevent.
##
## compute_cv = FALSE: inla.group.cv() is expensive and irrelevant to the
## question being asked here.
## ---------------------------------------------------------------------
STRATEGIES <- c("gaussian", "simplified.laplace", "laplace")

fits <- lapply(STRATEGIES, function(st) {
  message("strategy = ", st)
  fit_occupancy(f, dat, X2, strategy = st, compute_cv = FALSE)
})
names(fits) <- STRATEGIES

## ---------------------------------------------------------------------
## Occupancy fixed effects.
## ---------------------------------------------------------------------
fixed_tab <- do.call(rbind, lapply(STRATEGIES, function(st) {
  s <- fits[[st]]$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  data.frame(strategy = st, parameter = rownames(s), s, row.names = NULL)
}))

## Detection coefficients (hyperparameters of the occupancy family).
det_tab <- do.call(rbind, lapply(STRATEGIES, function(st) {
  dr <- det_draws(fits[[st]])
  data.frame(
    strategy  = st,
    parameter = colnames(dr),
    mean      = colMeans(dr),
    sd        = apply(dr, 2, sd),
    row.names = NULL
  )
}))

timing <- data.frame(strategy = STRATEGIES,
                     seconds  = vapply(fits, `[[`, numeric(1), "elapsed"),
                     waic     = vapply(fits, \(f) f$waic$waic, numeric(1)),
                     row.names = NULL)

## ---------------------------------------------------------------------
## The decision rule.
##
## Compare "laplace" (most accurate, slowest) against
## "simplified.laplace" (the default). Express the gap in units of the
## posterior SD -- a shift of 0.02 means nothing if sd = 0.5, and a great
## deal if sd = 0.01.
##
##   < 0.10 sd  approximation adequate; report as a positive check
##   0.10-0.25  borderline; run the MCMC benchmark
##   > 0.25 sd  the default approximation is not trustworthy here;
##              MCMC validation is required before any claim is made
## ---------------------------------------------------------------------
cmp <- merge(
  subset(fixed_tab, strategy == "simplified.laplace",
         select = c(parameter, mean, sd)),
  subset(fixed_tab, strategy == "laplace",
         select = c(parameter, mean, sd)),
  by = "parameter", suffixes = c("_sl", "_l")
)
cmp$shift_in_sd <- abs(cmp$mean_sl - cmp$mean_l) / cmp$sd_l
cmp$sd_ratio    <- cmp$sd_sl / cmp$sd_l

cmp_det <- merge(
  subset(det_tab, strategy == "simplified.laplace",
         select = c(parameter, mean, sd)),
  subset(det_tab, strategy == "laplace",
         select = c(parameter, mean, sd)),
  by = "parameter", suffixes = c("_sl", "_l")
)
cmp_det$shift_in_sd <- abs(cmp_det$mean_sl - cmp_det$mean_l) / cmp_det$sd_l
cmp_det$sd_ratio    <- cmp_det$sd_sl / cmp_det$sd_l

all_cmp <- rbind(cbind(block = "occupancy", cmp),
                 cbind(block = "detection", cmp_det))

write.csv(fixed_tab, file.path(OUTDIR, "fixed_by_strategy.csv"), row.names = FALSE)
write.csv(det_tab,   file.path(OUTDIR, "detection_by_strategy.csv"), row.names = FALSE)
write.csv(timing,    file.path(OUTDIR, "timing.csv"), row.names = FALSE)
write.csv(all_cmp,   file.path(OUTDIR, "strategy_gap.csv"), row.names = FALSE)

print(timing, digits = 4)
cat("\n-- simplified.laplace vs laplace, gap in posterior SD units --\n")
print(all_cmp[, c("block", "parameter", "mean_sl", "mean_l",
                  "shift_in_sd", "sd_ratio")], digits = 3)

worst <- max(all_cmp$shift_in_sd, na.rm = TRUE)
cat(sprintf("\nLargest shift: %.3f posterior SD -> %s\n", worst,
            if (worst < 0.10) "approximation adequate; report as a check."
            else if (worst < 0.25) "borderline; run the MCMC benchmark."
            else "NOT trustworthy; MCMC validation required."))

## ---------------------------------------------------------------------
## Optional: overlay the marginals. Visual disagreement in the tails can
## matter even when the means agree, and the tails are what the credible
## intervals depend on.
## ---------------------------------------------------------------------
plot_marginal <- function(param = "Int_occ") {
  ms <- lapply(fits, \(f) f$marginals.fixed[[param]])
  xr <- range(sapply(ms, \(m) range(m[, 1])))
  yr <- range(sapply(ms, \(m) range(m[, 2])))
  plot(NA, xlim = xr, ylim = yr, xlab = param, ylab = "density",
       main = paste("Posterior marginal:", param))
  for (i in seq_along(ms)) lines(ms[[i]], col = i, lwd = 2)
  legend("topright", legend = STRATEGIES, col = seq_along(ms), lwd = 2, bty = "n")
}

pdf(file.path(OUTDIR, "marginals.pdf"), width = 7, height = 5)
for (p in rownames(fits[[1]]$summary.fixed)) plot_marginal(p)
dev.off()
