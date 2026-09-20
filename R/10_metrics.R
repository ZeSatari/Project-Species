## =====================================================================
## 10_metrics.R -- frequentist evaluation of Bayesian estimators
## =====================================================================
## These replace the quantity previously reported as "MSE", which was in
## fact the squared posterior standard deviation of a SINGLE fit
## (Table S2, comment: "based on posterior sd 0.059"). That is a measure
## of posterior spread, not of distance from the truth: it is small even
## when the estimator is badly biased.
##
## Everything here is computed ACROSS replicates. With R = 1 none of it
## is defined, which is the point.
## =====================================================================

## ---------------------------------------------------------------------
## CRPS of a posterior sample against a scalar truth.
##
##   CRPS(F, y) = E|X - y| - 0.5 * E|X - X'|,   X, X' ~ F iid
##
## For the empirical distribution of n draws, with x_(1) <= ... <= x_(n):
##   E|X - X'| = (2 / n^2) * sum_i (2i - n - 1) * x_(i)
## (verified for n = 2: x = (0,1) gives 1, matching the direct double sum)
## ---------------------------------------------------------------------
crps_from_sample <- function(x, truth) {
  x <- sort(as.numeric(x[is.finite(x)]))
  n <- length(x)
  if (n < 2L) return(NA_real_)
  term1 <- mean(abs(x - truth))
  term2 <- (2 / n^2) * sum((2 * seq_len(n) - n - 1) * x)
  term1 - 0.5 * term2
}

## Cross-check against scoringRules when it is available. Run once.
check_crps_implementation <- function(tol = 1e-8) {
  if (!requireNamespace("scoringRules", quietly = TRUE)) {
    message("scoringRules not installed; skipping CRPS cross-check.")
    return(invisible(NA))
  }
  set.seed(1)
  x <- rnorm(5000, 0.3, 1.2); y <- 0.5
  a <- crps_from_sample(x, y)
  b <- scoringRules::crps_sample(y = y, dat = x)
  stopifnot(abs(a - b) < tol)
  message("CRPS implementation matches scoringRules.")
  invisible(TRUE)
}

## ---------------------------------------------------------------------
## Per-replicate summary computed directly from an INLA marginal.
##
## Preferred over sampling the full latent field. inla.posterior.sample()
## reconstructs every spatial node on every draw -- profiling showed
## ~66 s per model for n = 2000, against ~86 s for the fit itself -- and
## none of that work is needed here: bias, RMSE, coverage and CRPS are
## all functions of the MARGINAL posterior of a single parameter.
##
## The interval comes from inla.qmarginal(), which is exact, so coverage
## carries no Monte Carlo error of its own. Only CRPS needs draws, and
## sampling a stored 1-D marginal is cheap.
##
## Joint draws are still required for derived quantities such as
## beta0 + gamma_i; see species_intercepts().
## ---------------------------------------------------------------------
summarise_marginal <- function(marg, truth, level = 0.95, n_crps = 4000) {
  a  <- (1 - level) / 2
  qs <- INLA::inla.qmarginal(c(a, 1 - a), marg)
  mu <- INLA::inla.emarginal(function(x) x, marg)
  s2 <- INLA::inla.emarginal(function(x) x^2, marg)
  draws <- INLA::inla.rmarginal(n_crps, marg)

  data.frame(
    truth    = truth,
    est      = mu,
    post_sd  = sqrt(max(s2 - mu^2, 0)),
    lwr      = qs[1],
    upr      = qs[2],
    covered  = as.integer(truth >= qs[1] && truth <= qs[2]),
    ci_width = qs[2] - qs[1],
    crps     = crps_from_sample(draws, truth)
  )
}

## ---------------------------------------------------------------------
## Per-replicate summary for one parameter.
## `draws` are posterior draws; `truth` is the data-generating value.
## ---------------------------------------------------------------------
summarise_param <- function(draws, truth, level = 0.95) {
  a  <- (1 - level) / 2
  qs <- unname(quantile(draws, c(a, 1 - a), na.rm = TRUE))
  data.frame(
    truth    = truth,
    est      = mean(draws, na.rm = TRUE),
    post_sd  = sd(draws, na.rm = TRUE),
    lwr      = qs[1],
    upr      = qs[2],
    covered  = as.integer(truth >= qs[1] && truth <= qs[2]),
    ci_width = qs[2] - qs[1],
    crps     = crps_from_sample(draws, truth)
  )
}

## ---------------------------------------------------------------------
## Aggregate across replicates.
## `per_rep` must have columns: rep, parameter, truth, est, covered,
## ci_width, crps.
##
## Reported quantities:
##   bias     = mean(est) - truth
##   rmse     = sqrt(mean((est - truth)^2))          <- the honest MSE^(1/2)
##   mc_se    = Monte Carlo standard error of the bias
##   coverage = empirical coverage of the nominal 95% CrI
## Coverage far below 0.95 is the diagnostic the manuscript currently
## lacks: it, not MSE, is what distinguishes "correctly widened
## uncertainty" from "biased estimate".
## ---------------------------------------------------------------------
aggregate_reps <- function(per_rep, level = 0.95) {
  per_rep |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      n_rep     = dplyr::n(),
      truth     = dplyr::first(truth),
      mean_est  = mean(est),
      bias      = mean(est) - dplyr::first(truth),
      mc_se     = sd(est) / sqrt(dplyr::n()),
      rmse      = sqrt(mean((est - truth)^2)),
      coverage  = mean(covered),
      cov_se    = sqrt(mean(covered) * (1 - mean(covered)) / dplyr::n()),
      mean_width = mean(ci_width),
      mean_crps = mean(crps, na.rm = TRUE),
      .groups   = "drop"
    )
}

## ---------------------------------------------------------------------
## Species-wise predictive comparison MSOM vs. SSOM.
##
## The previous pipeline compared model_msom$waic$waic / 120000 against
## model_single_i$waic$waic / 10000. Those are computed on DIFFERENT
## response sets; dividing by n does not make them comparable, and the
## resulting column mostly tracks prevalence (a near all-zero series has
## low deviance) rather than model quality.
##
## The comparable quantity is, for each species i, the log predictive
## density evaluated ONLY on that species' observations, under each of
## the two models:
##     delta_elpd_i = elpd_i(MSOM) - elpd_i(SSOM)
## together with its standard error over the pointwise differences.
## ---------------------------------------------------------------------
elpd_pointwise <- function(fit, idx = NULL) {
  cv <- INLA::inla.group.cv(result = fit, num.level.sets = 3)
  lp <- log(cv$cv)
  if (!is.null(idx)) lp <- lp[idx]
  lp[is.finite(lp)]
}

## ---------------------------------------------------------------------
## WITHDRAWN: predictive comparison between the community and the
## single-species models.
##
## The intention was to score both models on the same observations, as a
## principled replacement for the per-observation WAIC comparison (which
## is invalid across different response sets). It does not work here, and
## the reason is worth recording.
##
## inla.group.cv() builds its groups from the model's own dependence
## structure. Under the community model, which carries a latent field
## shared by all species, those groups extend ACROSS species: predicting
## one observation of species i also withholds observations of other
## species. No comparable withholding happens in a model fitted to a
## single species. The two scores therefore answer different questions,
## and the comparison is systematically unfavourable to the community
## model.
##
## The empirical run made this plain: the difference was negative for all
## twelve species, between -55 and -730 nats over 3,357 observations, and
## its magnitude tracked prevalence rather than anything about model
## quality -- the signature of an artefact.
##
## What replaces it: the shrinkage comparison (species-level intercepts
## under each model, reported alongside naive occupancy) and the
## replicated MSOM simulation, where the truth is known.
##
## A valid comparison would need cross-validation with the folds defined
## by the analyst -- withholding whole sites and predicting them under
## both models -- rather than automatic group construction. That is a
## larger exercise and is noted as future work.
## ---------------------------------------------------------------------

## ---------------------------------------------------------------------
## Naive occupancy: proportion of site-years with at least one detection.
## Replaces the "Detection%" column, which was previously computed from
## the FIRST visit only, with missing visits left in the denominator, and
## was never defined in the manuscript.
## ---------------------------------------------------------------------
naive_occupancy <- function(Y) {
  stopifnot(is.matrix(Y))
  ok <- rowSums(!is.na(Y)) > 0
  mean(rowSums(Y[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0)
}

detection_rate <- function(Y) {
  stopifnot(is.matrix(Y))
  mean(Y == 1, na.rm = TRUE)   # over all non-missing site-visit records
}
