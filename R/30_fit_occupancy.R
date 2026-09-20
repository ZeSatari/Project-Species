## =====================================================================
## 30_fit_occupancy.R -- ONE fitting entry point for every model
## =====================================================================
## The previous pipeline fitted models in five scripts with silently
## divergent settings: control.fixed prec = 1 vs 1/2.72, detection priors
## specified in one script and defaulted in another, constr = TRUE on the
## SPDE for the real data but not the simulation, and int.strategy = "eb"
## throughout the simulations. Routing everything through one function
## makes that class of divergence impossible.
## =====================================================================

## ---------------------------------------------------------------------
## Detection design matrix for family = "occupancy".
## inla.mdata(Y, X) expects X laid out visit-major:
##   [1, cov1_v1, cov2_v1, | 1, cov1_v2, cov2_v2, | ...]
##
## `na_action` makes the missing-data convention EXPLICIT and testable.
## The manuscript states that missing observations of the co-occurring
## species are treated as non-detections, but the code passed raw NAs.
## Whichever convention is adopted must be stated and its sensitivity
## reported.
## ---------------------------------------------------------------------
build_det_X <- function(cov_list, K, na_action = c("zero", "keep")) {
  na_action <- match.arg(na_action)
  stopifnot(is.list(cov_list), length(cov_list) >= 1)

  n_coef <- length(cov_list) + 1L      # intercept + covariates
  if (exists("MAX_DET_COEF") && n_coef > MAX_DET_COEF) {
    stop("build_det_X(): ", n_coef, " detection coefficients requested but ",
         "the occupancy family provides only ", MAX_DET_COEF, " slots.")
  }
  if (exists("MAX_DET_COEF") && n_coef > MAX_DET_COEF - 2L) {
    warning("build_det_X(): ", n_coef, " of ", MAX_DET_COEF,
            " detection slots in use. The Laplace approximation degrades ",
            "as the hyperparameter dimension grows; validate against MCMC.")
  }

  dims <- vapply(cov_list, \(m) c(nrow(m), ncol(m)), numeric(2))
  if (length(unique(dims[1, ])) != 1L) stop("build_det_X(): differing n rows.")
  if (any(dims[2, ] != K))             stop("build_det_X(): each covariate needs K columns.")

  cov_list <- lapply(cov_list, \(m) {
    if (na_action == "zero") m[is.na(m)] <- 0
    m
  })

  blocks <- lapply(seq_len(K), \(j)
    cbind(1, do.call(cbind, lapply(cov_list, \(m) m[, j])))
  )
  X <- do.call(cbind, blocks)
  attr(X, "n_det_coef") <- length(cov_list) + 1L   # intercept + covariates
  attr(X, "na_action")  <- na_action
  X
}

## ---------------------------------------------------------------------
## Detection-coefficient hyperparameter priors.
## A prior must be supplied for EVERY detection coefficient. The previous
## code set priors for beta1 and beta2 only, so in model M2 -- which has
## three detection coefficients -- the coefficient on y_B silently took
## the INLA default. M1 and M2 therefore had non-comparable priors, which
## contaminates the WAIC comparison between them.
## ---------------------------------------------------------------------
det_hyper <- function(n_coef, prior = PRIORS$det_beta) {
  stopifnot(n_coef >= 1)
  h <- lapply(seq_len(n_coef), \(i)
    list(param = c(prior$mean, prior$prec), initial = 0))
  names(h) <- paste0("beta", seq_len(n_coef))
  h
}

## ---------------------------------------------------------------------
## Fit.
## ---------------------------------------------------------------------
fit_occupancy <- function(formula, data_list, det_X,
                          det_prior   = PRIORS$det_beta,
                          fixed_prior = PRIORS$fixed,
                          int_strategy = OPTS$int_strategy,
                          strategy     = "simplified.laplace",
                          compute_cv  = TRUE) {

  n_coef <- attr(det_X, "n_det_coef")
  if (is.null(n_coef)) stop("det_X must come from build_det_X().")

  t0 <- proc.time()[["elapsed"]]

  fit <- INLA::inla(
    formula, data = data_list, family = "occupancy",
    control.fixed   = fixed_prior,
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE,
                           cpo = FALSE),
    control.inla    = list(int.strategy = int_strategy,
                           strategy     = strategy),
    control.family  = list(
      control.link = list(model = "logit"),
      link.simple  = "logit",
      hyper        = det_hyper(n_coef, det_prior)
    ),
    verbose = OPTS$verbose
  )

  fit$elapsed <- proc.time()[["elapsed"]] - t0

  if (compute_cv) {
    cv <- INLA::inla.group.cv(result = fit, num.level.sets = 3)
    fit$ulogcv <- mean(log(cv$cv), na.rm = TRUE)
    fit$cv_pointwise <- log(cv$cv)
  }
  ## Set LAST, and always: det_draws() depends on it. Any code path that
  ## calls inla() directly instead of going through this function will
  ## leave it NULL and fail downstream.
  fit$n_det_coef <- n_coef
  fit
}

## ---------------------------------------------------------------------
## Detection-coefficient MARGINALS, keyed alpha0, alpha1, ...
##
## Cheap: these are already stored on the fitted object. Use this rather
## than det_draws() wherever only per-parameter summaries are needed.
## ---------------------------------------------------------------------
det_marginals <- function(fit, pattern = NULL) {
  if (is.null(pattern)) pattern <- "^beta\\[[0-9]+\\]|^beta[0-9]+"
  nm  <- rownames(fit$summary.hyperpar)
  idx <- grep(pattern, nm)
  if (is.null(fit$n_det_coef)) {
    stop("det_marginals(): fit$n_det_coef is NULL; fit via fit_occupancy().")
  }
  if (length(idx) != fit$n_det_coef) {
    stop("det_marginals(): expected ", fit$n_det_coef,
         " detection coefficients, matched ", length(idx),
         ".\n  rownames = ", paste(nm, collapse = " | "))
  }
  out <- fit$marginals.hyperpar[nm[idx]]
  names(out) <- paste0("alpha", seq_along(out) - 1L)
  out
}

## Occupancy fixed-effect marginals, by name.
fixed_marginals <- function(fit, names_wanted) {
  miss <- setdiff(names_wanted, names(fit$marginals.fixed))
  if (length(miss)) {
    stop("fixed_marginals(): not found: ", paste(miss, collapse = ", "),
         "\n  available: ", paste(names(fit$marginals.fixed), collapse = ", "))
  }
  fit$marginals.fixed[names_wanted]
}

## ---------------------------------------------------------------------
## Posterior draws for the DETECTION coefficients.
##
## These are hyperparameters of the occupancy family, not fixed effects,
## so they live in $summary.hyperpar / inla.hyperpar.sample(), NOT in
## $summary.fixed.
##
## !! I have not been able to execute R-INLA to confirm the exact row
## !! naming for this family in your build. The function therefore FAILS
## !! LOUDLY rather than silently returning the wrong columns. On first
## !! run, inspect rownames(fit$summary.hyperpar) and adjust `pattern`.
## ---------------------------------------------------------------------
det_draws <- function(fit, n = OPTS$n_post_samples, pattern = NULL) {
  ## R-INLA uses TWO conventions for these hyperparameters, differing by
  ## one in the index:
  ##   control.family hyper spec : beta1, beta2, beta3, ...  (1-based)
  ##   summary.hyperpar rownames : "beta[0] for occupancy observations",
  ##                               "beta[1] ...", "beta[2] ..."  (0-based)
  ## In both, the FIRST one is the detection intercept.
  ##
  ## No trailing \\b: the character after "]" is a space, and both are
  ## non-word characters, so there is no word boundary there and the
  ## match fails.
  if (is.null(pattern)) pattern <- "^beta\\[[0-9]+\\]|^beta[0-9]+"

  nm  <- rownames(fit$summary.hyperpar)
  idx <- grep(pattern, nm)

  ## Guard: a NULL here produces `if (logical(0))` and the opaque error
  ## "argument is of length zero". It means the fit did not come from
  ## fit_occupancy().
  if (is.null(fit$n_det_coef)) {
    stop("det_draws(): fit$n_det_coef is NULL. This model was not fitted ",
         "through fit_occupancy(); route it through that function so the ",
         "detection-coefficient count is recorded.")
  }
  if (length(idx) != fit$n_det_coef) {
    stop("det_draws(): expected ", fit$n_det_coef,
         " detection coefficients, matched ", length(idx),
         ".\n  rownames(summary.hyperpar) = ",
         paste(nm, collapse = " | "),
         "\n  Adjust the `pattern` argument to match this INLA build.")
  }
  s <- INLA::inla.hyperpar.sample(n = n, result = fit)
  ## Match by NAME, not by position: the column order of
  ## inla.hyperpar.sample() is not contractually tied to the row order of
  ## summary.hyperpar, and a silent mismatch here would swap the
  ## detection intercept with a covariate effect.
  want <- nm[idx]
  if (!all(want %in% colnames(s))) {
    stop("det_draws(): hyperpar sample columns do not match summary rows.\n",
         "  wanted: ", paste(want, collapse = " | "), "\n",
         "  have:   ", paste(colnames(s), collapse = " | "))
  }
  out <- s[, want, drop = FALSE]
  ## In the occupancy family beta1 is the DETECTION INTERCEPT and
  ## beta2, beta3, ... are the covariate effects, so alpha0 = beta1.
  colnames(out) <- paste0("alpha", seq_len(ncol(out)) - 1L)
  out
}

## ---------------------------------------------------------------------
## Posterior draws for OCCUPANCY fixed effects.
## ---------------------------------------------------------------------
fixed_draws <- function(fit, names_wanted, n = OPTS$n_post_samples) {
  sel <- setNames(vector("list", length(names_wanted)), names_wanted)
  sel[] <- 1
  s <- INLA::inla.posterior.sample(n = n, result = fit, selection = sel)
  out <- vapply(s, \(d) as.numeric(d$latent), numeric(length(names_wanted)))
  out <- if (length(names_wanted) == 1L) matrix(out, ncol = 1) else t(out)
  colnames(out) <- names_wanted
  out
}

## ---------------------------------------------------------------------
## Species-level intercept beta0 + gamma_i, with correct uncertainty.
##
## beta0 and gamma_i are NOT separately identified unless a sum-to-zero
## constraint is imposed on gamma. In the real-data script the species
## effect was fitted as f(species_id, model = "iid") with no constraint,
## which is why the posterior SD of every gamma_i came out at ~0.49 --
## essentially equal to the posterior SD of beta0 itself (0.51). That
## pattern is the signature of the two components trading off against
## each other, not of genuine species-level precision.
##
## Fit with constr = TRUE, and report this sum, which is identified.
## ---------------------------------------------------------------------
species_intercepts <- function(fit, n_species, n = OPTS$n_post_samples,
                               intercept_name = "Int_occ",
                               species_effect = "species_id") {
  sel <- list(1); names(sel) <- intercept_name
  sel[[species_effect]] <- seq_len(n_species)

  s <- INLA::inla.posterior.sample(n = n, result = fit, selection = sel)

  ## Index the latent vector BY NAME. inla.posterior.sample() does not
  ## guarantee that components appear in the order given in `selection`,
  ## and an earlier version of this function assumed it did: the first
  ## element was taken as the intercept and the rest as the species
  ## effects in order. That produced values with no relationship to
  ## species identity -- the community intercept plus a species effect
  ## correlated 0.98 with naive occupancy when computed directly, and not
  ## at all when computed by position.
  nm <- rownames(s[[1]]$latent)
  i_int <- grep(paste0("^", intercept_name, ":"), nm)
  i_sp  <- vapply(seq_len(n_species),
                  \(i) grep(paste0("^", species_effect, ":", i, "$"), nm),
                  integer(1))

  if (length(i_int) != 1L || any(is.na(i_sp))) {
    stop("species_intercepts(): could not locate components by name.\n",
         "  latent rownames: ", paste(head(nm, 20), collapse = " | "))
  }

  draws <- t(vapply(s, \(d) {
    v <- as.numeric(d$latent)
    v[i_int] + v[i_sp]
  }, numeric(n_species)))

  out <- data.frame(
    species = seq_len(n_species),
    mean    = colMeans(draws),
    sd      = apply(draws, 2, sd),
    lwr     = apply(draws, 2, quantile, 0.025),
    upr     = apply(draws, 2, quantile, 0.975)
  )

  ## Sanity check: the constrained species effects sum to zero, so the
  ## mean of beta0 + gamma_i must equal beta0.
  b0 <- fit$summary.fixed[intercept_name, "mean"]
  if (abs(mean(out$mean) - b0) > 0.05) {
    warning(sprintf(
      "species_intercepts(): mean of beta0 + gamma_i is %.3f but beta0 is %.3f. Check the indexing.",
      mean(out$mean), b0))
  }
  out
}
