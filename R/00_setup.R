## =====================================================================
## 00_setup.R -- environment, packages, reproducibility bookkeeping
## =====================================================================
## Sourced by every driver script. Does NOT fit anything.
## =====================================================================

suppressPackageStartupMessages({
  library(INLA)
  library(fmesher)
  library(sf)
  library(terra)
  library(dplyr)
  library(MASS)        # mvrnorm
  library(future.apply)
})

## ---------------------------------------------------------------------
## INLA version gate.
## The 'occupancy' likelihood family is comparatively recent and its
## hyperparameter naming has changed between releases. Pin and record.
##
## NOTE: do NOT use inla.list.models("likelihood") for this test. That
## function PRINTS the catalogue and returns NULL invisibly, so
##   "occupancy" %in% inla.list.models("likelihood")
## is always FALSE regardless of the build. Query the model list itself.
## ---------------------------------------------------------------------
.INLA_VERSION <- as.character(utils::packageVersion("INLA"))
message("R-INLA version: ", .INLA_VERSION)

.LIK <- INLA::inla.models()$likelihood
if (!"occupancy" %in% names(.LIK)) {
  stop("This INLA build does not expose family = 'occupancy'. ",
       "Update from https://inla.r-inla-download.org/R/stable")
}

## ---------------------------------------------------------------------
## Hard ceiling on the number of detection coefficients.
## The occupancy family declares a fixed number of hyperparameters
## (theta1 ... thetaN). Every detection coefficient -- intercept plus
## covariates -- consumes one.
##
## The ceiling itself belongs to this implementation, not to the
## INLA-SPDE approach, and could be raised. The constraint it reflects
## would remain: detection coefficients enter as hyperparameters rather
## than as latent field elements, and integration over that vector
## becomes infeasible as its dimension grows (Belmont et al. 2024), so
## the approximation degrades well before the ceiling is reached.
##
## Relevant to Contribution 1: each additional co-occurring species
## costs one slot.
## ---------------------------------------------------------------------
MAX_DET_COEF <- length(.LIK$occupancy$hyper)
message("Detection-coefficient ceiling (occupancy family): ", MAX_DET_COEF)

## Resolve the reported hyperparameter labels once, so det_draws() does
## not have to guess. Falls back to the "beta<i>" convention used in the
## official occupancy documentation.
DET_HYPER_NAMES <- vapply(.LIK$occupancy$hyper, function(x) {
  if (!is.null(x$short.name)) x$short.name else NA_character_
}, character(1))
if (all(is.na(DET_HYPER_NAMES))) {
  DET_HYPER_NAMES <- paste0("beta", seq_len(MAX_DET_COEF))
}
message("Detection hyperparameter labels: ",
        paste(head(DET_HYPER_NAMES, 4), collapse = ", "), ", ...")

## ---------------------------------------------------------------------
## Global inference settings.
## int.strategy: "ccd" integrates over the hyperparameters.
##   The previous pipeline used "eb" (Empirical Bayes), which FIXES the
##   hyperparameters at their posterior mode. Under "eb" every reported
##   95% credible interval understates uncertainty, so coverage is not
##   interpretable. Do not revert to "eb" for any result that is reported.
## ---------------------------------------------------------------------
OPTS <- list(
  int_strategy   = "ccd",
  num_threads    = "1:1",   # required for reproducible inla.qsample()
  n_post_samples = 2000,    # posterior draws for CRPS / derived quantities
  verbose        = FALSE
)

inla.setOption(num.threads = OPTS$num_threads)

## ---------------------------------------------------------------------
## future::multisession workers are FRESH R sessions. They inherit
## exported globals but NOT package attachments or inla.setOption()
## settings. Without this, workers silently run multi-threaded and
## inla.qsample(seed = ...) stops being reproducible.
## Call at the top of every function passed to future_lapply().
## ---------------------------------------------------------------------
init_worker <- function() {
  suppressPackageStartupMessages({
    library(INLA); library(fmesher); library(sf); library(terra)
    library(dplyr); library(MASS)
  })
  INLA::inla.setOption(num.threads = "1:1")
  invisible(NULL)
}

## ---------------------------------------------------------------------
## Prior specification -- ONE object, used by every model in the paper.
## The previous pipeline used prec = 1 in some scripts and prec = 1/2.72
## in others, and different PC priors for simulation vs. real data, with
## none of it documented in the manuscript. Anything reported in the
## paper must be traceable to this list.
##
## NOTE ON prior.range: in the simulation the PC prior median must NOT be
## placed at the true range. Doing so flatters the method. Set
## prior.range deliberately away from truth and run PRIORS_SENS as well.
## ---------------------------------------------------------------------
PRIORS <- list(
  ## fixed effects on the occupancy linear predictor: N(0, 1/prec)
  fixed        = list(prec = 1, prec.intercept = 1),

  ## detection coefficients enter as hyperparameters of the occupancy
  ## family: normal(mean, precision)
  det_beta     = list(mean = 0, prec = 1 / 3),

  ## PC prior for the SPDE: P(range < r0) = pr, P(sigma > s0) = ps
  spde_sim     = list(range = c(150, 0.5), sigma = c(1.0, 0.5)),
  spde_real    = list(range = c(5,   0.7), sigma = c(1.0, 0.5)),

  ## PC prior for AR(1) correlation, base model rho = 0
  ar1          = list(rho = list(prior = "pc.cor0", param = c(0.5, 0.3))),

  ## PC prior for the species-effect standard deviation
  sigma_gamma  = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
)

## An alternative prior set for the sensitivity analysis (Section on
## prior sensitivity). Same structure, deliberately different.
PRIORS_SENS <- modifyList(PRIORS, list(
  fixed    = list(prec = 1 / 4, prec.intercept = 1 / 4),
  det_beta = list(mean = 0, prec = 1 / 10),
  spde_sim = list(range = c(50, 0.5), sigma = c(0.5, 0.5))
))

## ---------------------------------------------------------------------
## Write the session record next to the results. Every table in the
## manuscript must be reproducible from this file plus a seed.
## ---------------------------------------------------------------------
## `tag` distinguishes runs that share an output directory. Without it,
## the three multi-species runs -- which all write to results/msom --
## each overwrote the previous run's session record, so only the last
## one survived. The manuscript states that session information is
## recorded for EACH set of results, so every script that shares a
## directory with another must pass a tag.
record_session <- function(outdir, tag = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sfx <- if (is.null(tag) || !nzchar(tag)) "" else paste0("_", tag)

  session_file <- file.path(outdir, paste0("sessionInfo", sfx, ".txt"))
  priors_file  <- file.path(outdir, paste0("priors_and_options", sfx, ".rds"))

  ## Overwriting a session record silently is the failure this argument
  ## exists to prevent, so say so rather than doing it quietly.
  if (file.exists(session_file)) {
    message("record_session(): overwriting ", session_file)
  }

  writeLines(
    c(paste("date:", Sys.time()),
      paste("INLA:", .INLA_VERSION),
      paste("int.strategy:", OPTS$int_strategy),
      if (nzchar(sfx)) paste("run tag:", tag) else NULL,
      capture.output(sessionInfo())),
    session_file
  )
  saveRDS(list(PRIORS = PRIORS, OPTS = OPTS, tag = tag), priors_file)
  invisible(NULL)
}

## ---------------------------------------------------------------------
## Spatio-temporal SPDE field for use with f(..., A.local = ).
##
## Convention determined EMPIRICALLY (R/99_diag_group.R), against a
## simulated AR(1) field with rho = 0.80:
##
##   A grouped + obs-level group   GroupRho 0.794  WAIC 473.4   CORRECT
##   A grouped + ngroup=, no group GroupRho 0.794  WAIC 473.4   equivalent
##   A grouped + make.index        fails: index object not found
##   A spatial + obs-level group   GroupRho 0.000  WAIC 949.4   BROKEN
##
## So: the grouping is carried by the COLUMNS of A, which must be built
## with group= and have spde$n.spde * n_time columns. The f() index is a
## length-n_obs placeholder. The inla.stack idiom (inla.spde.make.index)
## does NOT apply to A.local.
##
## The last variant is the dangerous one: it runs without error, returns
## a GroupRho, and has silently switched the temporal structure off.
##
## !! INLA emits
## !!   "Number of unused groups >= the number of groups used: k >= 1"
## !! for EVERY variant above, including the correct one. It is a FALSE
## !! POSITIVE with A.local: INLA counts groups from the index vector,
## !! which is all NA by design. Do not "fix" it. Note it in the
## !! manuscript's implementation subsection so a reader who reruns the
## !! code is not misled.
## ---------------------------------------------------------------------
make_st_field <- function(mesh, spde, coords, time, n_time,
                          name = "spatialfield") {
  stopifnot(ncol(coords) == 2, length(time) == nrow(coords),
            all(time %in% seq_len(n_time)))

  A <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords),
                              group = time, n.group = n_time)
  if (ncol(A) != spde$n.spde * n_time) {
    stop(sprintf("make_st_field(): A has %d columns, expected %d = n.spde * n_time.",
                 ncol(A), spde$n.spde * n_time))
  }

  idx <- list(rep(NA_real_, nrow(A)), time)
  names(idx) <- c(name, paste0(name, ".group"))

  list(A = A, idx = idx, n_field = ncol(A))
}

## ---------------------------------------------------------------------
## Run a fit while recording INLA's warnings. The unused-groups warning
## is expected and benign (above); anything else should be looked at.
## ---------------------------------------------------------------------
BENIGN_WARNINGS <- "unused groups"

with_warning_log <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  unexpected <- w[!grepl(paste(BENIGN_WARNINGS, collapse = "|"), w)]
  if (length(unexpected)) {
    for (u in unexpected) warning("INLA: ", u, call. = FALSE)
  }
  attr(val, "warnings") <- w
  val
}

## ---------------------------------------------------------------------
## Defensive helper: assert that two vectors are row-aligned.
## The single-species pipeline previously assigned a covariate with
##   x_s_vector <- as.vector(x_covariate); mutate(x_s = x_s_vector)
## where length() happened to match but the ordering did not. Every
## covariate join in this rewrite goes through extract_at() below.
## ---------------------------------------------------------------------
extract_at <- function(rast_layer, coords, name) {
  stopifnot(is.matrix(coords) || is.data.frame(coords), ncol(coords) == 2)
  v <- terra::extract(rast_layer, as.matrix(coords))[, 1]
  if (anyNA(v)) {
    stop(sprintf("extract_at('%s'): %d NA values -- points fall outside the raster.",
                 name, sum(is.na(v))))
  }
  v
}
