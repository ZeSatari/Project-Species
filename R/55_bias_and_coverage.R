## =====================================================================
## 55_bias_and_coverage.R -- two open weaknesses, tested cheaply
## =====================================================================
## PaperB reports two results that are weaknesses rather than findings:
##
##  (1) the community slope is displaced by about -0.076 with coverage
##      0.68-0.76 (Section 4.2). A stripped logistic mixed model with
##      random slopes showed no displacement (54_beta1_mechanism.R), so
##      the cause lies in what that model omits: the occupancy likelihood
##      or the spatial field. This script adds the occupancy layer, with
##      detection homogeneous and known in form but estimated, and still
##      no field. If the displacement appears here, the occupancy
##      likelihood is implicated and the field is not needed for it.
##
##  (2) coverage of the species-level intercepts beta0 + gamma_i falls to
##      0.70-0.82 for species far from the community mean (Table 5).
##      That is read as over-shrinkage. There is a second possible cause:
##      the reported intervals are built from the MARGINAL summaries of
##      beta0 and of gamma_i, which ignore their posterior covariance and
##      are therefore too narrow. This script computes both --- the
##      marginal-sum interval and the interval from joint posterior
##      samples --- on the same fits, so the two causes can be separated.
##
## Both run on a stripped model: no spatial field, detection either
## absent (the direct-observation case) or homogeneous across species.
## A few minutes for 200 replicates.
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(INLA))

OUTDIR <- "results/beta1"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "bias_and_coverage")

N_REP     <- 200L
N_SP      <- 12L
N_SITE    <- 200L
N_VISIT   <- 3L
BETA0     <- -0.5
BETA1     <-  0.7
SD_GAMMA  <-  0.8
SD_DELTA  <-  0.5
ALPHA0    <- -0.9          # homogeneous detection, as in run (c)
ALPHA1    <- -1.0          # on one survey covariate, as in run (c)
N_DRAWS   <- 500L
SEED_BASE <- 20260201L

## ---------------------------------------------------------------------
## One replicate: occupancy likelihood, homogeneous detection, no field
## ---------------------------------------------------------------------
one_rep <- function(r) {
  set.seed(SEED_BASE + r)
  gamma <- rnorm(N_SP, 0, SD_GAMMA); gamma <- gamma - mean(gamma)
  delta <- rnorm(N_SP, 0, SD_DELTA); delta <- delta - mean(delta)
  x     <- as.numeric(scale(rnorm(N_SITE)))

  eta <- outer(x, BETA1 + delta) +
         matrix(BETA0 + gamma, N_SITE, N_SP, byrow = TRUE)
  psi <- plogis(eta)
  z   <- matrix(rbinom(length(psi), 1, psi), N_SITE, N_SP)

  ## detection: one survey covariate, its coefficient common to all
  ## species, exactly as in run (c). build_det_X() needs at least one
  ## covariate, so an intercept-only detection model is not available.
  n_row <- N_SITE * N_SP
  g   <- matrix(runif(n_row * N_VISIT, -1, 1), n_row, N_VISIT)
  p   <- plogis(ALPHA0 + ALPHA1 * g)
  zv  <- as.vector(z)
  Y   <- matrix(rbinom(length(p), 1, as.vector(zv * p)), n_row, N_VISIT)

  dat <- data.frame(Int_occ = 1,
                    x  = rep(x, times = N_SP),
                    sp = rep(seq_len(N_SP), each = N_SITE))
  dat$sp2 <- dat$sp
  Xd <- build_det_X(list(g), K = N_VISIT)         # intercept + covariate

  f <- inla.mdata(Y, X) ~ -1 + Int_occ + x +
    f(sp,  model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(sp2, x, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma)
  fit <- fit_occupancy(f, c(as.list(dat), list(Y = Y, X = Xd)), Xd,
                       compute_cv = FALSE)

  ## --- (1) the community slope
  b1 <- fit$summary.fixed["x", ]

  ## --- (2) species intercepts, two ways
  fx <- fit$summary.fixed["Int_occ", ]
  sr <- fit$summary.random$sp
  truth_i <- BETA0 + gamma

  ## marginal-sum interval: add the two marginals, ignoring covariance,
  ## which is how the reported intervals were built
  sd_marg <- sqrt(fx$sd^2 + sr$sd^2)
  lo_marg <- (fx$mean + sr$mean) - 1.96 * sd_marg
  hi_marg <- (fx$mean + sr$mean) + 1.96 * sd_marg

  ## joint-sample interval: the posterior of the sum itself
  s  <- inla.posterior.sample(n = N_DRAWS, result = fit,
                              seed = SEED_BASE + r, num.threads = "1:1")
  nm <- rownames(s[[1]]$latent)
  i_int <- grep("^Int_occ:", nm); stopifnot(length(i_int) == 1L)
  i_sp  <- vapply(seq_len(N_SP), \(i) {
    j <- grep(sprintf("^sp:%d$", i), nm); stopifnot(length(j) == 1L); j
  }, 1L)
  draws <- vapply(s, \(d) d$latent[i_int] + d$latent[i_sp], numeric(N_SP))
  lo_joint <- apply(draws, 1, quantile, 0.025)
  hi_joint <- apply(draws, 1, quantile, 0.975)

  data.frame(
    rep      = r,
    species  = seq_len(N_SP),
    truth    = truth_i,
    dist     = abs(gamma),                       # distance from the mean
    beta1_hat = b1$mean, beta1_lo = b1$`0.025quant`, beta1_hi = b1$`0.975quant`,
    beta1_true = mean(BETA1 + delta),
    cov_marg  = as.integer(truth_i >= lo_marg  & truth_i <= hi_marg),
    cov_joint = as.integer(truth_i >= lo_joint & truth_i <= hi_joint),
    width_marg  = hi_marg - lo_marg,
    width_joint = hi_joint - lo_joint)
}

res <- do.call(rbind, lapply(seq_len(N_REP), function(r) {
  if (r %% 20 == 0) message("replicate ", r)
  tryCatch(one_rep(r), error = function(e) {
    message("  replicate ", r, " failed: ", conditionMessage(e)); NULL })
}))
saveRDS(res, file.path(OUTDIR, "bias_and_coverage.rds"))

## ---------------------------------------------------------------------
## (1) the community slope
## ---------------------------------------------------------------------
if (is.null(res) || !nrow(res))
  stop("every replicate failed; fix the error printed above before reading on")
b <- res[!duplicated(res$rep), ]
cat("\n", strrep("=", 74), "\n(1) COMMUNITY SLOPE, occupancy likelihood, no field\n", sep = "")
cat(sprintf("replicates %d\n", nrow(b)))
cat(sprintf("bias      %+.4f (MC SE %.4f)\n",
            mean(b$beta1_hat - b$beta1_true),
            sd(b$beta1_hat - b$beta1_true) / sqrt(nrow(b))))
cat(sprintf("coverage  %.3f (SE %.3f)\n",
            mean(b$beta1_lo <= b$beta1_true & b$beta1_true <= b$beta1_hi),
            sqrt(0.25 / nrow(b))))
cat("Compare with run (a): bias -0.076, coverage 0.68; and with the\n")
cat("field-free logistic model of 54_beta1_mechanism.R: bias +0.002.\n")
cat("A displacement here implicates the occupancy likelihood; its absence\n")
cat("leaves the spatial field as the remaining candidate.\n")

## ---------------------------------------------------------------------
## (2) species-level coverage, marginal vs joint intervals
## ---------------------------------------------------------------------
cat("\n", strrep("=", 74), "\n(2) SPECIES INTERCEPTS: how the interval is built\n", sep = "")
cat(sprintf("coverage, marginal-sum intervals: %.3f\n", mean(res$cov_marg)))
cat(sprintf("coverage, joint-sample intervals: %.3f\n", mean(res$cov_joint)))
cat(sprintf("mean width, marginal %.3f vs joint %.3f (ratio %.2f)\n",
            mean(res$width_marg), mean(res$width_joint),
            mean(res$width_joint) / mean(res$width_marg)))

res$terc <- cut(res$dist, quantile(res$dist, c(0, 1/3, 2/3, 1)),
                include.lowest = TRUE,
                labels = c("near mean", "middle", "far from mean"))
tab <- do.call(rbind, lapply(split(res, res$terc), function(d) data.frame(
  tercile = as.character(d$terc[1]), n = nrow(d),
  cov_marginal = mean(d$cov_marg), cov_joint = mean(d$cov_joint),
  width_marginal = mean(d$width_marg), width_joint = mean(d$width_joint))))
print(tab, row.names = FALSE, digits = 3)
cat("\nIf the joint intervals cover at the nominal rate while the marginal\n")
cat("ones do not, the under-coverage reported in Table 5 is an artefact of\n")
cat("how the intervals were formed, not evidence of over-shrinkage.\n")

write.csv(res, file.path(OUTDIR, "bias_and_coverage.csv"), row.names = FALSE)
cat("\nwritten to ", OUTDIR, "\n", sep = "")
