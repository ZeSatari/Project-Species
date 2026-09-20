## =====================================================================
## 54_beta1_mechanism.R -- what is the community slope actually estimating?
## =====================================================================
## PaperB, Section 4.2. The community slope is displaced downward
## (-0.076 in run (a), -0.074 with homogeneous detection, -0.050 with
## homogeneous slopes). Section 4.2 offers a conjecture: with the species
## slopes delta_i shrunk toward zero, the posterior for beta1 behaves like
## an information-weighted rather than an unweighted average of the
## species slopes, and since the Fisher information for a logistic slope,
## sum psi(1-psi)x^2, falls as the slope grows in magnitude, species with
## steeper slopes get less weight and the average is pulled down.
##
## The conjecture concerns the ESTIMATOR, not occupancy, detection or the
## field, so it can be tested in a stripped model: a logistic mixed model
## with random intercepts and random slopes, the states observed
## directly, no detection process, no spatial field. Each fit takes
## seconds rather than the half hour a full occupancy replicate needs.
##
## Per replicate we record
##   beta1_hat      the fitted community slope
##   mean_unw       the unweighted mean of the generating species slopes
##                  (0.7 by construction: the slopes are recentred)
##   mean_inf       the information-weighted mean, weights
##                  w_i = sum_s psi_is (1 - psi_is) x_s^2 at the TRUE psi
##   sd_delta_hat   posterior SD of the species-slope effect, a measure of
##                  how much the delta_i were shrunk
##
## Reading the output:
##  * if beta1_hat tracks mean_inf and not mean_unw, the conjecture holds
##    and the displacement is a property of the estimator in hierarchical
##    logistic models, with nothing specific to occupancy;
##  * if beta1_hat is unbiased for mean_unw here, the displacement needs
##    the occupancy likelihood, the detection process or the field, and
##    the conjecture in Section 4.2 should be dropped.
##
## The stripped model has no shrinkage-free comparison built in, so we
## also fit each replicate with the species slopes as FIXED effects, which
## estimates the unweighted mean without shrinkage.
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT. Expect a few minutes.
## =====================================================================

source("R/00_setup.R")
suppressPackageStartupMessages(library(INLA))

OUTDIR  <- "results/beta1"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "beta1_mechanism")

N_REP     <- 200L
N_SP      <- 12L
N_SITE    <- 200L          # as in the multi-species simulation
BETA0     <- -0.5
BETA1     <-  0.7
SD_GAMMA  <-  0.8
SD_DELTA  <-  0.5          # as in run (a); 0 reproduces run (b)
SEED_BASE <- 20260201L

one_rep <- function(r, sd_delta = SD_DELTA) {
  set.seed(SEED_BASE + r)
  gamma <- rnorm(N_SP, 0, SD_GAMMA); gamma <- gamma - mean(gamma)
  delta <- rnorm(N_SP, 0, sd_delta); delta <- delta - mean(delta)
  x     <- as.numeric(scale(rnorm(N_SITE)))

  ## true occupancy probabilities, and the state observed without error
  eta <- outer(x, BETA1 + delta) +
         matrix(BETA0 + gamma, N_SITE, N_SP, byrow = TRUE)
  psi <- plogis(eta)
  z   <- matrix(rbinom(length(psi), 1, psi), N_SITE, N_SP)

  ## the two candidate targets
  w        <- colSums(psi * (1 - psi) * x^2)      # Fisher information per species
  mean_unw <- mean(BETA1 + delta)
  mean_inf <- sum(w * (BETA1 + delta)) / sum(w)

  dat <- data.frame(y = as.vector(z),
                    x = rep(x, times = N_SP),
                    sp = rep(seq_len(N_SP), each = N_SITE))
  dat$sp2 <- dat$sp

  ## (i) the fitted model of the paper, in stripped form
  f1 <- y ~ 1 + x +
    f(sp,  model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(sp2, x, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma)
  m1 <- inla(f1, family = "binomial", data = dat,
             control.fixed = list(mean = 0, prec = 1),
             control.compute = list(config = FALSE))

  ## (ii) species slopes as fixed effects: no shrinkage, so the estimand
  ##      is the unweighted mean by construction
  dat$spf <- factor(dat$sp)
  m2 <- inla(y ~ 1 + x + spf + x:spf, family = "binomial", data = dat,
             control.fixed = list(mean = 0, prec = 1))
  fx <- m2$summary.fixed
  slope_i <- fx["x", "mean"] +
    c(0, fx[grep("^x:spf", rownames(fx)), "mean"])

  data.frame(rep = r,
             beta1_hat    = m1$summary.fixed["x", "mean"],
             beta1_sd     = m1$summary.fixed["x", "sd"],
             mean_unw     = mean_unw,
             mean_inf     = mean_inf,
             fixed_mean   = mean(slope_i),
             sd_delta_hat = 1 / sqrt(m1$summary.hyperpar[
               grep("sp2", rownames(m1$summary.hyperpar))[1], "mean"]))
}

res <- do.call(rbind, lapply(seq_len(N_REP), function(r) {
  if (r %% 20 == 0) message("replicate ", r)
  tryCatch(one_rep(r), error = function(e) {
    message("  replicate ", r, " failed: ", conditionMessage(e)); NULL })
}))
saveRDS(res, file.path(OUTDIR, "beta1_mechanism.rds"))

cat("\n", strrep("=", 74), "\n", sep = "")
cat(sprintf("%d replicates; random slopes with SD %.2f\n", nrow(res), SD_DELTA))
cat(sprintf("\nbias against the unweighted mean (the estimand of the fitted model): %+.4f (SE %.4f)\n",
            mean(res$beta1_hat - res$mean_unw),
            sd(res$beta1_hat - res$mean_unw) / sqrt(nrow(res))))
cat(sprintf("bias against the information-weighted mean:                          %+.4f (SE %.4f)\n",
            mean(res$beta1_hat - res$mean_inf),
            sd(res$beta1_hat - res$mean_inf) / sqrt(nrow(res))))
cat(sprintf("\nthe two targets differ by %+.4f on average (unweighted minus weighted)\n",
            mean(res$mean_unw - res$mean_inf)))
cat(sprintf("fixed-effect fit, bias against the unweighted mean:                  %+.4f (SE %.4f)\n",
            mean(res$fixed_mean - res$mean_unw),
            sd(res$fixed_mean - res$mean_unw) / sqrt(nrow(res))))

cat("\nDoes the displacement follow the gap between the two targets?\n")
d <- res$beta1_hat - res$mean_unw
g <- res$mean_inf  - res$mean_unw
print(summary(lm(d ~ g))$coefficients, digits = 3)
cat("A slope near 1 with a small intercept supports the conjecture;\n")
cat("a slope near 0 with a clearly negative intercept refutes it.\n")

cat(sprintf("\nposterior SD of the species-slope effect: median %.3f (generating %.2f)\n",
            median(res$sd_delta_hat, na.rm = TRUE), SD_DELTA))
write.csv(res, file.path(OUTDIR, "beta1_mechanism.csv"), row.names = FALSE)
cat("\nwritten to ", OUTDIR, "\n", sep = "")
