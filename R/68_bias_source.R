## =====================================================================
## 68_bias_source.R -- where does the negative slope bias come from?
## =====================================================================
## Both simulations show the environmental slope biased downward: about
## -0.024 in the single-species study and -0.050 in the multi-species
## one, neither diminishing with sample size. Every diagnostic so far has
## used the full model, so the likelihood itself has never been tested in
## isolation.
##
## This script strips everything away. One species, no latent field, no
## random effects, no constraint:
##     logit(psi_s) = beta0 + beta1 * x_s
##     logit(p_sj)  = alpha0
## If the bias survives that, it belongs to the occupancy mixture
## likelihood and not to any structure built on top of it.
##
## Four axes, varied one at a time from a baseline:
##
##   n         does it shrink with sample size? Finite-sample bias in a
##             mixture likelihood should; an approximation error or a
##             pseudo-true limit should not.
##   detection how does it depend on p? At p = 1 the occupancy state is
##             observed and the model reduces to logistic regression, so
##             any bias there is not about imperfect detection at all.
##   visits    more visits per site means the state is better resolved.
##   strategy  Gaussian, simplified Laplace, full Laplace. If the bias
##             moves with the approximation it is INLA's; if not, it is
##             the estimator's.
##
## The p = 1 cell carries a further check: the same data fitted by glm(),
## which has no approximation and no mixture. Agreement there separates
## "INLA is approximating badly" from "this estimator is biased".
##
## R replicates per cell, no spatial structure: minutes, not hours.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

OUTDIR <- "results/bias"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR)

R_REP  <- 200
SEED0  <- 20260501L

TRUTH <- list(beta0 = -0.5, beta1 = 0.7, alpha0 = qlogis(0.30))

## ---------------------------------------------------------------------
## Generator: the occupancy model with nothing else in it.
## ---------------------------------------------------------------------
sim_flat <- function(n, K, p_det, seed, beta1 = TRUTH$beta1) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, plogis(TRUTH$beta0 + beta1 * x))
  Y <- matrix(rbinom(n * K, 1, z * p_det), n, K)
  list(Y = Y, x = x, z = z)
}

## build_det_X() requires at least one covariate, and that guard has
## caught real errors elsewhere, so it is not relaxed. The detection model
## here has an intercept and nothing else, so the design is assembled
## directly in the same layout build_det_X() produces: K blocks of one
## column each, with the attributes fit_occupancy() reads.
det_X_intercept_only <- function(n, K) {
  X <- matrix(1, nrow = n, ncol = K)
  attr(X, "n_det_coef") <- 1L
  attr(X, "na_action")  <- "zero"
  X
}

fit_inla <- function(d, strategy = "simplified.laplace") {
  n  <- nrow(d$Y); K <- ncol(d$Y)
  Xd <- det_X_intercept_only(n, K)
  dat <- list(Y = d$Y, X = Xd, Int_occ = rep(1, n), xx = d$x)
  f <- inla.mdata(Y, X) ~ -1 + Int_occ + xx
  fit <- tryCatch(fit_occupancy(f, dat, Xd, strategy = strategy,
                                compute_cv = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(beta1 = NA_real_, covered = NA_real_))
  s <- fit$summary.fixed["xx", ]
  c(beta1 = s[["mean"]],
    covered = as.numeric(s[["0.025quant"]] <= TRUTH$beta1 &
                         TRUTH$beta1 <= s[["0.975quant"]]))
}

## With p = 1 the latent state is observed: occupancy is exactly
## "detected at least once", and ordinary logistic regression applies.
fit_glm <- function(d) {
  y <- as.integer(rowSums(d$Y) > 0)
  m <- glm(y ~ d$x, family = binomial)
  ci <- suppressMessages(confint(m))
  c(beta1 = unname(coef(m)[2]),
    covered = as.numeric(ci[2, 1] <= TRUTH$beta1 & TRUTH$beta1 <= ci[2, 2]))
}

run_cell <- function(label, n, K, p_det, strategy = "simplified.laplace",
                     use_glm = FALSE, R = R_REP) {
  message(sprintf("  %-34s n=%5d K=%d p=%.2f %s", label, n, K, p_det,
                  if (use_glm) "[glm]" else strategy))
  est <- vapply(seq_len(R), function(r) {
    d <- sim_flat(n, K, p_det, SEED0 + r)
    if (use_glm) fit_glm(d) else fit_inla(d, strategy)
  }, numeric(2))
  b <- est["beta1", ]; ok <- !is.na(b)
  data.frame(cell = label, n = n, K = K, p_det = p_det,
             method = if (use_glm) "glm" else strategy,
             n_ok = sum(ok),
             bias = mean(b[ok]) - TRUTH$beta1,
             mc_se = sd(b[ok]) / sqrt(sum(ok)),
             coverage = mean(est["covered", ok]))
}

## A single-replicate smoke test before committing to 3200 fits: if the
## intercept-only design is not accepted, or alpha0 is not recovered, the
## rest of the run is worthless.
{
  d0 <- sim_flat(2000, 3, 0.30, SEED0)
  Xd0 <- det_X_intercept_only(nrow(d0$Y), ncol(d0$Y))
  dat0 <- list(Y = d0$Y, X = Xd0, Int_occ = rep(1, nrow(d0$Y)), xx = d0$x)
  f0 <- inla.mdata(Y, X) ~ -1 + Int_occ + xx
  fit0 <- fit_occupancy(f0, dat0, Xd0, compute_cv = FALSE)
  a0 <- det_draws(fit0, n = 500)
  message(sprintf("smoke test: beta1 %.3f (truth %.3f), alpha0 %.3f (truth %.3f)",
                  fit0$summary.fixed["xx", "mean"], TRUTH$beta1,
                  mean(a0[, 1]), TRUTH$alpha0))
  if (abs(mean(a0[, 1]) - TRUTH$alpha0) > 0.25)
    stop("smoke test: alpha0 not recovered; the detection design is wrong.")
  rm(d0, Xd0, dat0, f0, fit0, a0)
}

cells <- list()

## --- axis 1: sample size, baseline p = 0.30, K = 3 --------------------
message("axis 1: sample size")
for (n in c(250, 500, 1000, 2000, 5000))
  cells[[length(cells) + 1]] <- run_cell(sprintf("n = %d", n), n, 3, 0.30)

## --- axis 2: detection probability, n = 1000, K = 3 -------------------
message("axis 2: detection probability")
for (p in c(0.15, 0.30, 0.60, 0.90))
  cells[[length(cells) + 1]] <- run_cell(sprintf("p = %.2f", p), 1000, 3, p)

## p = 1: the state is observed. INLA and glm on identical data.
cells[[length(cells) + 1]] <- run_cell("p = 1.00 (INLA)", 1000, 3, 1.00)
cells[[length(cells) + 1]] <- run_cell("p = 1.00 (glm)",  1000, 3, 1.00,
                                       use_glm = TRUE)

## --- axis 3: visits per site, n = 1000, p = 0.30 ----------------------
message("axis 3: visits per site")
for (K in c(2, 3, 5, 10))
  cells[[length(cells) + 1]] <- run_cell(sprintf("K = %d", K), 1000, K, 0.30)

## --- axis 4: approximation strategy, n = 1000, K = 3, p = 0.30 --------
message("axis 4: approximation strategy")
for (s in c("gaussian", "simplified.laplace", "laplace"))
  cells[[length(cells) + 1]] <- run_cell(paste("strategy:", s), 1000, 3, 0.30,
                                         strategy = s)

out <- do.call(rbind, cells); rownames(out) <- NULL
cat("\n", strrep("=", 92), "\n", sep = "")
print(out, row.names = FALSE, digits = 4)
write.csv(out, file.path(OUTDIR, "bias_source.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## Read the axes
## ---------------------------------------------------------------------
cat("\n", strrep("-", 92), "\n", sep = "")

nn <- out[grepl("^n = ", out$cell), ]
cat("sample size: bias goes",
    paste(sprintf("%+.4f", nn$bias), collapse = " -> "), "\n")
cat(if (abs(nn$bias[nrow(nn)]) < 0.5 * abs(nn$bias[1]))
      "  -> shrinks with n: consistent with finite-sample bias.\n"
    else
      "  -> does not shrink with n: not finite-sample bias.\n")

pp <- out[grepl("^p = ", out$cell) & out$method != "glm", ]
cat("\ndetection: bias goes",
    paste(sprintf("%.2f:%+.4f", pp$p_det, pp$bias), collapse = "  "), "\n")

g <- out[out$method == "glm", ]; i1 <- out[out$cell == "p = 1.00 (INLA)", ]
if (nrow(g) && nrow(i1)) {
  cat(sprintf("\nat p = 1, INLA %+.4f vs glm %+.4f\n", i1$bias, g$bias))
  if (abs(g$bias) > 0.01)
    cat("  -> glm is biased too: the estimator, not the approximation.\n")
  else if (abs(i1$bias) > 0.01)
    cat("  -> glm is clean but INLA is not: the approximation.\n")
  else
    cat("  -> both clean at p = 1: the bias needs imperfect detection.\n")
}

ss <- out[grepl("^strategy:", out$cell), ]
cat("\nstrategy: ",
    paste(sprintf("%s %+.4f", sub("strategy: ", "", ss$cell), ss$bias),
          collapse = " | "), "\n")
cat(if (diff(range(ss$bias)) < 0.005)
      "  -> unchanged across strategies: not an artefact of the approximation.\n"
    else
      "  -> moves with the strategy: the approximation is implicated.\n")

cat("\nFor reference: the single-species simulation gave -0.024 and the\n")
cat("multi-species one -0.050, both with the full model.\n")
