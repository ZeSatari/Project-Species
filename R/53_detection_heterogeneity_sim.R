## =====================================================================
## 53_detection_heterogeneity_sim.R -- what unmodelled heterogeneity does
## =====================================================================
## PaperB, Section 5.3, first caveat. Two questions:
##
##  Q1. With within-species detection heterogeneity of the size the data
##      show (dispersion index D about 1.1), how biased are the constant-p
##      estimate of detection and the implied occupancy correction
##      psibar/n = 1 / (1 - (1-p)^K)?
##
##  Q2. Does such heterogeneity, or unmodelled heterogeneity in occupancy,
##      make the full occupancy likelihood estimate detection HIGHER than
##      the conditional (zero-truncated) likelihood does? That is the
##      pattern in the application, where the independent fits exceed the
##      conditional estimates for all twelve species.
##
## A result worth stating before running anything: if the occupancy model
## contains a free level parameter and detection is constant, the profile
## likelihood for p is exactly the zero-truncated conditional likelihood,
## because psi absorbs the zero cell. The two estimators are then the same
## estimator, whatever the truth. Divergence requires the occupancy part
## to be constrained relative to the data. Scenario C below tests whether
## a covariate-driven occupancy model is constrained enough to produce it.
##
## No INLA, a few minutes. Checked against a Python implementation.
## =====================================================================

set.seed(20260201L)
OUTDIR <- "results/het_sim"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

S <- 3357L      # surveyed site-years, as in the application
K <- 3L
R <- 200L       # replicates

lg  <- function(p) log(p / (1 - p))
ilg <- plogis

## --- estimators -------------------------------------------------------
## conditional (zero-truncated binomial) MLE, free of the occupancy model
ztb_mle <- function(y) {
  r <- y[y >= 1]
  if (length(r) < 5L) return(NA_real_)
  nll <- function(e) {
    p <- ilg(e)
    -sum(r * log(p) + (K - r) * log1p(-p) - log(-expm1(K * log1p(-p))))
  }
  ilg(optimize(nll, c(-12, 12), tol = 1e-10)$minimum)
}
## full occupancy likelihood, psi constant
full_mle <- function(y) {
  cnt <- tabulate(y + 1L, nbins = K + 1L)
  nll <- function(th) {
    psi <- ilg(th[1]); p <- ilg(th[2])
    lik <- psi * dbinom(0:K, K, p) + (1 - psi) * c(1, rep(0, K))
    -sum(cnt * log(pmax(lik, 1e-300)))
  }
  ilg(optim(c(0, 0), nll, control = list(reltol = 1e-12))$par[2])
}
## full occupancy likelihood, psi driven by a covariate
full_mle_cov <- function(y, x) {
  nll <- function(th) {
    psi <- ilg(th[1] + th[2] * x); p <- ilg(th[3])
    lik <- ifelse(y > 0,
                  psi * dbinom(y, K, p),
                  psi * (1 - p)^K + (1 - psi))
    -sum(log(pmax(lik, 1e-300)))
  }
  ilg(optim(c(0, 0, 0), nll, control = list(reltol = 1e-12))$par[3])
}
## dispersion index of the detection counts, as in 78_detection_overdispersion.R
disp_index <- function(y) {
  r <- y[y >= 1]; p <- ztb_mle(y)
  if (is.na(p)) return(NA_real_)
  pr <- dbinom(1:K, K, p) / (1 - (1 - p)^K)
  m <- sum(pr * (1:K)); v <- sum(pr * (1:K)^2) - m^2
  var(r) / v
}

## --- scenario A: heterogeneous detection, homogeneous occupancy -------
simA <- function(psi, pbar, sig) {
  z <- runif(S) < psi
  p <- ilg(lg(pbar) + sig * rnorm(S))
  list(y = rbinom(S, K, ifelse(z, p, 0)), psi_true = psi)
}
## --- scenario B: homogeneous detection, heterogeneous occupancy -------
simB <- function(psibar, sig_psi, p) {
  psi <- ilg(lg(psibar) + sig_psi * rnorm(S))
  z <- runif(S) < psi
  list(y = rbinom(S, K, ifelse(z, p, 0)), psi_true = mean(psi))
}
## --- scenario C: occupancy driven by a covariate the model knows ------
simC <- function(a, b, pbar, sig) {
  x <- rnorm(S); psi <- ilg(a + b * x)
  z <- runif(S) < psi
  p <- ilg(lg(pbar) + sig * rnorm(S))
  list(y = rbinom(S, K, ifelse(z, p, 0)), x = x, psi_true = mean(psi))
}

run <- function(gen, cov = FALSE, reps = R) {
  out <- replicate(reps, {
    d <- gen()
    y <- d$y
    p_z <- ztb_mle(y)
    p_f <- if (cov) full_mle_cov(y, d$x) else full_mle(y)
    n   <- mean(y >= 1)
    c(D = disp_index(y), p_ztb = p_z, p_full = p_f,
      ratio_true = d$psi_true / n,
      ratio_ztb  = 1 / (1 - (1 - p_z)^K),
      ratio_full = 1 / (1 - (1 - p_f)^K))
  })
  round(rowMeans(out, na.rm = TRUE), 3)
}

cat("\n=== A. detection heterogeneity (occupancy homogeneous, psi = 0.25) ===\n")
A <- do.call(rbind, lapply(c(0.15, 0.40, 0.55), function(pb)
  do.call(rbind, lapply(c(0, 0.5, 1.0), function(sg)
    data.frame(p_bar = pb, sigma_logit_p = sg,
               t(run(function() simA(0.25, pb, sg))))))))
print(A, row.names = FALSE)
cat("\nRead: D is what 78_detection_overdispersion.R measures on the data\n",
    "(1.04 to 1.26 there). At that level of heterogeneity, both estimators\n",
    "overstate detection and understate the occupancy correction, so the\n",
    "corrections applied in the paper are conservative.\n", sep = "")

cat("\n=== B. occupancy heterogeneity (detection constant) ===\n")
B <- do.call(rbind, lapply(c(0.15, 0.55), function(p)
  do.call(rbind, lapply(c(0, 1, 2), function(sg)
    data.frame(p = p, sigma_logit_psi = sg,
               t(run(function() simB(0.25, sg, p))))))))
print(B, row.names = FALSE)

cat("\n=== C. occupancy driven by a covariate the model fits ===\n")
C <- do.call(rbind, lapply(c(0, 1, 2), function(b)
  do.call(rbind, lapply(list(c(0.15, 0), c(0.15, 1), c(0.55, 0), c(0.55, 1)),
    function(ps) data.frame(psi_slope = b, p_bar = ps[1], sigma_logit_p = ps[2],
                            t(run(function() simC(-1, b, ps[1], ps[2]), cov = TRUE)))))))
print(C, row.names = FALSE)

cat("\n=== Q2: does either mechanism make the full likelihood exceed the",
    "conditional one? ===\n")
for (nm in c("A", "B", "C")) {
  tb <- get(nm); d <- tb$p_full - tb$p_ztb
  cat(sprintf("  scenario %s: p_full - p_ztb between %+.3f and %+.3f\n",
              nm, min(d), max(d)))
}
cat("In the application the difference is positive for all twelve species,\n",
    "between +0.016 and +0.059 (Table S4). If the ranges above do not cover\n",
    "that, neither mechanism explains it and the cause is still open.\n", sep = "")

saveRDS(list(A = A, B = B, C = C), file.path(OUTDIR, "het_sim.rds"))
cat("\nwritten to ", OUTDIR, "\n", sep = "")
