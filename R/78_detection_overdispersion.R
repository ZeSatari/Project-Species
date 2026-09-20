## =====================================================================
## 78_detection_overdispersion.R -- are detection counts binomial?
## =====================================================================
## PaperB, Section 5.3, first caveat, and Appendix S3.6.
##
## Scripts 68 and 69 show that the conditional (zero-truncated binomial)
## estimate of per-visit detection is lower than the independent fits'
## estimate for all twelve species. Both are consistent under a correctly
## specified constant-detection model, so a common sign suggests a shared
## departure from it. Detection that varies among sites, or occupancy that
## changes within a season, both make the number of detections per recorded
## site-year OVERDISPERSED relative to the zero-truncated binomial: too many
## 1s and 3s, too few 2s.
##
## Test, per species, on site-years with all three visits made and at
## least one detection:
##   Y+ ~ ZTB(3, p),  P(Y+ = y) = C(3,y) p^y (1-p)^(3-y) / (1 - (1-p)^3)
##   mean = 3p / (1 - q^3),  E[Y+^2] = (3pq + 9p^2) / (1 - q^3),  q = 1 - p
## Dispersion index D = sample variance / ZTB variance at the fitted p.
## D ~ 1 under the model; D > 1 indicates overdispersion. The p-value is a
## parametric bootstrap (refitting p in each replicate), because expected
## cell counts are small for the rare species and chi-square is unreliable.
##
## Checked in Python: the ZTB moments above, and that D averages 1.00 under
## the model and 1.12 under site-level heterogeneity (logit SD 1).
## NOT RUN IN R BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/real"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
B    <- 2000L
SEED <- 20260201L

y  <- hbefTrends$y                        # [species, site, year, visit]
SP <- dimnames(y)[[1]]
K  <- dim(y)[4]
stopifnot(K == 3L)

ztb_probs <- function(p, K = 3L) {
  yy <- seq_len(K)
  dbinom(yy, K, p) / (1 - (1 - p)^K)
}

## MLE of p = moment match of the ZTB mean (one-parameter exponential family)
fit_p <- function(ybar, K = 3L) {
  if (ybar <= 1 + 1e-10) return(NA_real_)   # all singletons: p -> 0, unidentified
  if (ybar >= K - 1e-10) return(1)
  uniroot(function(p) K * p / (1 - (1 - p)^K) - ybar,
          c(1e-10, 1 - 1e-10), tol = 1e-12)$root
}

disp_index <- function(yv, K = 3L) {
  p <- fit_p(mean(yv), K)
  if (is.na(p) || p >= 1) return(c(p = p, D = NA_real_))
  pr <- ztb_probs(p, K); yy <- seq_len(K)
  m  <- sum(pr * yy); v <- sum(pr * yy^2) - m^2
  c(p = p, D = var(yv) / v)
}

set.seed(SEED)
res <- do.call(rbind, lapply(SP, function(sp) {
  a    <- y[sp, , , , drop = TRUE]                    # [site, year, visit]
  Kst  <- apply(!is.na(a), c(1, 2), sum)
  Ypl  <- apply(a, c(1, 2), function(v) sum(v, na.rm = TRUE))
  yv   <- Ypl[Kst == 3L & Ypl >= 1]                   # recorded, 3 visits made
  n    <- length(yv)
  obs  <- tabulate(yv, nbins = 3L)
  st   <- disp_index(yv)

  if (is.na(st["D"]) || n < 10) {
    boot_p <- NA_real_; ex <- rep(NA_real_, 3)
  } else {
    pr <- ztb_probs(st["p"])
    ex <- n * pr
    Db <- replicate(B, disp_index(sample.int(3L, n, TRUE, prob = pr))["D"])
    Db <- Db[is.finite(Db)]
    boot_p <- (1 + sum(Db >= st["D"])) / (1 + length(Db))
  }
  data.frame(species = sp, n_rec3 = n, p_hat3 = unname(st["p"]),
             obs1 = obs[1], obs2 = obs[2], obs3 = obs[3],
             exp1 = ex[1], exp2 = ex[2], exp3 = ex[3],
             D = unname(st["D"]), p_boot = boot_p)
}))

res <- res[order(res$n_rec3), ]
cat("\n", strrep("=", 90), "\n", sep = "")
print(res, row.names = FALSE, digits = 3)

ok <- is.finite(res$D)
cat(sprintf("\nSpecies with D > 1: %d of %d; with bootstrap p < 0.05: %d\n",
            sum(res$D[ok] > 1), sum(ok), sum(res$p_boot[ok] < 0.05)))
cat("Overdispersion appears as obs1 and obs3 above exp1 and exp3, obs2 below exp2.\n")
cat("p-values are per species and not adjusted for 12 tests; D is the summary to report.\n")

write.csv(res, file.path(OUTDIR, "detection_overdispersion.csv"), row.names = FALSE)
