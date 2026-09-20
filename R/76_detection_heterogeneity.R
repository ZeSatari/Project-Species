## =====================================================================
## 76_detection_heterogeneity.R -- is pooled detection tenable?
## =====================================================================
## PaperB, Section 2.4 and Appendix S3.5.
##
## Among site-years at which species i was recorded, the number of
## detections Y+ is zero-truncated binomial(K_st, p_i) under constant p_i
## and conditionally independent visits. psi cancels, so p_i is identified
## WITHOUT an occupancy model:
##
##   l_i(p) = sum_{recorded st} [ Y+ log p + (K_st - Y+) log(1-p)
##                                - log(1 - (1-p)^K_st) ]
##
## With K_st = K for all st this MLE equals the moment estimator
## m_K^{-1}(d_i / n_i) used in the main text.
##
## Outputs, per species: n_i, d_i, p_hat_i (conditional MLE) and SE, the
## implied ratio psibar/n = 1 / mean(1 - (1 - p_hat)^K_st), observed and
## pooled-p predicted mean detections per recorded site-year, and a
## likelihood-ratio test of a common p. The LRT treats species as
## independent given occupancy, which is only approximately true (PaperA
## shows detections are associated within surveys); read its p-value as
## approximate.
##
## Run on hbefTrends; its output populates Table 5 (tab:detection) of PaperB.
## The arithmetic was also checked in Python against simulated data.
## =====================================================================

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/real"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## Optional: pooled per-visit detection probability from the fitted
## hierarchical model, plogis(alpha0) at standardised covariates = 0.
## Leave NA to use the pooled conditional MLE below instead.
P_POOLED_FIT <- NA_real_

y  <- hbefTrends$y                      # [species, site, year, visit]
stopifnot(length(dim(y)) == 4L)
SP <- dimnames(y)[[1]]

## Ratios as reported in Table S1 of PaperB, for comparison only.
reported <- data.frame(
  species = c("NAWA","AMRE","CAWA","BAWW","BLPW","MAWA",
              "BHVI","REVI","OVEN","BLBW","BTBW","BTNW"),
  ratio_hier  = c(1.22,1.17,1.17,1.16,1.18,1.17,1.17,1.16,1.15,1.16,1.14,1.14),
  ratio_indep = c(2.74,2.03,1.74,2.98,1.17,1.24,2.01,1.07,1.07,1.10,1.07,1.10)
)

## ---------------------------------------------------------------------
## Per species: visits made, detections, recorded site-years.
## ---------------------------------------------------------------------
per_species <- lapply(SP, function(sp) {
  a    <- y[sp, , , , drop = TRUE]                 # [site, year, visit]
  Kst  <- apply(!is.na(a), c(1, 2), sum)           # visits made
  Ypl  <- apply(a, c(1, 2), function(v) sum(v, na.rm = TRUE))
  surveyed <- Kst > 0
  list(K = Kst[surveyed], Y = Ypl[surveyed])
})
names(per_species) <- SP

## Missing-visit pattern should be identical across species (PaperA).
Kmat <- sapply(SP, function(sp) {
  a <- y[sp, , , , drop = TRUE]; as.vector(apply(!is.na(a), c(1, 2), sum))
})
if (any(apply(Kmat, 1, function(r) length(unique(r)) > 1)))
  warning("Visit counts differ between species; check before interpreting.")

## ---------------------------------------------------------------------
## Conditional (zero-truncated binomial) log-likelihood on the logit scale.
## ---------------------------------------------------------------------
ztb_ll <- function(eta, K, Y) {
  p   <- plogis(eta)
  rec <- Y >= 1
  if (!any(rec)) return(NA_real_)
  K <- K[rec]; Y <- Y[rec]
  ## log(1 - (1-p)^K) computed stably
  log1mq <- log(-expm1(K * log1p(-p)))
  sum(Y * log(p) + (K - Y) * log1p(-p) - log1mq)
}

fit_p <- function(K, Y, lo = -12, hi = 12) {
  o <- optimize(function(e) -ztb_ll(e, K, Y), c(lo, hi))
  eta <- o$minimum
  ## observed information by central difference on the logit scale
  h  <- 1e-4
  d2 <- (ztb_ll(eta + h, K, Y) - 2 * ztb_ll(eta, K, Y) +
         ztb_ll(eta - h, K, Y)) / h^2
  se_eta <- if (is.finite(d2) && d2 < 0) sqrt(-1 / d2) else NA_real_
  p <- plogis(eta)
  list(eta = eta, p = p, se_p = se_eta * p * (1 - p),
       ll = -o$objective,
       at_bound = abs(eta - lo) < 0.01 || abs(eta - hi) < 0.01)
}

res <- do.call(rbind, lapply(SP, function(sp) {
  d  <- per_species[[sp]]
  f  <- fit_p(d$K, d$Y)
  rec <- d$Y >= 1
  data.frame(
    species   = sp,
    n_st      = length(d$K),
    n_rec     = sum(rec),
    naive     = mean(rec),
    det_freq  = sum(d$Y) / sum(d$K),
    ratio_dn  = (sum(d$Y) / sum(d$K)) / mean(rec),
    p_hat     = f$p,
    se_p      = f$se_p,
    at_bound  = f$at_bound,
    implied_ratio = 1 / mean(1 - (1 - f$p)^d$K),
    obs_det_per_rec = mean(d$Y[rec]),
    ll        = f$ll
  )
}))

## ---------------------------------------------------------------------
## Pooled p: conditional MLE under p_1 = ... = p_N, and the LRT.
## ---------------------------------------------------------------------
K_all <- unlist(lapply(per_species, `[[`, "K"))
Y_all <- unlist(lapply(per_species, `[[`, "Y"))
pooled <- fit_p(K_all, Y_all)
lrt  <- 2 * (sum(res$ll) - pooled$ll)
df   <- length(SP) - 1
p_lrt <- pchisq(lrt, df, lower.tail = FALSE)

p_pool <- if (is.na(P_POOLED_FIT)) pooled$p else P_POOLED_FIT
res$pred_det_per_rec_pooled <- vapply(SP, function(sp) {
  d <- per_species[[sp]]; rec <- d$Y >= 1
  mean(d$K[rec] * p_pool / (1 - (1 - p_pool)^d$K[rec]))
}, numeric(1))
res$implied_ratio_pooled <- vapply(SP, function(sp) {
  d <- per_species[[sp]]; 1 / mean(1 - (1 - p_pool)^d$K)
}, numeric(1))

res <- merge(res, reported, by = "species", sort = FALSE)
res <- res[order(res$naive), ]

cat("\n", strrep("=", 90), "\n", sep = "")
print(res[, c("species", "naive", "det_freq", "ratio_dn", "p_hat", "se_p",
              "at_bound", "implied_ratio", "ratio_indep", "ratio_hier")],
      row.names = FALSE, digits = 3)

cat(sprintf("\nPooled conditional MLE p = %.3f (se %.3f); used for prediction: %.3f\n",
            pooled$p, pooled$se_p, p_pool))
cat(sprintf("LRT of a common p: %.1f on %d df, p = %.2g (approximate; see header)\n",
            lrt, df, p_lrt))
cat(sprintf("cor(implied ratio, independent-fit ratio)  = %.3f\n",
            cor(res$implied_ratio, res$ratio_indep)))
cat(sprintf("cor(implied ratio, hierarchical-fit ratio) = %.3f\n",
            cor(res$implied_ratio, res$ratio_hier)))

cat("\nPosterior-predictive-style check (plug-in): detections per recorded site-year\n")
print(res[, c("species", "obs_det_per_rec", "pred_det_per_rec_pooled")],
      row.names = FALSE, digits = 3)
cat("A pooled-detection model predicts nearly the same value for every species;\n")
cat("the observed column should be compared against it species by species.\n")

if (any(res$at_bound))
  warning("p_hat at the optimisation bound for: ",
          paste(res$species[res$at_bound], collapse = ", "),
          ". Occupancy is then weakly identified for that species (PaperB, Sec. 2.4).")

write.csv(res, file.path(OUTDIR, "detection_heterogeneity.csv"), row.names = FALSE)
