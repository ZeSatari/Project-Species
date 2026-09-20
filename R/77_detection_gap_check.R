## =====================================================================
## 77_detection_gap_check.R -- why do implied ratios exceed ratio_indep?
## =====================================================================
## PaperB, Section 5.3, first caveat.
##
## The conditional estimate p_hat_i (script 68) implies
##     psibar/n = 1 / mean_st{ 1 - (1 - p_hat_i)^K_st },
## which exceeds the independent fits' psibar/n for every species.
## Hypothesis: each independent fit's detection probability, averaged over
## the visits actually made, is higher than p_hat_i. If so, recomputing the
## implied ratio with the fitted detection probabilities should reproduce
## ratio_indep.
##
## Self-contained: loads fits.rds and rebuilds only the helpers it needs,
## exactly as 63_avg_occupancy.R does. det_draws() comes from the same
## sourced files that 63 uses. No posterior draws of psi are taken.
##
## Requires: results/real/fits.rds (from the application pipeline) and
##           results/real/detection_heterogeneity.csv (from script 68).
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

if (!exists("det_draws", mode = "function"))
  stop("det_draws() not found after sourcing R/00_setup.R, R/10_metrics.R ",
       "and R/30_fit_occupancy.R. Locate its definition with\n",
       "  grep -n 'det_draws <-' R/*.R\n",
       "and source that file above.")

OUTDIR  <- "results/real"
N_DRAWS <- 1000
set.seed(20260201L)   # in case det_draws() uses R's RNG; see note at end

fits <- readRDS(file.path(OUTDIR, "fits.rds"))
ssom <- fits$ssom
dh   <- read.csv(file.path(OUTDIR, "detection_heterogeneity.csv"))

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K        <- dim(hbefTrends$y)[4]
N_ROW    <- N_SITE * N_YEAR

## --- helpers, identical to 63_avg_occupancy.R -------------------------
stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
}
p_detect_mean <- function(fit, Xday, Xtod) {
  dr <- det_draws(fit, n = N_DRAWS)          # columns: intercept, day, tod
  stopifnot(ncol(dr) >= 3)
  pm <- matrix(NA_real_, N_ROW, K)
  for (j in seq_len(K)) {
    eta <- outer(rep(1, N_ROW), dr[, 1]) +
           outer(Xday[, j], dr[, 2]) +
           outer(Xtod[, j], dr[, 3])
    pm[, j] <- rowMeans(plogis(eta))
  }
  pm
}

Xday <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")
Xday[is.na(Xday)] <- 0;       Xtod[is.na(Xtod)] <- 0

## --- the check --------------------------------------------------------
gap <- do.call(rbind, lapply(SP_NAMES, function(sp) {
  message("species ", sp)
  Y    <- stack_species(sp)                 # N_ROW x K, NA = visit not made
  made <- !is.na(Y)
  surv <- rowSums(made) > 0

  p_s <- p_detect_mean(ssom[[sp]], Xday, Xtod)
  p_s[!made] <- 0                           # missing visits contribute nothing

  P_rec <- 1 - apply(1 - p_s, 1, prod)      # P(recorded | occupied), visits made
  data.frame(
    species            = sp,
    p_fit_indep        = mean(p_s[made]),
    implied_ratio_fitp = 1 / mean(P_rec[surv])
  )
}))

gap <- merge(gap, dh[, c("species", "naive", "p_hat", "se_p",
                         "implied_ratio", "ratio_indep")],
             by = "species", sort = FALSE)
gap <- gap[order(gap$naive), ]
gap$z_p <- (gap$p_fit_indep - gap$p_hat) / gap$se_p

cat("\n", strrep("=", 90), "\n", sep = "")
print(gap[, c("species", "p_hat", "p_fit_indep", "z_p",
              "implied_ratio", "implied_ratio_fitp", "ratio_indep")],
      row.names = FALSE, digits = 3)

cat("\nIf implied_ratio_fitp ~ ratio_indep while implied_ratio > ratio_indep,\n")
cat("the gap is explained by the independent fits estimating higher detection\n")
cat("than the conditional likelihood.\n")
cat("Note: if det_draws() uses inla.hyperpar.sample(), set.seed() above may not\n")
cat("control it; rerun once to confirm the table is stable to 2 decimals.\n")

write.csv(gap, file.path(OUTDIR, "detection_gap_check.csv"), row.names = FALSE)
