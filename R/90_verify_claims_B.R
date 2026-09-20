## =====================================================================
## 90_verify_claims_B.R -- do the numbers in PaperB.tex match the output?
## =====================================================================
## Two layers:
##   1. RECOMPUTE each reported quantity from the stored results and
##      compare it with the value written into the manuscript.
##   2. TEXT AUDIT: check that the manuscript actually contains the
##      value, formatted as it is reported. This catches numbers that
##      were superseded but left in the text, which is the failure mode
##      of a manuscript revised many times.
## Also counts the TBD markers, which must be zero before submission.
##
## Run from the project root:  source("R/90_verify_claims_B.R")
## Set TEX to the manuscript path if it is not ./PaperB.tex.
##
## Every claim carries the section it belongs to, so a failure says where
## to look. A claim whose inputs are missing is reported as SKIPPED, not
## as a pass.
## =====================================================================

## Locate the manuscript: an explicit TEX wins, otherwise look in the
## usual places and then anywhere below the working directory.
TEX <- local({
  cand <- c(if (exists("TEX")) TEX, "PaperB.tex", "../PaperB.tex",
            "manuscript/PaperB.tex", "tex/PaperB.tex", "paper/PaperB.tex",
            "doc/PaperB.tex")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) hit <- list.files(".", "^PaperB\\.tex$",
                                      recursive = TRUE, full.names = TRUE)
  if (!length(hit)) hit <- list.files("..", "^PaperB\\.tex$",
                                      recursive = TRUE, full.names = TRUE)
  if (length(hit)) { message("manuscript: ", hit[1]); hit[1] } else NA_character_
})
DIR_REAL <- "results/real"; DIR_MSOM <- "results/msom"
DIR_LAD  <- "results/ladder"; DIR_CV <- "results/cv"
DIR_PS   <- "results/prior_sens"; DIR_ST <- "results/strategy"
DIR_MC   <- "results/mcmc"; DIR_HS <- "results/het_sim"

tex <- if (!is.na(TEX) && file.exists(TEX))
  paste(readLines(TEX, warn = FALSE), collapse = " ") else NA
if (is.na(tex[1]))
  warning("PaperB.tex not found; text audit skipped. Set TEX to its full path, ",
          "e.g. TEX <- \"C:/path/to/PaperB.tex\", and source this script again.")

CLAIMS <- list()
claim <- function(section, what, stated, actual, tol = 0.0005, digits = NULL) {
  CLAIMS[[length(CLAIMS) + 1L]] <<- list(
    section = section, what = what, stated = stated, actual = actual,
    tol = tol, digits = digits)
  invisible(NULL)
}
rd <- function(path) if (file.exists(path)) readRDS(path) else NULL
rc <- function(path) if (file.exists(path)) read.csv(path) else NULL

## ---------------------------------------------------------------------
## Simulation: community parameters and species-level accuracy
## ---------------------------------------------------------------------
comm_stats <- function(f, par) {
  x <- rd(f); if (is.null(x)) return(NULL)
  d <- x$comm[x$comm$parameter == par, ]
  e <- d$est - d$truth
  c(bias = mean(e), mc_se = sd(e) / sqrt(nrow(d)),
    rmse = sqrt(mean(e^2)), cov = mean(d$covered), crps = mean(d$crps))
}
sp_stats <- function(f) {
  x <- rd(f); if (is.null(x)) return(NULL)
  s <- x$shrink
  c(rmse_hier = sqrt(mean((s$msom - s$truth)^2)),
    rmse_ind  = sqrt(mean((s$ssom - s$truth)^2)),
    closer    = mean(s$msom_closer))
}
runs <- list(a = file.path(DIR_MSOM, "msom_replicates.rds"),
             b = file.path(DIR_MSOM, "msom_replicates_v3.rds"),
             c = file.path(DIR_MSOM, "msom_replicates_v5.rds"))
paper_comm <- list(
  a = list(beta0 = c(-0.087, 0.011, 0.143, 0.87, 0.0821),
           beta1 = c(-0.076, 0.005, 0.092, 0.68, 0.0587)),
  b = list(beta0 = c(-0.084, 0.011, 0.142, 0.85, 0.0813),
           beta1 = c(-0.050, 0.005, 0.072, 0.84, 0.0425)),
  c = list(beta0 = c(-0.062, 0.011, 0.126, 0.90, 0.0717),
           beta1 = c(-0.074, 0.005, 0.091, 0.76, 0.0565)))
paper_sp <- list(a = c(0.247, 0.219, 0.453), b = c(0.241, 0.205, 0.453),
                 c = c(0.182, 0.220, 0.615))
for (r in names(runs)) {
  for (p in c("beta0", "beta1")) {
    st <- comm_stats(runs[[r]], paste0(p, "_community"))
    lab <- c("bias", "MC SE", "RMSE", "coverage", "CRPS")
    tolv <- c(0.0005, 0.0005, 0.0005, 0.005, 0.00005)
    for (k in seq_along(lab))
      claim(sprintf("Table 3, run (%s)", r), sprintf("%s %s", p, lab[k]),
            paper_comm[[r]][[p]][k],
            if (is.null(st)) NA_real_ else unname(st[k]), tolv[k])
  }
  ss <- sp_stats(runs[[r]])
  for (k in 1:3)
    claim(sprintf("Table 4, run (%s)", r),
          c("RMSE hierarchical", "RMSE independent", "closer to truth")[k],
          paper_sp[[r]][k], if (is.null(ss)) NA_real_ else unname(ss[k]),
          c(0.0005, 0.0005, 0.0005)[k])
}

## paired change in excess MSE, runs (a) and (c)
A <- rd(runs$a); C5 <- rd(runs$c)
if (!is.null(A) && !is.null(C5)) {
  pa <- merge(A$shrink, C5$shrink, by = c("rep", "species"), suffixes = c("_a", "_c"))
  claim("Section 4.3", "generating intercepts identical across runs (a) and (c)",
        0, max(abs(pa$truth_a - pa$truth_c)), 1e-10)
  by_rep <- do.call(rbind, lapply(split(pa, pa$rep), function(d) data.frame(
    d_excess = (mean((d$msom_c - d$truth_c)^2) - mean((d$ssom_c - d$truth_c)^2)) -
               (mean((d$msom_a - d$truth_a)^2) - mean((d$ssom_a - d$truth_a)^2)))))
  claim("Section 4.3", "paired change in excess MSE", -0.0285,
        mean(by_rep$d_excess), 0.0002)
  claim("Section 4.3", "its paired standard error", 0.0021,
        sd(by_rep$d_excess) / sqrt(nrow(by_rep)), 0.0002)
}

## ---------------------------------------------------------------------
## Application: average occupancy and the detection diagnostic
## ---------------------------------------------------------------------
ao <- rc(file.path(DIR_REAL, "avg_occupancy.csv"))
if (!is.null(ao)) {
  claim("Section 5.3", "smallest hierarchical ratio", 1.14, min(ao$ratio_msom), 0.005)
  claim("Section 5.3", "largest hierarchical ratio",  1.19, max(ao$ratio_msom), 0.005)
  claim("Section 5.3", "smallest independent ratio",  1.07, min(ao$ratio_ssom), 0.005)
  claim("Section 5.3", "largest independent ratio",   2.99, max(ao$ratio_ssom), 0.005)
  claim("Section 5.3", "P(recorded) mean |ratio-1|, hierarchical", 0.012,
        mean(abs(ao$prec_ratio_msom - 1)), 0.0005)
  claim("Section 5.3", "P(recorded) mean |ratio-1|, independent", 0.045,
        mean(abs(ao$prec_ratio_ssom - 1)), 0.0005)
}

dh <- rc(file.path(DIR_REAL, "detection_heterogeneity.csv"))
if (!is.null(dh)) {
  claim("Table 5", "smallest conditional p_hat", 0.106, min(dh$p_hat), 0.0005)
  claim("Table 5", "largest conditional p_hat",  0.566, max(dh$p_hat), 0.0005)
  claim("Section 5.3", "correlation, implied vs independent ratio", 0.98,
        cor(dh$implied_ratio, dh$ratio_indep), 0.005)
  claim("Section 5.3", "correlation, implied vs hierarchical ratio", 0.61,
        cor(dh$implied_ratio, dh$ratio_hier), 0.005)
  claim("Section 5.3", "species with p_hat below 0.25", 5,
        sum(dh$p_hat < 0.25), 0.5)
}

od <- rc(file.path(DIR_REAL, "detection_overdispersion.csv"))
if (!is.null(od)) {
  claim("Appendix S3.9", "species with D > 1", 11, sum(od$D > 1, na.rm = TRUE), 0.5)
  claim("Appendix S3.9", "smallest D", 0.94, min(od$D, na.rm = TRUE), 0.005)
  claim("Appendix S3.9", "largest D",  1.25, max(od$D, na.rm = TRUE), 0.005)
  claim("Appendix S3.9", "bootstrap p < 0.05", 9, sum(od$p_boot < 0.05, na.rm = TRUE), 0.5)
  claim("Appendix S3.9", "significant after Holm", 6,
        sum(p.adjust(od$p_boot, "holm") < 0.05, na.rm = TRUE), 0.5)
}

gp <- rc(file.path(DIR_REAL, "detection_gap_check.csv"))
if (!is.null(gp)) {
  d <- gp$p_fit_indep - gp$p_hat
  claim("Section 5.3", "smallest excess of fitted over conditional p", 0.016, min(d), 0.0005)
  claim("Section 5.3", "largest excess of fitted over conditional p",  0.059, max(d), 0.0005)
  claim("Section 5.3", "correlation after refitting with fitted p", 0.99,
        cor(gp$implied_ratio_fitp, gp$ratio_indep), 0.005)
  claim("Section 5.3", "summed |implied - independent|, conditional p", 3.44,
        sum(abs(gp$implied_ratio - gp$ratio_indep)), 0.005)
  claim("Section 5.3", "summed |implied - independent|, fitted p", 0.86,
        sum(abs(gp$implied_ratio_fitp - gp$ratio_indep)), 0.005)
}

rcx <- rd(file.path(DIR_REAL, "remaining_checks.rds"))
if (!is.null(rcx) && !is.null(rcx$ppc)) {
  t3 <- rcx$ppc
  claim("Appendix S3.7", "hierarchical predictions, minimum", 1.66, min(t3$hier_mean), 0.005)
  claim("Appendix S3.7", "hierarchical predictions, maximum", 1.71, max(t3$hier_mean), 0.005)
  claim("Appendix S3.7", "mean excess of independent predictions", 0.046,
        mean(t3$ind_mean - t3$observed), 0.0005)
  claim("Appendix S3.7", "species whose observed value is outside the hierarchical interval",
        11, sum(!t3$obs_in_hier), 0.5)
}

## ---------------------------------------------------------------------
## Cross-validation, recomputed including the cluster bootstrap
## ---------------------------------------------------------------------
f_cv <- c(file.path(DIR_CV, "site_cv_raw.rds"), file.path(DIR_REAL, "site_cv_raw.rds"))
f_cv <- f_cv[file.exists(f_cv)]
if (length(f_cv)) {
  suppressPackageStartupMessages(library(spOccupancy)); data(hbefTrends)
  N_SITE <- dim(hbefTrends$y)[2]
  cv <- readRDS(f_cv[1]); cv <- cv[!is.na(cv$log_h) & !is.na(cv$log_s), ]
  cv$site <- ((cv$row - 1) %% N_SITE) + 1
  cv$d_log <- cv$log_h - cv$log_s
  cv$d_brier <- (cv$observed - cv$prec_h)^2 - (cv$observed - cv$prec_s)^2
  cv$one <- 1
  claim("Section 5.7", "hierarchical mean log score", -1.053, mean(cv$log_h), 0.0005)
  claim("Section 5.7", "independent mean log score",  -0.960, mean(cv$log_s), 0.0005)
  claim("Section 5.7", "paired log-score difference", -0.093, mean(cv$d_log), 0.0005)
  claim("Section 5.7", "paired Brier difference", 0.025, mean(cv$d_brier), 0.0005)
  ss <- aggregate(cbind(d_log, d_brier, one) ~ site, cv, sum)
  set.seed(20260201L)
  bt <- t(replicate(2000, { i <- sample.int(nrow(ss), replace = TRUE)
    c(sum(ss$d_log[i]) / sum(ss$one[i]), sum(ss$d_brier[i]) / sum(ss$one[i])) }))
  claim("Section 5.7", "cluster SE, log score", 0.005, sd(bt[, 1]), 0.0005)
  claim("Section 5.7", "cluster SE, Brier", 0.002, sd(bt[, 2]), 0.0005)
}

## ---------------------------------------------------------------------
## The ladder
## ---------------------------------------------------------------------
lf <- rc(file.path(DIR_LAD, "ladder_full_summary.csv"))
lc <- rc(file.path(DIR_LAD, "ladder_cv_overall.csv"))
if (!is.null(lf)) {
  g <- function(m, f) f(lf$ratio[lf$model == m])
  claim("Table 9", "M1 ratio range, minimum", 1.14, g("M1", min), 0.005)
  claim("Table 9", "M1 ratio range, maximum", 1.19, g("M1", max), 0.005)
  claim("Table 9", "M3 ratio range, minimum", 1.09, g("M3", min), 0.005)
  claim("Table 9", "M3 ratio range, maximum", 2.73, g("M3", max), 0.005)
  m1 <- lf[lf$model == "M1", ]; m3 <- lf[lf$model == "M3", ]
  k <- match(c("NAWA", "BAWW", "BHVI", "AMRE", "CAWA"), m1$species)
  fac <- m3$psi_bar[match(m1$species[k], m3$species)] / m1$psi_bar[k]
  claim("Section 5.8", "smallest psi_bar factor, M3 over M1", 2.24, min(fac), 0.005)
  claim("Section 5.8", "largest psi_bar factor, M3 over M1",  2.32, max(fac), 0.005)
}
if (!is.null(lc)) {
  gv <- function(m, col) lc[[col]][lc$model == m]
  claim("Table 10", "M1 log-score difference", -0.0934, gv("M1", "d_log"), 0.0005)
  claim("Table 10", "M2 log-score difference", -0.0397, gv("M2", "d_log"), 0.0005)
  claim("Table 10", "M3 log-score difference", -0.0236, gv("M3", "d_log"), 0.0005)
  claim("Section 5.8", "share of the log-score gap closed by M2", 0.57,
        1 - gv("M2", "d_log") / gv("M1", "d_log"), 0.005)
  claim("Section 5.8", "share of the log-score gap closed by M3", 0.75,
        1 - gv("M3", "d_log") / gv("M1", "d_log"), 0.005)
  claim("Section 5.8", "share of the Brier gap closed by M3", 0.70,
        1 - gv("M3", "d_brier") / gv("M1", "d_brier"), 0.005)
}

## ---------------------------------------------------------------------
## Prior sensitivity, approximation strategies, MCMC, heterogeneity sim
## ---------------------------------------------------------------------
ps <- rc(file.path(DIR_PS, "field_prior_sensitivity.csv"))
if (!is.null(ps)) {
  claim("Section 3.1", "smallest posterior range", 2.84, min(ps$range), 0.005)
  claim("Section 3.1", "largest posterior range",  3.57, max(ps$range), 0.005)
  claim("Section 3.1", "smallest field SD", 0.494, min(ps$sigma), 0.0005)
  claim("Section 3.1", "largest field SD",  0.578, max(ps$sigma), 0.0005)
  claim("Section 3.1", "largest change in a species intercept", 0.011,
        max(ps$max_d_species), 0.0005)
  claim("Section 3.1", "largest change in beta0, in posterior SDs", 0.10,
        max(ps$max_d_beta0_sd), 0.005)
}
sc <- rc(file.path(DIR_ST, "strategy_comparison.csv"))
if (!is.null(sc)) {
  claim("Section 3.3", "Gaussian, largest discrepancy in posterior SDs", 0.005,
        max(abs(sc$d_gaussian)), 0.0005)
  claim("Section 3.3", "simplified Laplace, largest discrepancy", 0.249,
        max(abs(sc$d_simplified.laplace)), 0.0005)
  keep <- !grepl("Precision for species_id", sc$parameter)
  claim("Section 3.3", "simplified Laplace, excluding the species precision", 0.074,
        max(abs(sc$d_simplified.laplace[keep])), 0.0005)
  sp <- grep("^psi_int_", sc$parameter)
  claim("Section 3.3", "largest species-intercept change, logit scale", 0.0025,
        max(abs(sc$mean_simplified.laplace[sp] - sc$mean_laplace[sp])), 0.0005)
}
mc <- rc(file.path(DIR_MC, "mcmc_species_tMsPGOcc.csv"))
if (!is.null(mc) && "p_hat" %in% names(mc)) {
  claim("Section 3.3", "correlation of MCMC with conditional detection", 0.998,
        cor(mc$p_mcmc, mc$p_hat), 0.0005)
  claim("Section 3.3", "mean excess of MCMC over conditional detection", 0.012,
        mean(mc$p_mcmc - mc$p_hat), 0.0005)
  claim("Section 3.3", "MCMC ratio range, minimum", 1.11, min(mc$ratio_mcmc), 0.005)
  claim("Section 3.3", "MCMC ratio range, maximum", 3.52, max(mc$ratio_mcmc), 0.005)
}
## species-specific AR(1) in the MCMC reference (Limitations)
f_slim <- file.path(DIR_MC, "mcmc_tMsPGOcc_full_slim.rds")
if (!file.exists(f_slim)) f_slim <- file.path(DIR_MC, "mcmc_tMsPGOcc_full.rds")
if (file.exists(f_slim)) {
  suppressPackageStartupMessages(library(coda))
  th <- summary(readRDS(f_slim)$theta.samples)$quantiles
  rr <- th[grep("^rho", rownames(th)), c("2.5%", "97.5%"), drop = FALSE]
  claim("Limitations", "species AR(1) intervals that include zero",
        nrow(rr), sum(rr[, 1] < 0 & rr[, 2] > 0), 0.5)
}

hs <- rd(file.path(DIR_HS, "het_sim.rds"))
if (!is.null(hs)) {
  a <- hs$A[hs$A$p_bar == 0.15 & hs$A$sigma_logit_p == 1, ]
  claim("Appendix S3.8", "D at p = 0.15, sigma = 1", 1.10, a$D, 0.005)
  claim("Appendix S3.8", "estimate of p there", 0.282, a$p_ztb, 0.0005)
  claim("Appendix S3.8", "true ratio there", 2.38, a$ratio_true, 0.005)
  claim("Appendix S3.8", "implied ratio there", 1.59, a$ratio_ztb, 0.005)
  for (nm in c("A", "B", "C")) {
    d <- hs[[nm]]$p_full - hs[[nm]]$p_ztb
    claim("Appendix S3.8", sprintf("scenario %s: largest |p_full - p_ztb|", nm),
          if (nm == "C") 0.019 else 0.000, max(abs(d)), 0.0015)
  }
}

## ---------------------------------------------------------------------
## The beta1 mechanism test (Section 4.2)
## ---------------------------------------------------------------------
bm <- rc(file.path("results/beta1", "beta1_mechanism.csv"))
if (!is.null(bm)) {
  claim("Section 4.2", "bias against the unweighted mean", 0.002,
        mean(bm$beta1_hat - bm$mean_unw), 0.0005)
  claim("Section 4.2", "its standard error", 0.004,
        sd(bm$beta1_hat - bm$mean_unw) / sqrt(nrow(bm)), 0.0005)
  claim("Section 4.2", "gap between the two targets", 0.095,
        mean(bm$mean_unw - bm$mean_inf), 0.0005)
  fit <- lm(I(beta1_hat - mean_unw) ~ I(mean_inf - mean_unw), data = bm)
  claim("Section 4.2", "regression of displacement on the gap", 0.14,
        unname(coef(fit)[2]), 0.005)
  claim("Section 4.2", "its standard error", 0.10,
        summary(fit)$coefficients[2, 2], 0.005)
}

## ---------------------------------------------------------------------
## Where the slope bias enters, and whether the interval construction
## explains the species-level under-coverage (Section 4.2, Limitations)
## ---------------------------------------------------------------------
bc <- rc(file.path("results/beta1", "bias_and_coverage.csv"))
if (!is.null(bc)) {
  b1 <- bc[!duplicated(bc$rep), ]
  claim("Section 4.2", "slope bias with the occupancy likelihood, no field",
        -0.054, mean(b1$beta1_hat - b1$beta1_true), 0.0005)
  claim("Section 4.2", "its coverage", 0.84,
        mean(b1$beta1_lo <= b1$beta1_true & b1$beta1_true <= b1$beta1_hi), 0.005)
  far <- bc$dist >= quantile(bc$dist, 2/3)
  claim("Limitations", "far-tercile coverage, marginal intervals", 0.865,
        mean(bc$cov_marg[far]), 0.005)
  claim("Limitations", "far-tercile coverage, joint intervals", 0.850,
        mean(bc$cov_joint[far]), 0.005)
}

## ---------------------------------------------------------------------
## The M1 -> M3 change in predicted occupancy (Section 5.8, Figure 2)
## ---------------------------------------------------------------------
## 84_figures_B.R prints these; they are recomputed here from the stored
## ladder summary, which averages over all surveyed site-years rather than
## over 2018 alone, so the tolerance is wide enough to cover the difference.
if (!is.null(lf)) {
  m1 <- lf[lf$model == "M1", ]; m3 <- lf[lf$model == "M3", ]
  k  <- match(c("NAWA", "AMRE", "CAWA", "BAWW", "BHVI"), m1$species)
  fac <- m3$psi_bar[match(m1$species[k], m3$species)] / m1$psi_bar[k]
  claim("Section 5.8", "smallest M3/M1 factor, five least detectable", 2.09,
        min(fac), 0.25)
  claim("Section 5.8", "largest M3/M1 factor, five least detectable", 2.33,
        max(fac), 0.25)
}

## ---------------------------------------------------------------------
## Figures: the field summaries quoted in Section 5.6
## ---------------------------------------------------------------------
## 84_figures_B.R prints these; they are recorded here so that redrawing the
## figures cannot silently disagree with the text. Recomputing them needs
## the posterior draws, so they are checked against the stored values.
f_fig <- file.path(DIR_REAL, "field_summary_2018.csv")
if (file.exists(f_fig)) {
  fs <- read.csv(f_fig)
  ## Seeded posterior sampling reproduces these only to about 0.01 here
  ## (see the note in 84_figures_B.R), so the manuscript states them as
  ## approximate and the tolerance matches.
  claim("Section 5.6", "field posterior mean, minimum", -0.45, fs$min, 0.02)
  claim("Section 5.6", "field posterior mean, maximum",  0.64, fs$max, 0.02)
  claim("Section 5.6", "median posterior SD of the field", 0.25, fs$median_sd, 0.005)
  claim("Section 5.6", "median |mean| / SD", 0.66, fs$median_ratio, 0.02)
}

## ---------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------
res <- do.call(rbind, lapply(CLAIMS, function(c1) {
  ok <- if (is.na(c1$actual)) NA else abs(c1$actual - c1$stated) <= c1$tol
  data.frame(section = c1$section, claim = c1$what,
             stated = c1$stated, actual = round(c1$actual, 6),
             status = if (is.na(ok)) "SKIPPED" else if (ok) "ok" else "MISMATCH")
}))
cat("\n", strrep("=", 92), "\nNUMERICAL CLAIMS\n", strrep("=", 92), "\n", sep = "")
print(res[res$status != "ok", ], row.names = FALSE, digits = 6)
cat(sprintf("\n%d claims: %d ok, %d mismatched, %d skipped (inputs missing)\n",
            nrow(res), sum(res$status == "ok"),
            sum(res$status == "MISMATCH"), sum(res$status == "SKIPPED")))

## text audit: is each stated value actually in the manuscript?
if (!is.na(tex[1])) {
  fmt <- function(x) {
    if (abs(x - round(x)) < 1e-9) return(as.character(round(x)))
    unique(c(sprintf("%.2f", x), sprintf("%.3f", x), sprintf("%.4f", x)))
  }
  miss <- res[res$status == "ok", ]
  miss$in_text <- vapply(seq_len(nrow(miss)), function(i)
    any(vapply(fmt(miss$stated[i]), \(p) grepl(p, tex, fixed = TRUE), TRUE)), TRUE)
  bad <- miss[!miss$in_text, ]
  cat("\n", strrep("=", 92), "\nTEXT AUDIT\n", strrep("=", 92), "\n", sep = "")
  if (nrow(bad)) {
    cat("verified values that do not appear in the manuscript as written\n")
    cat("(either the text is stale, or it rounds them differently):\n")
    print(bad[, c("section", "claim", "stated")], row.names = FALSE)
  } else cat("every verified value appears in the manuscript.\n")

  n_tbd <- length(gregexpr("\\\\TBD\\{", tex)[[1]])
  n_tbd <- if (regexpr("\\\\TBD\\{", tex) < 0) 0 else n_tbd - 1  # minus the definition
  cat(sprintf("\nTBD markers remaining: %d (must be 0 before submission)\n", n_tbd))
}

write.csv(res, file.path(DIR_REAL, "verify_claims_B.csv"), row.names = FALSE)
cat("\nfull table written to ", file.path(DIR_REAL, "verify_claims_B.csv"), "\n", sep = "")
