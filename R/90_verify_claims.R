## =====================================================================
## 90_verify_claims.R -- check every numerical claim in the manuscript
## =====================================================================
## Each claim asserted in the text is restated here as a testable
## expression evaluated against the saved results. Run before submission,
## and again after any rerun: if a number changes, this reports which
## sentence is now wrong instead of leaving it to a reviewer.
##
## Adding a claim to the manuscript means adding a line here.
##
## Exit behaviour: prints a table of PASS/FAIL and stops with an error if
## anything fails, so it can be wired into a build.
## =====================================================================

suppressPackageStartupMessages({ library(dplyr) })

CLAIMS <- list()
claim <- function(id, where, text, value, expected, tol = 0.0006) {
  CLAIMS[[length(CLAIMS) + 1]] <<- list(
    id = id, where = where, text = text,
    value = tryCatch(as.numeric(value), error = function(e) NA_real_),
    expected = expected, tol = tol
  )
  invisible(NULL)
}

## ---------------------------------------------------------------------
## Load results
## ---------------------------------------------------------------------
rec   <- read.csv("results/single/TableS2_recovery.csv")
## The confounded regime was run separately; its recovery table is
## written by aggregate_reps() on the M2 rows of confounded_params.rds.
f_conf <- "results/single/confounded_recovery.csv"
rec_conf <- if (file.exists(f_conf)) read.csv(f_conf) else NULL
cmp   <- read.csv("results/single/TableS1_comparison.csv")
shr   <- read.csv("results/msom/shrinkage.csv")
sp    <- read.csv("results/real/Table1_species.csv")
det   <- read.csv("results/real/detection_coefficients_M2.csv")
ssc   <- read.csv("results/real/single_species_comparison.csv")
conf  <- read.csv("results/real/spatial_confounding.csv")
msom  <- readRDS("results/msom/msom_replicates.rds")

g  <- function(d, ...) { f <- list(...); r <- d
         for (k in names(f)) r <- r[r[[k]] == f[[k]], , drop = FALSE]
         stopifnot(nrow(r) == 1); r }
inv <- function(x) 1 / (1 + exp(-x))

## =====================================================================
## Results, single-species simulation
## =====================================================================
claim("R1.1", "Results 3.1.1",
      "baseline alpha0 bias is 0.10 to 0.12 across regimes",
      max(abs(rec$bias[rec$model == "M1" & rec$parameter == "alpha0"])),
      expected = 0.117, tol = 0.001)

for (s in c("independent", "shared", "observer")) {
  claim(paste0("R1.2.", s), "Results 3.1.1 / tab:s2-recovery",
        paste("M1 alpha0 coverage,", s),
        g(rec, scenario = s, model = "M1", parameter = "alpha0")$coverage,
        expected = c(independent = 0.655, shared = 0.550, observer = 0.595)[[s]])
  claim(paste0("R1.3.", s), "Results 3.1.1 / tab:s2-recovery",
        paste("M2 alpha0 coverage,", s),
        g(rec, scenario = s, model = "M2", parameter = "alpha0")$coverage,
        expected = c(independent = 0.940, shared = 0.930, observer = 0.855)[[s]])
}

## The text originally asserted 0.005, taken from a single replicate and
## never checked against the full run. The true maximum is 0.0117, driven
## by the observer regime. The text now states 0.012.
claim("R1.4", "Results 3.1.1",
      "largest M1-M2 difference in occupancy bias is 0.012",
      max(abs(
        sapply(c("beta0", "beta1"), function(p)
          sapply(c("independent", "shared", "observer"), function(s)
            g(rec, scenario = s, model = "M1", parameter = p)$bias -
            g(rec, scenario = s, model = "M2", parameter = p)$bias))
      )),
      expected = 0.012, tol = 0.001)

## The abstract, the highlights and the Discussion all quote a rounded
## coverage range rather than the three separate figures. Lock the
## endpoints, so that a rerun which moves either one flags the rounded
## range as well. An earlier draft quoted 0.86-0.95, whose upper end
## coincided with the nominal level and so read as though coverage had
## been fully restored.
a0_m2 <- rec$coverage[rec$model == "M2" & rec$parameter == "alpha0"]
a0_m1 <- rec$coverage[rec$model == "M1" & rec$parameter == "alpha0"]

claim("R1.7", "Abstract / Highlights / Results / Discussion",
      "M2 alpha0 coverage: maximum is 0.940, quoted as 0.94",
      max(a0_m2), expected = 0.940)
claim("R1.8", "Abstract / Highlights / Results / Discussion",
      "M2 alpha0 coverage: minimum is 0.855, quoted as 0.86",
      min(a0_m2), expected = 0.855)
claim("R1.9", "Abstract / Highlights / Discussion",
      "M1 alpha0 coverage: minimum is 0.550, quoted as 0.55",
      min(a0_m1), expected = 0.550)
claim("R1.10", "Abstract / Highlights / Discussion",
      "M1 alpha0 coverage: maximum is 0.655, quoted as 0.66",
      max(a0_m1), expected = 0.655)
claim("R1.11", "Results 3.1.1 / Discussion",
      "M2 alpha0 bias range 0.016 to 0.060",
      max(rec$bias[rec$model == "M2" & rec$parameter == "alpha0"]),
      expected = 0.060, tol = 0.001)

## Detection-parameter qualifications now stated in Limitations and
## Conclusions, replacing a bare claim that they were "unaffected".
claim("R1.12", "Limitations / Conclusions",
      "observer regime: M2 alpha0 bias 0.060",
      g(rec, scenario = "observer", model = "M2", parameter = "alpha0")$bias,
      expected = 0.060, tol = 0.001)
claim("R1.13", "Limitations / Conclusions",
      "observer regime: M2 alpha0 coverage 0.855",
      g(rec, scenario = "observer", model = "M2", parameter = "alpha0")$coverage,
      expected = 0.855)

## The confounded-regime figures are two distinct quantities, and the
## text now separates them: alpha_2's bias under `confounded` (0.053, or
## 7.5% of the generating value) and its paired displacement from
## `independent` (0.062). They reconcile because alpha_2 carried a bias
## of -0.010 under `independent`. All three are locked, so the
## arithmetic cannot drift apart in a rerun.
if (!is.null(rec_conf)) {
  gc2 <- function(par) {
    r <- rec_conf[rec_conf$parameter == par, ]
    stopifnot(nrow(r) == 1); r
  }
  claim("R1.14", "Results 3.1.1 / Discussion / Abstract",
        "confounded regime: alpha2 bias is 7.5% of the generating value",
        100 * gc2("alpha2")$bias / 0.700, expected = 7.5, tol = 0.2)
  claim("R1.15", "Results 3.1.1",
        "independent regime: alpha2 bias -0.010",
        g(rec, scenario = "independent", model = "M2",
          parameter = "alpha2")$bias, expected = -0.010, tol = 0.002)
  claim("R1.16", "Results 3.1.1",
        "the two reconcile: 0.053 minus (-0.010) is the 0.062 displacement",
        gc2("alpha2")$bias -
          g(rec, scenario = "independent", model = "M2",
            parameter = "alpha2")$bias,
        expected = 0.062, tol = 0.002)
  claim("R1.17", "Results 3.1.1 / tab:s2-confounded",
        "confounded regime: beta0 bias 0.049",
        gc2("beta0")$bias, expected = 0.049, tol = 0.001)
  claim("R1.18", "Results 3.1.1 / tab:s2-confounded",
        "confounded regime: alpha0 coverage 0.950",
        gc2("alpha0")$coverage, expected = 0.950)
  claim("R1.19", "Results 3.1.1 / tab:s2-confounded",
        "confounded regime: beta1 coverage 0.895",
        gc2("beta1")$coverage, expected = 0.895)
} else {
  message("confounded_recovery.csv not found; confounded-regime claims skipped.")
}

## The shared regime is now described by numbers rather than as
## "unaffected", so those numbers are locked.
claim("R1.20", "Discussion 4.2",
      "shared regime: M2 beta0 bias 0.011",
      g(rec, scenario = "shared", model = "M2", parameter = "beta0")$bias,
      expected = 0.011, tol = 0.001)
claim("R1.21", "Discussion 4.2",
      "shared regime: M2 beta0 coverage 0.870",
      g(rec, scenario = "shared", model = "M2", parameter = "beta0")$coverage,
      expected = 0.870)
claim("R1.22", "Discussion 4.2",
      "independent regime: M2 beta0 bias 0.022",
      g(rec, scenario = "independent", model = "M2", parameter = "beta0")$bias,
      expected = 0.022, tol = 0.001)
claim("R1.23", "Discussion 4.2",
      "independent regime: M2 beta0 coverage 0.905",
      g(rec, scenario = "independent", model = "M2", parameter = "beta0")$coverage,
      expected = 0.905)

claim("R1.5", "Results 3.1.1",
      "M2 achieved lower WAIC in every replicate",
      min(cmp$prop_M2_better), expected = 1.0)

claim("R1.6", "Results 3.1.1",
      "ULOOCV gains 0.0058 / 0.0062 / 0.0024",
      g(cmp, scenario = "independent")$d_ulogcv_mean,
      expected = 0.0058, tol = 0.0001)

## ---------------------------------------------------------------------
## Field-scaling diagnostic (71_bias_field_scaling.R). These three back
## the claim in Appendix S1.1 that the slope bias is not attenuation
## induced by the latent field. Skipped with a message rather than an
## error if the diagnostic has not been run, since it is not needed to
## reproduce any table.
## ---------------------------------------------------------------------
f_bs <- "results/bias_scaling/bias_scaling.csv"
if (file.exists(f_bs)) {
  bs <- read.csv(f_bs)
  ctrl <- bs[which.min(bs$sd_field), ]

  claim("R1.24", "Appendix S1.1",
        "bias at sigma ~ 0 is -0.018, so the field is not the mechanism",
        ctrl$bias, expected = -0.018, tol = 0.002)

  claim("R1.25", "Appendix S1.1",
        "Monte Carlo SE of the control bias is 0.008",
        ctrl$mc_se, expected = 0.008, tol = 0.001)

  ## Weighted least squares of bias on the generating sigma, weights
  ## 1/mc_se^2. The manuscript reports the slope and its standard error.
  w  <- 1 / bs$mc_se^2
  mx <- sum(w * bs$sd_field) / sum(w)
  my <- sum(w * bs$bias)     / sum(w)
  sxx <- sum(w * (bs$sd_field - mx)^2)
  slope <- sum(w * (bs$sd_field - mx) * (bs$bias - my)) / sxx

  claim("R1.26", "Appendix S1.1",
        "weighted slope of bias on sigma is -0.019",
        slope, expected = -0.019, tol = 0.002)

  claim("R1.27", "Appendix S1.1",
        "slope SE is 0.012, so the trend is not resolved at R = 40",
        sqrt(1 / sxx), expected = 0.012, tol = 0.002)

  claim("R1.28", "Appendix S1.1",
        "coverage stays within 0.925-0.975 across field variances",
        as.numeric(min(bs$coverage) >= 0.925 && max(bs$coverage) <= 0.975),
        expected = 1)
} else {
  message("skipping R1.24-R1.28: ", f_bs, " not found ",
          "(run R/71_bias_field_scaling.R)")
}

## ---------------------------------------------------------------------
## MCMC benchmark (72_mcmc_benchmark.R). These back the Limitations
## paragraph reporting that the detection component agrees with an MCMC
## fit while the occupancy covariates do not, and Appendix S4.
## ---------------------------------------------------------------------
f_mc <- "results/mcmc_benchmark/mcmc_estimates.csv"
if (file.exists(f_mc)) {
  mc <- read.csv(f_mc, stringsAsFactors = FALSE)
  gm <- function(m, p) mc$mean[mc$model == m & mc$parameter == p]

  claim("R9.1", "Limitations / Appendix S4",
        "MCMC biotic coefficient 0.403",
        gm("M2", "oven"), expected = 0.403, tol = 0.002)

  claim("R9.2", "Limitations / Appendix S4",
        "MCMC linear elevation -1.46 under M2",
        gm("M2", "elev_s"), expected = -1.458, tol = 0.005)

  claim("R9.3", "Limitations / Appendix S4",
        "MCMC linear elevation -1.48 under M1",
        gm("M1", "elev_s"), expected = -1.476, tol = 0.005)

  ## The comparison the text makes is in units of the INLA posterior SD,
  ## so it is recomputed here rather than taken as a stated number.
  f_in <- "results/real/m1_m2_summaries.rds"
  if (file.exists(f_in)) {
    ins <- readRDS(f_in)
    ib  <- ins$m2$hyper[grep("^beta\\[", rownames(ins$m2$hyper)), ]
    a2_mean <- ib[4, "mean"]; a2_sd <- ib[4, "sd"]

    claim("R9.4", "Limitations / Appendix S4",
          "biotic coefficient agrees to within 0.26 posterior SD",
          abs(a2_mean - gm("M2", "oven")) / a2_sd,
          expected = 0.26, tol = 0.04)

    el <- ins$m2$fixed["elev_s", ]
    claim("R9.5", "Limitations / Appendix S4",
          "linear elevation differs by nearly three posterior SD",
          abs(el[["mean"]] - gm("M2", "elev_s")) / el[["sd"]],
          expected = 2.95, tol = 0.10)
  } else {
    message("skipping R9.4-R9.5: ", f_in, " not found (run 82_slim_outputs.R)")
  }
} else {
  message("skipping R9.1-R9.5: ", f_mc, " not found ",
          "(run R/72_mcmc_benchmark.R)")
}

## ---------------------------------------------------------------------
## Prior sensitivity for the field (73_prior_sensitivity_field.R).
## ---------------------------------------------------------------------
f_pf <- "results/prior_field/prior_field.csv"
if (file.exists(f_pf)) {
  pf <- read.csv(f_pf, stringsAsFactors = FALSE)
  m2 <- pf[pf$model == "M2" & !is.na(pf$alpha2), ]
  sd_ref <- m2$alpha2_sd[m2$prior == "baseline"]

  claim("R9.6", "Limitations",
        "alpha_2 spans 0.4168 to 0.4179 across six prior specifications",
        min(m2$alpha2), expected = 0.4168, tol = 0.0005)

  claim("R9.7", "Limitations",
        "that spread is 0.02 of one posterior SD",
        diff(range(m2$alpha2)) / sd_ref, expected = 0.02, tol = 0.01)

  claim("R9.8", "Limitations",
        "practical range stays between 2.96 and 3.42 km",
        as.numeric(min(m2$range) >= 2.95 && max(m2$range) <= 3.43),
        expected = 1)

  claim("R9.9", "Limitations",
        "marginal standard deviation stays between 2.64 and 3.17",
        as.numeric(min(m2$sigma) >= 2.63 && max(m2$sigma) <= 3.18),
        expected = 1)

  claim("R9.10", "Limitations",
        "linear elevation moves 0.39 of a posterior SD",
        diff(range(m2$elev)) / 0.3112, expected = 0.39, tol = 0.03)
} else {
  message("skipping R9.6-R9.10: ", f_pf, " not found ",
          "(run R/73_prior_sensitivity_field.R)")
}

## =====================================================================
## Results, empirical single-species
## =====================================================================
b <- det[grepl("biotic", det$coefficient), ]
claim("R2.1", "Results 3.1.2 / tab:det-emp",
      "biotic coefficient 0.419", b$mean, expected = 0.419, tol = 0.001)
claim("R2.2", "Results 3.1.2", "biotic 95% CrI lower 0.310",
      b$lwr, expected = 0.310, tol = 0.001)
claim("R2.3", "Results 3.1.2", "biotic 95% CrI upper 0.523",
      b$upr, expected = 0.523, tol = 0.001)

claim("R2.4", "Results 3.1.2",
      "detection probability rises from 0.505 to 0.608",
      inv(det$mean[det$coefficient == "intercept"]),
      expected = 0.505, tol = 0.002)
claim("R2.5", "Results 3.1.2",
      "detection probability with biotic detection",
      inv(det$mean[det$coefficient == "intercept"] + b$mean),
      expected = 0.608, tol = 0.002)

claim("R2.6", "Results 3.1.2", "day and tod credible intervals span zero",
      as.numeric(all(det$lwr[det$coefficient %in% c("day", "tod")] < 0 &
                     det$upr[det$coefficient %in% c("day", "tod")] > 0)),
      expected = 1)

claim("R2.7", "Results 3.1.2", "WAIC M1 9936.5",
      ssc$waic[ssc$model == "M1"], expected = 9936.5, tol = 0.5)
claim("R2.8", "Results 3.1.2", "WAIC M2 9886.0",
      ssc$waic[ssc$model == "M2"], expected = 9886.0, tol = 0.5)
claim("R2.9", "Results 3.1.2", "ULOOCV improves by 0.008",
      ssc$ulogcv[ssc$model == "M2"] - ssc$ulogcv[ssc$model == "M1"],
      expected = 0.008, tol = 0.001)

## Spatial confounding
cl <- conf[conf$term == "elev_s", ]; cq <- conf[conf$term == "elev_s2", ]
claim("R2.10", "Results 3.1.2", "elevation shift 1.7 posterior SD",
      cl$shift_in_sd, expected = 1.7, tol = 0.05)
claim("R2.11", "Results 3.1.2", "elevation^2 shift 3.2 posterior SD",
      cq$shift_in_sd, expected = 3.2, tol = 0.05)
claim("R2.12", "Results 3.1.2", "SD inflation 5.9 (linear)",
      cl$sd_inflation, expected = 5.9, tol = 0.1)
claim("R2.13", "Results 3.1.2", "SD inflation 3.8 (quadratic)",
      cq$sd_inflation, expected = 3.8, tol = 0.1)

## =====================================================================
## Results, MSOM simulation
## =====================================================================
cm <- msom$comm
b0 <- cm[cm$parameter == "beta0_community", ]
b1 <- cm[cm$parameter == "beta1_community", ]
##[PaperB] claim("R3.1", "Results 3.2.1 / tab:s3-community", "community beta0 bias -0.087",
##[PaperB]       mean(b0$est) - b0$truth[1], expected = -0.087, tol = 0.001)
##[PaperB] claim("R3.2", "Results 3.2.1 / tab:s3-community", "community beta0 coverage 0.87",
##[PaperB]       mean(b0$covered), expected = 0.87)
##[PaperB] claim("R3.3", "Results 3.2.1 / tab:s3-community", "community beta1 bias -0.076",
##[PaperB]       mean(b1$est) - b1$truth[1], expected = -0.076, tol = 0.001)
##[PaperB] claim("R3.4", "Results 3.2.1 / tab:s3-community", "community beta1 coverage 0.68",
##[PaperB]       mean(b1$covered), expected = 0.68)

sk <- msom$shrink
##[PaperB] claim("R3.5", "Results 3.2.1", "MSOM closer to truth in 45.3% of pairs",
##[PaperB]       100 * mean(sk$msom_closer), expected = 45.3, tol = 0.1)
##[PaperB] claim("R3.6", "Results 3.2.1", "RMSE 0.247 hierarchical",
##[PaperB]       sqrt(mean((sk$msom - sk$truth)^2)), expected = 0.247, tol = 0.001)
##[PaperB] claim("R3.7", "Results 3.2.1", "RMSE 0.219 independent",
##[PaperB]       sqrt(mean((sk$ssom - sk$truth)^2)), expected = 0.219, tol = 0.001)
##[PaperB] claim("R3.8", "Results 3.2.1",
##[PaperB]       "RMSE higher under pooling for 11 of 12 species",
##[PaperB]       sum(shr$rmse_msom > shr$rmse_ssom), expected = 11, tol = 0.5)

## ---------------------------------------------------------------------
## Homogeneous-slope diagnostic run (v3)
## ---------------------------------------------------------------------
f3 <- "results/msom/msom_replicates_v3.rds"
if (file.exists(f3)) {
  v3 <- readRDS(f3)
  b0_3 <- v3$comm[v3$comm$parameter == "beta0_community", ]
  b1_3 <- v3$comm[v3$comm$parameter == "beta1_community", ]
##[PaperB]   claim("R3.9",  "Results 3.2.1 / Appendix S1",
##[PaperB]         "homogeneous-slope run: beta1 coverage 0.84",
##[PaperB]         mean(b1_3$covered), expected = 0.84)
##[PaperB]   claim("R3.10", "Results 3.2.1 / Appendix S1",
##[PaperB]         "homogeneous-slope run: beta1 bias -0.050",
##[PaperB]         mean(b1_3$est) - b1_3$truth[1], expected = -0.050, tol = 0.001)
##[PaperB]   claim("R3.11", "Results 3.2.1 / Appendix S1",
##[PaperB]         "homogeneous-slope run: beta0 coverage 0.85",
##[PaperB]         mean(b0_3$covered), expected = 0.85)
##[PaperB]   claim("R3.12", "Results 3.2.1 / Appendix S1",
##[PaperB]         "homogeneous-slope run: RMSE 0.241 hierarchical",
##[PaperB]         sqrt(mean((v3$shrink$msom - v3$shrink$truth)^2)),
##[PaperB]         expected = 0.241, tol = 0.001)
##[PaperB]   claim("R3.13", "Results 3.2.1 / Appendix S1",
##[PaperB]         "homogeneous-slope run: RMSE 0.205 independent",
##[PaperB]         sqrt(mean((v3$shrink$ssom - v3$shrink$truth)^2)),
##[PaperB]         expected = 0.205, tol = 0.001)
} else {
  message("msom_replicates_v3.rds not found; v3 claims not checked.")
}

## =====================================================================
## Results, empirical MSOM
## =====================================================================
## An ordering check only: it shows the sum-to-zero constraint has not
## reordered the species, not that the estimates are accurate.
##[PaperB] claim("R4.1", "Results 3.2.2",
##[PaperB]       "constrained estimates preserve the ordering by naive frequency",
##[PaperB]       cor(sp$naive_occ, sp$msom_int), expected = 0.982, tol = 0.003)
##[PaperB] claim("R4.1b", "Results 3.2.2",
##[PaperB]       "independent fits preserve that ordering less closely",
##[PaperB]       cor(sp$naive_occ, sp$ssom_int), expected = 0.889, tol = 0.003)

##[PaperB] claim("R4.2", "Results 3.2.2", "nine of twelve species shrunk toward the mean",
##[PaperB]       sum(sp$shrunk_toward_mean), expected = 9, tol = 0.5)

## `blpw` was defined by the withdrawn inadmissibility block; take the
## value from `sp` directly instead. The stray `mawa` assignment that
## followed belonged to the same block and is dropped.
claim("R4.7", "Results 3.2.2", "BLPW naive occupancy 0.096",
      sp$naive_occ[sp$species == "BLPW"], expected = 0.096)

claim("R4.11", "Results 3.2.2",
      "prevalence ranges from 2.3% to 73.3%",
      100 * max(sp$naive_occ), expected = 73.3, tol = 0.1)

## ---------------------------------------------------------------------
## Convergence under the final mesh
## ---------------------------------------------------------------------
f_conv <- "results/real/ssom_convergence.rds"
if (file.exists(f_conv)) {
  cv <- readRDS(f_conv)
  claim("R2.14", "Results 3.1.2 / Limitations",
        "all twelve species converged under the final mesh",
        sum(cv$converged), expected = 12)
} else {
  message("ssom_convergence.rds not found; convergence claim not checked.")
}

## ---------------------------------------------------------------------
## Mesh resolution and identifiability of the spatial scale
##
## R7.1 is the one that matters. The argument in the Results is not that
## the fine mesh gave a large range but that it gave a range exceeding
## the diameter of the study area, so that correlation never decays
## within the data and the parameter is informed by the prior. If a
## refit ever brought the estimate below that threshold, the argument
## would need rebuilding, and this claim is what would say so.
## ---------------------------------------------------------------------
f_mesh <- "results/mesh/mesh_comparison.csv"
if (file.exists(f_mesh)) {
  mh <- read.csv(f_mesh)
  co <- mh[mh$mesh == "coarse", ]
  fi <- mh[grepl("^fine", mh$mesh), ]

  ## Maximum distance between survey locations, recomputed rather than
  ## taken on trust.
  suppressPackageStartupMessages(library(spOccupancy))
  data(hbefTrends)
  max_d <- max(dist(hbefTrends$coords / 1000))

  claim("R7.1", "Results 3.1.2 / Limitations",
        "fine-mesh range exceeds the maximum inter-site distance",
        as.numeric(all(fi$range > max_d)), expected = 1)
  claim("R7.2", "Results 3.1.2 / Limitations",
        "maximum inter-site distance 7.7 km",
        max_d, expected = 7.74, tol = 0.02)
  claim("R7.3", "Results 3.1.2",
        "coarse-fine range difference dwarfs the fine mesh's own spread",
        abs(co$range - mean(fi$range)) / diff(range(fi$range)),
        expected = 36, tol = 5)
  claim("R7.4", "Results 3.1.2 / tab:mesh",
        "coarse mesh range 3.27 km",
        co$range, expected = 3.27, tol = 0.02)
  claim("R7.5", "Results 3.1.2 / tab:mesh",
        "coarse mesh sigma 3.09",
        co$sigma, expected = 3.09, tol = 0.02)
  claim("R7.6", "Results 3.1.2 / tab:mesh",
        "fine mesh sigma about 14",
        mean(fi$sigma), expected = 14.0, tol = 0.1)
  claim("R7.7", "Results 3.1.2",
        "fine mesh achieved the lower WAIC",
        as.numeric(mean(fi$waic) < co$waic), expected = 1)
  ## The manuscript quotes 42% for the M1 coarse fit and 43% for the
  ## MSOM, whose range is 3.30 rather than 3.27. They are different fits;
  ## check each against its own source.
  claim("R7.8", "Results 3.1.2",
        "M1 coarse range is 42% of the domain diameter",
        100 * co$range / max_d, expected = 42.2, tol = 0.3)
##[PaperB]   claim("R7.9", "Results 3.2.2",
##[PaperB]         "MSOM range is 43% of the domain diameter",
##[PaperB]         100 * 3.30 / max_d, expected = 42.6, tol = 0.3)
} else {
  message("results/mesh/mesh_comparison.csv not found; mesh claims not checked.")
}

## ---------------------------------------------------------------------
## Site-level cross-validation
## ---------------------------------------------------------------------
f_cv <- "results/cv/site_cv_raw.rds"
if (file.exists(f_cv)) {
  cvd <- readRDS(f_cv)
  cvd <- cvd[!is.na(cvd$log_h) & !is.na(cvd$log_s), ]
##[PaperB]   claim("R8.1", "Results 3.2.2 / tab:cv",
##[PaperB]         "independent log score -0.960", mean(cvd$log_s),
##[PaperB]         expected = -0.960, tol = 0.005)
##[PaperB]   claim("R8.2", "Results 3.2.2 / tab:cv",
##[PaperB]         "hierarchical log score -1.053", mean(cvd$log_h),
##[PaperB]         expected = -1.053, tol = 0.005)
##[PaperB]   claim("R8.3", "Results 3.2.2",
##[PaperB]         "paired log-score difference 0.093",
##[PaperB]         mean(cvd$log_s - cvd$log_h), expected = 0.093, tol = 0.002)
##[PaperB]   claim("R8.4", "Results 3.2.2 / tab:cv",
##[PaperB]         "independent Brier 0.110",
##[PaperB]         mean((cvd$prec_s - cvd$observed)^2), expected = 0.110, tol = 0.002)
##[PaperB]   claim("R8.5", "Results 3.2.2 / tab:cv",
##[PaperB]         "hierarchical Brier 0.136",
##[PaperB]         mean((cvd$prec_h - cvd$observed)^2), expected = 0.136, tol = 0.002)
  ## The claim that matters: independent wins for every species, not on
  ## average. If a rerun ever reverses one, the text must change.
  by_sp_cv <- do.call(rbind, lapply(split(cvd, cvd$species), function(d)
    data.frame(better = mean((d$prec_s - d$observed)^2) <
                        mean((d$prec_h - d$observed)^2))))
##[PaperB]   claim("R8.6", "Results 3.2.2 / Abstract / Conclusions",
##[PaperB]         "independent fitting predicts better for all twelve species",
##[PaperB]         sum(by_sp_cv$better), expected = 12)
} else {
  message("results/cv/site_cv_raw.rds not found; CV claims not checked.")
}

## ---------------------------------------------------------------------
## Absent Ovenbird records
## ---------------------------------------------------------------------
ms <- read.csv("results/real/missing_sensitivity.csv")
claim("R6.1", "Results 3.1.2 / tab:missing",
      "biotic coefficient moves by at most 0.02 posterior SD across treatments",
      max(abs(ms$alpha_biotic - ms$alpha_biotic[1]) / ms$sd[1]),
      expected = 0.02, tol = 0.005)
claim("R6.2", "Results 3.1.2",
      "all three treatments use the same 8770 records",
      length(unique(ms$n_used)), expected = 1)
## The text originally said "the same WAIC to the nearest unit". The spread
## is 0.52, which that phrasing obscured while still passing a tolerance of
## 1. The text now states 0.5 and the table carries one decimal place.
claim("R6.3", "Results 3.1.2", "WAIC spread across treatments is 0.5",
      max(ms$waic) - min(ms$waic), expected = 0.5, tol = 0.05)

## ---------------------------------------------------------------------
## Posterior uncertainty paragraph
## ---------------------------------------------------------------------
fig  <- readRDS("results/real/figure_species_data.rds")
ssom <- readRDS("results/real/ssom_summaries.rds")
no   <- setNames(sp$naive_occ, sp$species)

a <- aggregate(cbind(mean, sd) ~ species, fig, mean)
a$naive <- no[as.character(a$species)]

##[PaperB] claim("R5.1", "Results 3.2.2", "probability-scale SD minimum 0.008",
##[PaperB]       min(a$sd), expected = 0.008, tol = 0.0006)
##[PaperB] claim("R5.2", "Results 3.2.2", "probability-scale SD maximum 0.048",
##[PaperB]       max(a$sd), expected = 0.048, tol = 0.0006)
##[PaperB] claim("R5.3", "Results 3.2.2", "SD peaks for BHVI",
##[PaperB]       as.numeric(as.character(a$species[which.max(a$sd)]) == "BHVI"),
##[PaperB]       expected = 1)
##[PaperB] claim("R5.4", "Results 3.2.2", "BHVI mean psi 0.29",
##[PaperB]       a$mean[a$species == "BHVI"], expected = 0.29, tol = 0.006)

## Inverse delta method: SD on the logit scale.
fig$sd_logit <- fig$sd / (fig$mean * (1 - fig$mean))
b <- aggregate(sd_logit ~ species, fig, mean)
b$naive <- no[as.character(b$species)]

##[PaperB] claim("R5.5", "Results 3.2.2", "correlation with prevalence -0.53 over twelve",
##[PaperB]       cor(b$sd_logit, b$naive), expected = -0.53, tol = 0.01)
##[PaperB] claim("R5.5b", "Results 3.2.2", "correlation -0.30 excluding NAWA and AMRE",
##[PaperB]       cor(b$sd_logit[!b$species %in% c("NAWA", "AMRE")],
##[PaperB]           b$naive   [!b$species %in% c("NAWA", "AMRE")]),
##[PaperB]       expected = -0.30, tol = 0.01)
##[PaperB] claim("R5.6", "Results 3.2.2", "NAWA logit-scale SD 0.269",
##[PaperB]       b$sd_logit[b$species == "NAWA"], expected = 0.269, tol = 0.001)
##[PaperB] claim("R5.6b", "Results 3.2.2", "AMRE logit-scale SD 0.262",
##[PaperB]       b$sd_logit[b$species == "AMRE"], expected = 0.262, tol = 0.001)
##[PaperB] claim("R5.7", "Results 3.2.2", "BTNW logit-scale SD 0.243",
##[PaperB]       b$sd_logit[b$species == "BTNW"], expected = 0.243, tol = 0.001)
##[PaperB] claim("R5.7b", "Results 3.2.2", "REVI logit-scale SD 0.239",
##[PaperB]       b$sd_logit[b$species == "REVI"], expected = 0.239, tol = 0.001)
## The minimum is BHVI, not the most common species. The manuscript says so
## because an earlier draft implied the low end of the range was BTNW.
##[PaperB] claim("R5.7c", "Results 3.2.2", "minimum logit-scale SD belongs to BHVI",
##[PaperB]       as.numeric(as.character(b$species[which.min(b$sd_logit)]) == "BHVI"),
##[PaperB]       expected = 1)
##[PaperB] claim("R5.8", "Results 3.2.2", "hierarchical logit-scale SD range 1.13-fold",
##[PaperB]       max(b$sd_logit) / min(b$sd_logit), expected = 1.13, tol = 0.02)
##[PaperB] claim("R5.8b", "Results 3.2.2", "ten-species band runs 0.237 to 0.247",
##[PaperB]       max(b$sd_logit[!b$species %in% c("NAWA", "AMRE")]),
##[PaperB]       expected = 0.247, tol = 0.001)

## Independent-model intercept standard deviations.
s2 <- vapply(ssom, function(f) f$fixed["Int_occ", "sd"], numeric(1))
##[PaperB] claim("R5.9",  "Results 3.2.2", "independent SD minimum 0.209",
##[PaperB]       min(s2), expected = 0.209, tol = 0.001)
##[PaperB] claim("R5.10", "Results 3.2.2", "independent SD maximum 0.560",
##[PaperB]       max(s2), expected = 0.560, tol = 0.001)
##[PaperB] claim("R5.11", "Results 3.2.2", "independent SD range 2.7-fold",
##[PaperB]       max(s2) / min(s2), expected = 2.68, tol = 0.02)
##[PaperB] claim("R5.12", "Results 3.2.2", "smallest independent SD belongs to BTNW",
##[PaperB]       as.numeric(names(s2)[which.min(s2)] == "BTNW"), expected = 1)
##[PaperB] claim("R5.13", "Results 3.2.2", "second smallest belongs to BLBW",
##[PaperB]       as.numeric(names(sort(s2))[2] == "BLBW"), expected = 1)
## The text first said seven; the count is nine of the eleven species
## recorded more often than NAWA. The corrected figure strengthens the
## point rather than weakening it.
##[PaperB] claim("R5.14", "Results 3.2.2",
##[PaperB]       "NAWA independent SD smaller than nine of eleven more frequent species",
##[PaperB]       sum(s2 > s2["NAWA"] & no[names(s2)] > no["NAWA"]), expected = 9, tol = 0.5)

## =====================================================================
## Report
## =====================================================================
res <- do.call(rbind, lapply(CLAIMS, function(c) data.frame(
  id = c$id, where = c$where, claim = c$text,
  value = c$value, expected = c$expected,
  pass = isTRUE(abs(c$value - c$expected) <= c$tol)
)))

cat("\n", strrep("=", 78), "\n", sep = "")
print(res[, c("id", "value", "expected", "pass", "claim")],
      row.names = FALSE, digits = 5, right = FALSE)
cat(strrep("=", 78), "\n", sep = "")
cat(sprintf("%d of %d claims verified.\n", sum(res$pass), nrow(res)))

if (any(!res$pass)) {
  cat("\nFAILED -- the manuscript states these numbers but the results do not:\n")
  print(res[!res$pass, c("id", "where", "claim", "value", "expected")],
        row.names = FALSE, digits = 5, right = FALSE)
  stop("Manuscript claims do not match results. Fix the text, not this script.")
}

write.csv(res, "results/claim_verification.csv", row.names = FALSE)
