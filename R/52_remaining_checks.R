## =====================================================================
## 52_remaining_checks.R -- the cheap outstanding checks for PaperB
## =====================================================================
## Five checks, each reading stored output only (no model is refitted):
##
##   1. Mesh and field constraint of the stored application fits
##      (PaperB Sec. 5.1 TBD): node count of the mesh rebuilt from the
##      script specification, node count implied by the fitted field, and
##      whether the field's posterior mean integrates to zero per year
##      (the signature of constr = TRUE in inla.spde2.pcmatern).
##   2. Pooled detection coefficients of the hierarchical model
##      (Table S3.4 TBD).
##   3. Posterior predictive check of detections per recorded site-year,
##      by species and model class, from posterior draws of the detection
##      coefficients (Sec. 5.3 TBD).
##   4. Cluster (site-level) bootstrap standard errors for the paired
##      cross-validation differences, from site_cv_raw.rds (Sec. 3.5 and
##      5.7 TBDs).
##   5. The v4 test replicate: -5.65 or -5.55? (Sec. 4.4 TBD). Looks for
##      stored output; reports clearly if there is none.
##
## Requires: results/real/fits.rds, site_cv_raw.rds (searched for under
## results/, it lives in results/cv/) and,
## for check 5, anything under results/msom/ from v4.
## Each check is wrapped so that a failure in one does not stop the rest.
## NOT RUN BY THE AUTHOR OF THIS DRAFT.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUT_REAL <- "results/real"
OUT_MSOM <- "results/msom"
N_DRAWS  <- 1000
B_BOOT   <- 2000
SEED     <- 20260201L

fits   <- readRDS(file.path(OUT_REAL, "fits.rds"))
m_msom <- fits$msom
ssom   <- fits$ssom

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K        <- dim(hbefTrends$y)[4]
N_ROW    <- N_SITE * N_YEAR

## helpers identical to 63_avg_occupancy.R / 65_site_cv.R
stack_species <- function(sp) {
  y <- hbefTrends$y[sp, , , ]
  do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
}
stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  M <- (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
  M[is.na(M)] <- 0
  M
}
coords <- hbefTrends$coords / 1000
Xday   <- stack_det_cov("day")
Xtod   <- stack_det_cov("tod")
Y_all  <- lapply(SP_NAMES, stack_species); names(Y_all) <- SP_NAMES

section <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n",
                             strrep("=", 78), "\n", sep = "")
safely <- function(expr, label) {
  tryCatch(expr, error = function(e) {
    cat("!! ", label, " failed: ", conditionMessage(e), "\n", sep = "")
    NULL
  })
}
results <- list()

## ---------------------------------------------------------------------
## 1. Mesh and field constraint
## ---------------------------------------------------------------------
section("1. Mesh and field constraint")
results$mesh <- safely({
  mesh <- fmesher::fm_mesh_2d(
    boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
    max.edge = c(0.35, 1.2), min.angle = 21,
    offset = c(0.05, 1.5), cutoff = 0.35)
  n_fit_h <- nrow(m_msom$summary.random$spatialfield) / N_YEAR
  n_fit_s <- nrow(ssom[[1]]$summary.random$spatialfield) / N_YEAR
  cat(sprintf("mesh rebuilt from script spec: %d nodes\n", mesh$n))
  cat(sprintf("nodes implied by fitted field: hierarchical %g, independent (%s) %g\n",
              n_fit_h, SP_NAMES[1], n_fit_s))

  ## integral of the posterior-mean field per year, relative to its scale
  c0 <- Matrix::diag(fmesher::fm_fem(mesh)$c0)
  fm <- m_msom$summary.random$spatialfield$mean
  stopifnot(length(fm) == mesh$n * N_YEAR)
  per_year <- t(vapply(seq_len(N_YEAR), function(tt) {
    v <- fm[((tt - 1) * mesh$n + 1):(tt * mesh$n)]
    c(weighted_mean = sum(c0 * v) / sum(c0),
      unweighted_mean = mean(v),
      sd = sd(v))
  }, numeric(3)))
  per_year <- data.frame(year = seq_len(N_YEAR), per_year)
  per_year$rel_weighted <- abs(per_year$weighted_mean) / per_year$sd
  print(per_year, digits = 3, row.names = FALSE)
  all_w <- sum(rep(c0, N_YEAR) * fm) / (N_YEAR * sum(c0))
  cat(sprintf("area-weighted mean over all years: %.2e (field SD ~ %.3f)\n",
              all_w, sd(fm)))
  cat("Reading: if the weighted mean is ~0 relative to the SD (say < 1e-3)\n",
      "for every year, or over all years jointly, the field was fitted with\n",
      "an integrate-to-zero constraint. Values of the order of the SD mean\n",
      "it was not.\n", sep = "")
  list(mesh_n = mesh$n, n_fit_h = n_fit_h, n_fit_s = n_fit_s,
       per_year = per_year, all_years_weighted_mean = all_w)
}, "check 1")

## ---------------------------------------------------------------------
## 2. Pooled detection coefficients of the hierarchical model
## ---------------------------------------------------------------------
section("2. Detection coefficients, hierarchical model (pooled)")
results$det_coef <- safely({
  hp <- m_msom$summary.hyperpar
  cat("summary.hyperpar rows matching 'beta':\n")
  print(hp[grep("beta", rownames(hp), ignore.case = TRUE),
           c("mean", "sd", "0.025quant", "0.975quant")], digits = 3)
  set.seed(SEED)
  dr <- det_draws(m_msom, n = N_DRAWS)
  stopifnot(ncol(dr) >= 3)
  tab <- data.frame(
    term = c("intercept", "survey day", "time of day"),
    mean = colMeans(dr[, 1:3]),
    q025 = apply(dr[, 1:3], 2, quantile, 0.025),
    q975 = apply(dr[, 1:3], 2, quantile, 0.975))
  cat("\nfrom det_draws() (as used for all derived quantities):\n")
  print(tab, digits = 3, row.names = FALSE)
  cat(sprintf("plogis(intercept) = %.3f\n", plogis(tab$mean[1])))
  tab
}, "check 2")

## ---------------------------------------------------------------------
## 3. Posterior predictive check: detections per recorded site-year
## ---------------------------------------------------------------------
section("3. Detections per recorded site-year: observed vs posterior predictive")
results$ppc <- safely({
  ## expected detections given at least one, for one draw and one set of rows
  cond_mean <- function(P, made, rows) {
    Pm <- P[rows, , drop = FALSE]; Mm <- made[rows, , drop = FALSE]
    Pm[!Mm] <- 0
    mean(rowSums(Pm) / (1 - exp(rowSums(log1p(-Pm)))))   # vectorised prod(1 - p)
  }
  p_mat <- function(b) plogis(b[1] + b[2] * Xday + b[3] * Xtod)

  set.seed(SEED)
  dr_h <- det_draws(m_msom, n = N_DRAWS)
  made_l <- lapply(Y_all, \(Y) !is.na(Y))
  rec_l  <- lapply(Y_all, \(Y) which(rowSums(Y == 1, na.rm = TRUE) >= 1))
  obs    <- vapply(SP_NAMES, \(sp) mean(rowSums(Y_all[[sp]][rec_l[[sp]], , drop = FALSE],
                                                  na.rm = TRUE)), numeric(1))

  ## hierarchical: one detection matrix per draw, shared by all species
  pred_h <- matrix(NA_real_, N_DRAWS, length(SP_NAMES), dimnames = list(NULL, SP_NAMES))
  for (d in seq_len(N_DRAWS)) {
    P <- p_mat(dr_h[d, 1:3])
    for (sp in SP_NAMES) pred_h[d, sp] <- cond_mean(P, made_l[[sp]], rec_l[[sp]])
  }
  ## independent: species-specific draws
  pred_s <- matrix(NA_real_, N_DRAWS, length(SP_NAMES), dimnames = list(NULL, SP_NAMES))
  for (sp in SP_NAMES) {
    set.seed(SEED)
    dr_s <- det_draws(ssom[[sp]], n = N_DRAWS)
    for (d in seq_len(N_DRAWS))
      pred_s[d, sp] <- cond_mean(p_mat(dr_s[d, 1:3]), made_l[[sp]], rec_l[[sp]])
  }
  q <- function(M, p) apply(M, 2, quantile, p)
  tab <- data.frame(
    species  = SP_NAMES,
    n_rec    = vapply(rec_l, length, integer(1)),
    observed = obs,
    hier_mean = colMeans(pred_h), hier_lo = q(pred_h, .025), hier_hi = q(pred_h, .975),
    ind_mean  = colMeans(pred_s), ind_lo  = q(pred_s, .025), ind_hi  = q(pred_s, .975))
  tab$obs_in_hier <- tab$observed >= tab$hier_lo & tab$observed <= tab$hier_hi
  tab$obs_in_ind  <- tab$observed >= tab$ind_lo  & tab$observed <= tab$ind_hi
  naive <- vapply(SP_NAMES, \(sp) { Y <- Y_all[[sp]]; ok <- rowSums(!is.na(Y)) > 0
                   mean(rowSums(Y[ok, , drop = FALSE] == 1, na.rm = TRUE) > 0) }, numeric(1))
  tab <- tab[order(naive[tab$species]), ]
  print(tab, digits = 3, row.names = FALSE)
  cat("Intervals are for the posterior expectation of the mean, not for a\n",
      "replicated data set, so they are narrower than a full predictive\n",
      "interval; with n_rec in the hundreds the difference is small.\n", sep = "")
  tab
}, "check 3")

## ---------------------------------------------------------------------
## 4. Cluster bootstrap SE for the paired CV differences
## ---------------------------------------------------------------------
section("4. Site-level cross-validation: cluster standard errors")
results$cv <- safely({
  f_cv <- c(file.path(OUT_REAL, "site_cv_raw.rds"),
            file.path("results", "cv", "site_cv_raw.rds"),
            list.files("results", pattern = "^site_cv_raw\\.rds$",
                       recursive = TRUE, full.names = TRUE))
  f_cv <- f_cv[file.exists(f_cv)]
  if (!length(f_cv)) stop("site_cv_raw.rds not found under results/")
  cat("reading ", f_cv[1], "\n", sep = "")
  cv <- readRDS(f_cv[1])
  cv <- cv[!is.na(cv$log_h) & !is.na(cv$log_s), ]
  cv$site    <- ((cv$row - 1) %% N_SITE) + 1
  cv$d_log   <- cv$log_h - cv$log_s
  cv$d_brier <- (cv$observed - cv$prec_h)^2 - (cv$observed - cv$prec_s)^2
  cv$one     <- 1

  cat(sprintf("records %d | sites %d | species %d\n",
              nrow(cv), length(unique(cv$site)), length(unique(cv$species))))
  cat(sprintf("means: log score hier %.3f indep %.3f | Brier hier %.3f indep %.3f\n",
              mean(cv$log_h), mean(cv$log_s),
              mean((cv$observed - cv$prec_h)^2), mean((cv$observed - cv$prec_s)^2)))

  site_sum <- aggregate(cbind(d_log, d_brier, one) ~ site, cv, sum)
  boot_stat <- function(idx) {
    s <- site_sum[idx, ]
    c(d_log = sum(s$d_log) / sum(s$one), d_brier = sum(s$d_brier) / sum(s$one))
  }
  set.seed(SEED)
  bt <- t(replicate(B_BOOT, boot_stat(sample.int(nrow(site_sum), replace = TRUE))))
  est <- c(d_log = mean(cv$d_log), d_brier = mean(cv$d_brier))
  naive_se <- c(sd(cv$d_log), sd(cv$d_brier)) / sqrt(nrow(cv))
  out <- data.frame(
    difference = c("log score (hier - indep)", "Brier (hier - indep)"),
    estimate   = est,
    se_naive   = naive_se,
    se_cluster = apply(bt, 2, sd),
    ci_lo      = apply(bt, 2, quantile, 0.025),
    ci_hi      = apply(bt, 2, quantile, 0.975))
  out$inflation <- out$se_cluster / out$se_naive
  print(out, digits = 3, row.names = FALSE)

  ## by species, cluster SE of the log-score difference
  by_sp <- do.call(rbind, lapply(split(cv, cv$species), function(d) {
    ss <- aggregate(cbind(d_log, one) ~ site, d, sum)
    b  <- replicate(B_BOOT, { i <- sample.int(nrow(ss), replace = TRUE)
                              sum(ss$d_log[i]) / sum(ss$one[i]) })
    data.frame(species = d$species[1], d_log = mean(d$d_log),
               se_cluster = sd(b), ci_lo = quantile(b, .025), ci_hi = quantile(b, .975))
  }))
  cat("\nby species (log-score difference, cluster bootstrap):\n")
  print(by_sp, digits = 3, row.names = FALSE)
  cat("\nNote: 65_site_cv.R computes the log score at the posterior mean of\n",
      "psi and of the detection coefficients (plug-in), not as a posterior\n",
      "predictive density; that is how it should be described.\n", sep = "")
  list(overall = out, by_species = by_sp)
}, "check 4")

## ---------------------------------------------------------------------
## 5. The v4 test replicate
## ---------------------------------------------------------------------
section("5. v4 test replicate: -5.65 or -5.55?")
results$v4 <- safely({
  cand <- c(file.path(OUT_MSOM, "msom_replicates_v4.rds"),
            list.files(file.path(OUT_MSOM, "checkpoints_v4"),
                       pattern = "\\.rds$", full.names = TRUE),
            list.files(OUT_MSOM, pattern = "v4.*\\.rds$|test.*\\.rds$",
                       full.names = TRUE, recursive = TRUE))
  cand <- unique(cand[file.exists(cand)])
  if (!length(cand)) {
    cat("No stored v4 output found under ", OUT_MSOM, ".\n", sep = "")
    cat("The test replicate was evidently run interactively and not saved.\n",
        "Options: (i) describe the result qualitatively in the paper, since\n",
        "run (d) is being redesigned; or (ii) rerun it (about 35 min):\n",
        "  RUN_FULL <- FALSE; source('R/50_run_msom_sim_v4.R')\n",
        "  RNGkind(\"L'Ecuyer-CMRG\"); r1 <- run_msom_rep(1)\n",
        "  r1$shrink[which.min(abs(r1$shrink$truth + 3.67)), ]\n",
        "and note that the original test may have used the default generator,\n",
        "so the rerun need not reproduce either number.\n", sep = "")
    NULL
  } else {
    cat("found:\n"); cat(paste0("  ", cand), sep = "\n")
    hits <- do.call(rbind, lapply(cand, function(f) {
      x  <- readRDS(f)
      sh <- if (is.list(x) && !is.null(x$shrink)) x$shrink else
            if (is.list(x)) do.call(rbind, lapply(x, `[[`, "shrink")) else NULL
      if (is.null(sh) || !"truth" %in% names(sh)) return(NULL)
      h <- sh[abs(sh$truth - (-3.67)) < 0.005, ]
      if (nrow(h)) cbind(file = basename(f), h) else NULL
    }))
    if (is.null(hits) || !nrow(hits)) {
      cat("No species with generating intercept -3.67 in the stored files.\n")
    } else {
      print(hits, digits = 4, row.names = FALSE)
    }
    hits
  }
}, "check 5")

saveRDS(results, file.path(OUT_REAL, "remaining_checks.rds"))
cat("\nAll results saved to ", file.path(OUT_REAL, "remaining_checks.rds"), "\n", sep = "")
