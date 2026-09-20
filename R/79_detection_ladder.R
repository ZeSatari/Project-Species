## =====================================================================
## 79_detection_ladder.R -- does relaxing the pooling change the picture?
## =====================================================================
## PaperB, Section 5.8. The hierarchical and independent classes differ in
## four respects (Sec. 5.1): (1) partial pooling of the occupancy
## intercepts, (2) pooled vs species-specific detection, (3) a common vs a
## species-specific elevation response, (4) one shared field vs one field
## per species. This script removes (3) and then (2), as far as the
## ten-hyperparameter ceiling allows, leaving (1) and (4).
##
##   M1  current model: pooled detection, common elevation response
##   M2  M1 + species-specific elevation slopes (iid, sum-to-zero)
##   M3  M2 + grouped detection intercepts (G = 3): 5 detection
##       coefficients, within the ceiling of 10
##
## Detection groups come from the conditional estimates p_hat_i of script
## 68, which are free of the occupancy model. They use the same data as
## the fits, so two alternative groupings are refitted on the full data as
## a sensitivity check (Sec. 5.8): G = 2, and G = 3 with the two borderline
## species moved up.
##
## For each model: average occupancy per species, the predictive check on
## detections per recorded site-year, and site-level cross-validation on
## the SAME folds as 65_site_cv.R. The independent-fit scores are read
## from site_cv_raw.rds rather than recomputed, so only the hierarchical
## models are fitted: 3 full fits + 3 x 5 fold fits (+2 sensitivity fits).
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT. Check the first fold's messages
## before leaving it running.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR      <- "results/ladder"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "ladder")

K_FOLDS     <- 5
SEED_FOLDS  <- 20260401L      # must match 65_site_cv.R
SEED        <- 20260201L
N_DRAWS     <- 500
## Stages. Set any of these to FALSE before sourcing to skip that stage,
## e.g. a first pass with only the full-data fits:
##     DO_CV <- FALSE; DO_SENS <- FALSE; source("R/79_detection_ladder.R")
## Then, once the fits look right, the whole thing:
##     rm(DO_CV, DO_SENS); source("R/79_detection_ladder.R")
if (!exists("DO_FULL")) DO_FULL <- TRUE    # 3 fits
if (!exists("DO_SENS")) DO_SENS <- TRUE    # 2 more fits
if (!exists("DO_CV"))   DO_CV   <- TRUE    # 15 more fits, the slow stage
F_CV_INDEP  <- "results/cv/site_cv_raw.rds"
F_PHAT      <- "results/real/detection_heterogeneity.csv"

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP     <- length(SP_NAMES)
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K_VIS    <- dim(hbefTrends$y)[4]
N_ROW    <- N_SITE * N_YEAR

## ---------------------------------------------------------------------
## Data, exactly as in 65_site_cv.R
## ---------------------------------------------------------------------
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
coords  <- hbefTrends$coords / 1000
elev_s  <- as.numeric(scale(hbefTrends$occ.covs$elev))
site_id <- rep(seq_len(N_SITE), times = N_YEAR)
time_id <- rep(seq_len(N_YEAR), each = N_SITE)
Xday    <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")
Y_all   <- lapply(SP_NAMES, stack_species); names(Y_all) <- SP_NAMES

mesh <- fmesher::fm_mesh_2d(
  boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
  max.edge = c(0.35, 1.2), min.angle = 21,
  offset = c(0.05, 1.5), cutoff = 0.35)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))
message("mesh: ", mesh$n, " nodes")

## Folds: the same partition as 65_site_cv.R (same seed, same generator).
set.seed(SEED_FOLDS)
fold_of_site <- sample(rep_len(seq_len(K_FOLDS), N_SITE))

## ---------------------------------------------------------------------
## Detection groups from the conditional estimates
## ---------------------------------------------------------------------
ph <- read.csv(F_PHAT)
ph <- ph[match(SP_NAMES, ph$species), ]
stopifnot(!any(is.na(ph$p_hat)))

grp3     <- cut(ph$p_hat, c(-Inf, 0.30, 0.50, Inf), labels = FALSE)
grp2     <- cut(ph$p_hat, c(-Inf, 0.35, Inf), labels = FALSE)
grp3_alt <- grp3
grp3_alt[SP_NAMES %in% c("BLPW", "MAWA")] <- 3L   # borderline species moved up

cat("\ndetection groups (G = 3):\n")
print(data.frame(species = SP_NAMES, p_hat = round(ph$p_hat, 3),
                 g3 = grp3, g2 = grp2, g3_alt = grp3_alt), row.names = FALSE)

## group indicator columns for the stacked (species-major) design:
## one column per non-reference group, constant across visits
det_group_cols <- function(g, n_row = N_ROW) {
  levs <- sort(unique(g))[-1]
  lapply(levs, function(l) {
    ind <- rep(as.numeric(g == l), each = n_row)   # species-major
    matrix(rep(ind, K_VIS), ncol = K_VIS)
  })
}

## detection ceiling, from 00_setup.R if it exports one
if (!exists("DET_CEILING")) DET_CEILING <- 10L

## ---------------------------------------------------------------------
## Model specification. `slopes`: species-specific elevation terms.
## `groups`: NULL for pooled detection, else a group vector of length N_SP.
## ---------------------------------------------------------------------
make_fit <- function(Y_list, slopes, groups, held_rows = integer(0)) {
  Y_h <- do.call(rbind, lapply(SP_NAMES, function(sp) {
    y <- Y_list[[sp]]; if (length(held_rows)) y[held_rows, ] <- NA; y
  }))
  det_list <- list(do.call(rbind, rep(list(Xday), N_SP)),
                   do.call(rbind, rep(list(Xtod), N_SP)))
  if (!is.null(groups)) det_list <- c(det_list, det_group_cols(groups))
  Xd <- build_det_X(det_list, K = K_VIS)
  ## build_det_X() lays out K columns per coefficient, so the design has
  ## K_VIS * n_coef columns. The ceiling applies to the COEFFICIENTS.
  n_coef <- ncol(Xd) / K_VIS
  stopifnot(n_coef == round(n_coef))
  if (n_coef > DET_CEILING)
    stop("detection predictor needs ", n_coef, " coefficients; ceiling is ",
         DET_CEILING)

  spde <- INLA::inla.spde2.pcmatern(mesh,
            prior.range = PRIORS$spde_real$range,
            prior.sigma = PRIORS$spde_real$sigma, constr = TRUE)
  Xocc <- data.frame(
    Int_occ = 1,
    elev_s  = rep(elev_s[site_id], times = N_SP),
    elev_s2 = rep(elev_s[site_id]^2, times = N_SP),
    species_id = rep(seq_len(N_SP), each = N_ROW),
    x = rep(coords[site_id, 1], times = N_SP),
    y = rep(coords[site_id, 2], times = N_SP),
    time = rep(time_id, times = N_SP))
  Xocc$sp_e1 <- Xocc$species_id
  Xocc$sp_e2 <- Xocc$species_id
  stf <- make_st_field(mesh, spde, Xocc[, c("x", "y")],
                       time = Xocc$time, n_time = N_YEAR)

  f <- if (slopes)
    inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
      f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(sp_e1, elev_s,  model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(sp_e2, elev_s2, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(spatialfield, model = spde, A.local = stf$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))
  else
    inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
      f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(spatialfield, model = spde, A.local = stf$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))

  fit <- fit_occupancy(f, c(as.list(Xocc), stf$idx, list(Y = Y_h, X = Xd)),
                       Xd, compute_cv = FALSE)
  list(fit = fit, n_det_coef = n_coef, n_det_col = ncol(Xd))
}

MODELS <- list(
  M1 = list(slopes = FALSE, groups = NULL,
            label = "pooled detection, common elevation"),
  M2 = list(slopes = TRUE,  groups = NULL,
            label = "+ species-specific elevation"),
  M3 = list(slopes = TRUE,  groups = grp3,
            label = "+ grouped detection (G = 3)"))

## ---------------------------------------------------------------------
## Derived quantities from a fitted model
## ---------------------------------------------------------------------
## Detection probability per row and visit for species i.
## `b` is one row of detection coefficients (a draw, or their mean, which
## is what 65_site_cv.R uses for scoring).
det_p <- function(b, groups, i) {
  eta <- b[1] + b[2] * Xday + b[3] * Xtod
  if (!is.null(groups)) {
    levs <- sort(unique(groups))[-1]
    j <- match(groups[i], levs)
    if (!is.na(j)) eta <- eta + b[3 + j]
  }
  plogis(eta)
}

## Posterior mean occupancy for every species, N_ROW x N_SP.
## The field is built once per draw and reused across species, which is
## what makes this affordable: 500 draws, not 500 x 12.
psi_mean_all <- function(s, slopes) {
  nm  <- rownames(s[[1]]$latent)
  i1  <- function(pat) { j <- grep(pat, nm); stopifnot(length(j) == 1); j }
  i_f <- grep("^spatialfield:", nm); n_spde <- length(i_f) / N_YEAR
  idx <- list(int = i1("^Int_occ:"), e1 = i1("^elev_s:"), e2 = i1("^elev_s2:"),
              sp = vapply(seq_len(N_SP), \(i) i1(sprintf("^species_id:%d$", i)), 1L))
  if (slopes) {
    idx$s1 <- vapply(seq_len(N_SP), \(i) i1(sprintf("^sp_e1:%d$", i)), 1L)
    idx$s2 <- vapply(seq_len(N_SP), \(i) i1(sprintf("^sp_e2:%d$", i)), 1L)
  }
  ev <- elev_s[site_id]; ev2 <- ev^2
  acc <- matrix(0, N_ROW, N_SP)
  for (d in seq_along(s)) {
    v <- s[[d]]$latent
    fld <- numeric(N_ROW)
    for (tt in seq_len(N_YEAR)) {
      sl <- i_f[((tt - 1) * n_spde + 1):(tt * n_spde)]
      fld[time_id == tt] <- as.numeric(A_site %*% v[sl])
    }
    base <- v[idx$int] + fld
    for (i in seq_len(N_SP)) {
      b1 <- v[idx$e1]; b2 <- v[idx$e2]
      if (slopes) { b1 <- b1 + v[idx$s1[i]]; b2 <- b2 + v[idx$s2[i]] }
      acc[, i] <- acc[, i] + plogis(base + v[idx$sp[i]] + b1 * ev + b2 * ev2)
    }
  }
  acc / length(s)
}

summarise_full <- function(m, spec, tag) {
  s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = m$fit,
                                    seed = SEED, num.threads = "1:1")
  psi <- psi_mean_all(s, spec$slopes)
  dr  <- det_draws(m$fit, n = N_DRAWS)

  do.call(rbind, lapply(seq_len(N_SP), function(i) {
    Y    <- Y_all[[SP_NAMES[i]]]
    made <- !is.na(Y)
    surv <- rowSums(made) > 0
    rec  <- which(rowSums(Y == 1, na.rm = TRUE) >= 1)
    naive <- mean(rowSums(Y[surv, , drop = FALSE] == 1, na.rm = TRUE) > 0)
    ## predictive check, over draws of the detection coefficients
    pred <- vapply(seq_len(N_DRAWS), function(d) {
      P <- det_p(dr[d, ], spec$groups, i)
      P[!made] <- 0
      Pr <- P[rec, , drop = FALSE]
      mean(rowSums(Pr) / (1 - exp(rowSums(log1p(-Pr)))))
    }, numeric(1))
    data.frame(model = tag, species = SP_NAMES[i], naive = naive,
               psi_bar = mean(psi[surv, i]),
               ratio   = mean(psi[surv, i]) / naive,
               obs_det_per_rec = mean(rowSums(Y[rec, , drop = FALSE], na.rm = TRUE)),
               pred_det_per_rec = mean(pred),
               pred_lo = quantile(pred, .025), pred_hi = quantile(pred, .975))
  }))
}

## ---------------------------------------------------------------------
## Full-data fits
## ---------------------------------------------------------------------
if (!DO_FULL && !file.exists(file.path(OUTDIR, "ladder_fits.rds")))
  stop("DO_FULL = FALSE but no stored fits at ", OUTDIR,
       "/ladder_fits.rds; run the full-data stage first.")

full <- list(); summ <- list()
if (DO_FULL) for (tag in names(MODELS)) {
  message("\n=== full fit: ", tag, " (", MODELS[[tag]]$label, ") ===")
  t0 <- Sys.time()
  full[[tag]] <- make_fit(Y_all, MODELS[[tag]]$slopes, MODELS[[tag]]$groups)
  message(sprintf("  %d detection coefficients (%d design columns) | %.1f min",
                  full[[tag]]$n_det_coef, full[[tag]]$n_det_col,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  summ[[tag]] <- summarise_full(full[[tag]], MODELS[[tag]], tag)
  print(summ[[tag]][, c("species", "naive", "psi_bar", "ratio",
                        "obs_det_per_rec", "pred_det_per_rec")],
        row.names = FALSE, digits = 3)
}
if (DO_FULL) {
  saveRDS(full, file.path(OUTDIR, "ladder_fits.rds"))
  summ_df <- do.call(rbind, summ)
  write.csv(summ_df, file.path(OUTDIR, "ladder_full_summary.csv"), row.names = FALSE)

  ## hyperparameters of each model, side by side
  cat("\n=== hyperparameters ===\n")
  for (tag in names(MODELS)) {
    cat("\n", tag, ": ", MODELS[[tag]]$label, "\n", sep = "")
    print(full[[tag]]$fit$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")],
          digits = 3)
  }
} else {
  full    <- readRDS(file.path(OUTDIR, "ladder_fits.rds"))
  summ_df <- read.csv(file.path(OUTDIR, "ladder_full_summary.csv"))
  message("full-data stage skipped; stored fits loaded.")
}

## ---------------------------------------------------------------------
## Sensitivity: alternative groupings, full data only
## ---------------------------------------------------------------------
if (DO_SENS) {
  for (tag in c("M3_G2", "M3_alt")) {
    g <- if (tag == "M3_G2") grp2 else grp3_alt
    message("\n=== sensitivity fit: ", tag, " ===")
    fm <- make_fit(Y_all, TRUE, g)
    sm <- summarise_full(fm, list(slopes = TRUE, groups = g), tag)
    summ_df <- rbind(summ_df, sm)
    print(sm[, c("species", "ratio", "obs_det_per_rec", "pred_det_per_rec")],
          row.names = FALSE, digits = 3)
  }
  write.csv(summ_df, file.path(OUTDIR, "ladder_full_summary.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------------
## Cross-validation on the folds of 65_site_cv.R
## ---------------------------------------------------------------------
log_score <- function(psi, p, yrow) {
  obs <- !is.na(yrow)
  if (!any(obs)) return(NA_real_)
  lp <- sum(yrow[obs] * log(p[obs]) + (1 - yrow[obs]) * log(1 - p[obs]))
  if (any(yrow[obs] == 1)) log(psi) + lp else log(psi * exp(lp) + (1 - psi))
}
p_recorded <- function(psi, p, yrow) {
  obs <- !is.na(yrow); psi * (1 - prod(1 - p[obs]))
}

if (!DO_CV) {
  cat("\nCV stage skipped (DO_CV = FALSE). Full-data results are in ",
      OUTDIR, ".\n", sep = "")
} else {

rows <- list()
for (k in seq_len(K_FOLDS)) {
  held_rows <- which(site_id %in% which(fold_of_site == k))
  message(sprintf("\nfold %d: %d site-years withheld", k, length(held_rows)))
  for (tag in names(MODELS)) {
    spec <- MODELS[[tag]]
    t0 <- Sys.time()
    m  <- make_fit(Y_all, spec$slopes, spec$groups, held_rows = held_rows)
    s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = m$fit,
                                      seed = SEED, num.threads = "1:1")
    psi <- psi_mean_all(s, spec$slopes)[held_rows, , drop = FALSE]
    b   <- colMeans(det_draws(m$fit, n = N_DRAWS))   # plug-in, as in 65
    message(sprintf("  %s fitted and scored in %.1f min", tag,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))

    for (i in seq_len(N_SP)) {
      P  <- det_p(b, spec$groups, i)
      Yh <- Y_all[[SP_NAMES[i]]][held_rows, , drop = FALSE]
      for (r in seq_along(held_rows)) {
        yr <- Yh[r, ]
        if (all(is.na(yr))) next
        p <- P[held_rows[r], ]
        rows[[length(rows) + 1]] <- data.frame(
          model = tag, fold = k, species = SP_NAMES[i], row = held_rows[r],
          observed = as.integer(any(yr == 1, na.rm = TRUE)),
          log_m  = log_score(psi[r, i], p, yr),
          prec_m = p_recorded(psi[r, i], p, yr))
      }
    }
  }
}
cv_lad <- do.call(rbind, rows)
saveRDS(cv_lad, file.path(OUTDIR, "ladder_cv_raw.rds"))

## ---------------------------------------------------------------------
## Compare with the independent fits of 65_site_cv.R
## ---------------------------------------------------------------------
cv_ind <- readRDS(F_CV_INDEP)
cmp <- merge(cv_lad, cv_ind[, c("fold", "species", "row", "observed",
                                "log_s", "prec_s", "log_h", "prec_h")],
             by = c("fold", "species", "row"))
if (nrow(cmp) != nrow(cv_lad))
  warning("merge lost rows: ", nrow(cv_lad) - nrow(cmp),
          ". Are the folds identical to 65_site_cv.R?")
stopifnot(all(cmp$observed.x == cmp$observed.y))

cmp$site    <- ((cmp$row - 1) %% N_SITE) + 1
cmp$d_log   <- cmp$log_m - cmp$log_s
cmp$d_brier <- (cmp$observed.x - cmp$prec_m)^2 - (cmp$observed.x - cmp$prec_s)^2
cmp$one     <- 1

boot_se <- function(d, col) {
  ss <- aggregate(as.formula(sprintf("cbind(%s, one) ~ site", col)), d, sum)
  b  <- replicate(2000, { i <- sample.int(nrow(ss), replace = TRUE)
                          sum(ss[[col]][i]) / sum(ss$one[i]) })
  c(est = mean(d[[col]]), se = sd(b), lo = quantile(b, .025), hi = quantile(b, .975))
}
set.seed(SEED)
cv_tab <- do.call(rbind, lapply(split(cmp, cmp$model), function(d) {
  data.frame(model = d$model[1],
             log_hier = mean(d$log_m), log_indep = mean(d$log_s),
             t(boot_se(d, "d_log")), t(boot_se(d, "d_brier")))
}))
names(cv_tab) <- c("model", "log_hier", "log_indep",
                   "d_log", "d_log_se", "d_log_lo", "d_log_hi",
                   "d_brier", "d_brier_se", "d_brier_lo", "d_brier_hi")
cat("\n=== site-level cross-validation, each model against the independent fits ===\n")
print(cv_tab, row.names = FALSE, digits = 3)
cat("\nM1 should reproduce the published difference (log score -0.093,\n",
    "Brier +0.025); if it does not, the folds or the draws differ.\n", sep = "")

by_sp <- do.call(rbind, lapply(split(cmp, list(cmp$model, cmp$species)), function(d)
  if (!nrow(d)) NULL else
    data.frame(model = d$model[1], species = d$species[1],
               d_log = mean(d$d_log), d_brier = mean(d$d_brier))))
cat("\n=== by species ===\n")
print(reshape(by_sp[, c("model", "species", "d_log")], direction = "wide",
              idvar = "species", timevar = "model"), row.names = FALSE, digits = 3)

write.csv(cv_tab, file.path(OUTDIR, "ladder_cv_overall.csv"), row.names = FALSE)
write.csv(by_sp,  file.path(OUTDIR, "ladder_cv_by_species.csv"), row.names = FALSE)

}  # end if (DO_CV)

cat("\nDone. Outputs in ", OUTDIR, "\n", sep = "")
