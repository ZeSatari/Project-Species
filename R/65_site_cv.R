## =====================================================================
## 65_site_cv.R -- site-held-out cross-validation, hierarchical vs
##                 independent
## =====================================================================
## The manuscript declines to compare the two model classes by
## inla.group.cv, because automatic group construction forms groups from
## the model's own dependence structure: under the hierarchical model the
## groups reach across species, so the two are not scored on the same
## withholding. The manuscript says the valid comparison would use folds
## defined by the analyst. This script does that.
##
## DESIGN. Sites are partitioned into K folds. For each fold, every
## observation at those sites is withheld -- for all twelve species at
## once -- and both model classes are refitted on what remains. Held-out
## detections are then predicted from the covariates and from the latent
## field interpolated to the held-out locations. Both classes see exactly
## the same folds and the same withheld records.
##
## WHAT THE COMPARISON DOES AND DOES NOT ISOLATE. Each independent fit
## sees one species; the hierarchical fit sees twelve. If the
## hierarchical model predicts better, part of that advantage comes from
## having more data with which to estimate the shared field, not from
## pooling the species intercepts. The two cannot be separated here, and
## the result should be read as a comparison of model classes as they
## would actually be used, not as an isolation of the pooling mechanism.
##
## SCORES. Two, because they answer different questions:
##   log score  mean log predictive density of the held-out detection
##              histories, marginalising over the latent occupancy state.
##              The quantity a statistician wants; hard to interpret.
##   Brier      mean squared error of P(detected at least once) against
##              what was observed. Interpretable on the probability
##              scale; less sensitive to confident errors.
## AUC is reported alongside as a rank measure.
##
## Runtime: K x 13 fits. At K = 5 that is 65 fits, several hours.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/cv"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR)

K_FOLDS <- 5
SEED    <- 20260401L
N_DRAWS <- 500

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP   <- length(SP_NAMES)
N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
K_VIS  <- dim(hbefTrends$y)[4]

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
elev_s <- as.numeric(scale(hbefTrends$occ.covs$elev))
site_id <- rep(seq_len(N_SITE), times = N_YEAR)
time_id <- rep(seq_len(N_YEAR), each = N_SITE)
N_ROW   <- N_SITE * N_YEAR

Xday <- stack_det_cov("day"); Xtod <- stack_det_cov("tod")
Y_all <- lapply(SP_NAMES, stack_species); names(Y_all) <- SP_NAMES

## Folds are over SITES, not rows: a site is in exactly one fold, and all
## of its site-years across all species go with it.
set.seed(SEED)
fold_of_site <- sample(rep_len(seq_len(K_FOLDS), N_SITE))
message("fold sizes (sites): ", paste(table(fold_of_site), collapse = " "))

mesh <- fmesher::fm_mesh_2d(
  boundary = fmesher::fm_nonconvex_hull(coords, convex = 0.3),
  max.edge = c(0.35, 1.2), min.angle = 21,
  offset = c(0.05, 1.5), cutoff = 0.35)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

## ---------------------------------------------------------------------
## Scores. Each is computed per held-out site-year and averaged.
##
## The log score marginalises over occupancy, so it is the density of the
## observed detection history under the fitted model:
##   psi * prod_j p^y (1-p)^(1-y)  +  (1-psi) * I{no detection}
## ---------------------------------------------------------------------
log_score <- function(psi, p, yrow) {
  obs <- !is.na(yrow)
  if (!any(obs)) return(NA_real_)
  lp <- sum(yrow[obs] * log(p[obs]) + (1 - yrow[obs]) * log(1 - p[obs]))
  if (any(yrow[obs] == 1)) log(psi) + lp
  else log(psi * exp(lp) + (1 - psi))
}
p_recorded <- function(psi, p, yrow) {
  obs <- !is.na(yrow)
  psi * (1 - prod(1 - p[obs]))
}

auc <- function(score, label) {
  if (length(unique(label)) < 2) return(NA_real_)
  r <- rank(score)
  n1 <- sum(label == 1); n0 <- sum(label == 0)
  (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---------------------------------------------------------------------
## Fit one model and predict at the held-out sites.
##
## Withholding is by setting the response to NA: INLA then treats those
## rows as unobserved but keeps them in the latent field, so the field is
## interpolated to the held-out locations from their neighbours. This is
## what makes site-level prediction possible at all.
## ---------------------------------------------------------------------
predict_fold <- function(fit, mesh, held_rows, species_index = NULL) {
  s  <- INLA::inla.posterior.sample(n = N_DRAWS, result = fit)
  nm <- rownames(s[[1]]$latent)
  i1 <- function(pat) { i <- grep(pat, nm); stopifnot(length(i) == 1); i }
  i_int <- i1("^Int_occ:"); i_e1 <- i1("^elev_s:"); i_e2 <- i1("^elev_s2:")
  i_sp  <- if (is.null(species_index)) NA_integer_ else
             i1(sprintf("^species_id:%d$", species_index))
  i_f   <- grep("^spatialfield:", nm)
  n_spde <- length(i_f) / N_YEAR

  if (max(held_rows) > N_ROW)
    stop("held_rows must index the single-species layout (max ", N_ROW,
         "); received a maximum of ", max(held_rows))

  psi <- matrix(NA_real_, length(held_rows), N_DRAWS)
  for (d in seq_len(N_DRAWS)) {
    v <- s[[d]]$latent
    fld <- numeric(N_ROW)
    for (tt in seq_len(N_YEAR)) {
      sl <- i_f[((tt - 1) * n_spde + 1):(tt * n_spde)]
      fld[time_id == tt] <- as.numeric(A_site %*% v[sl])
    }
    eta <- v[i_int] + v[i_e1] * elev_s[site_id] +
           v[i_e2] * elev_s[site_id]^2 + fld
    if (!is.na(i_sp)) eta <- eta + v[i_sp]
    psi[, d] <- plogis(eta[held_rows])
  }
  dr <- det_draws(fit, n = N_DRAWS)
  list(psi = rowMeans(psi), det = colMeans(dr))
}

## ---------------------------------------------------------------------
## Cross-validation
## ---------------------------------------------------------------------
rows <- list()

for (k in seq_len(K_FOLDS)) {
  held_sites <- which(fold_of_site == k)
  held_rows  <- which(site_id %in% held_sites)
  message(sprintf("fold %d: %d sites, %d site-years withheld",
                  k, length(held_sites), length(held_rows)))

  ## ---- hierarchical: one fit, all species ----
  Y_h <- do.call(rbind, lapply(SP_NAMES, function(sp) {
    y <- Y_all[[sp]]; y[held_rows, ] <- NA; y
  }))
  Xd_h <- build_det_X(list(do.call(rbind, rep(list(Xday), N_SP)),
                           do.call(rbind, rep(list(Xtod), N_SP))), K = K_VIS)
  spde <- INLA::inla.spde2.pcmatern(mesh,
            prior.range = PRIORS$spde_real$range,
            prior.sigma = PRIORS$spde_real$sigma, constr = TRUE)
  Xocc_h <- data.frame(
    Int_occ = 1,
    elev_s  = rep(elev_s[site_id], times = N_SP),
    elev_s2 = rep(elev_s[site_id]^2, times = N_SP),
    species_id = rep(seq_len(N_SP), each = N_ROW),
    x = rep(coords[site_id, 1], times = N_SP),
    y = rep(coords[site_id, 2], times = N_SP),
    time = rep(time_id, times = N_SP))
  stf_h <- make_st_field(mesh, spde, Xocc_h[, c("x", "y")],
                         time = Xocc_h$time, n_time = N_YEAR)
  f_h <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
    f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(spatialfield, model = spde, A.local = stf_h$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))
  fit_h <- fit_occupancy(f_h, c(as.list(Xocc_h), stf_h$idx,
                                list(Y = Y_h, X = Xd_h)),
                         Xd_h, compute_cv = FALSE)

  ## ---- independent: twelve fits ----
  for (i in seq_len(N_SP)) {
    sp <- SP_NAMES[i]
    Y_s <- Y_all[[sp]]; Y_s[held_rows, ] <- NA
    Xd_s <- build_det_X(list(Xday, Xtod), K = K_VIS)
    stf_s <- make_st_field(mesh, spde, data.frame(x = coords[site_id, 1],
                                                  y = coords[site_id, 2]),
                           time = time_id, n_time = N_YEAR)
    dat_s <- c(list(Int_occ = rep(1, N_ROW),
                    elev_s = elev_s[site_id], elev_s2 = elev_s[site_id]^2,
                    Y = Y_s, X = Xd_s), stf_s$idx)
    f_s <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
      f(spatialfield, model = spde, A.local = stf_s$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))
    fit_s <- tryCatch(fit_occupancy(f_s, dat_s, Xd_s, compute_cv = FALSE),
                      error = function(e) NULL)
    if (is.null(fit_s)) { message("  ", sp, ": independent fit failed"); next }

    ## held_rows indexes the single-species layout, which is what
    ## predict_fold works in: the field, elevation and site are shared
    ## across species, and the only species-specific term is v[i_sp],
    ## supplied through `species_index`. An earlier version offset these
    ## indices by (i-1)*N_ROW to reach into the stacked twelve-species
    ## design; that ran off the end of a length-N_ROW vector and returned
    ## NA for every species after the first.
    ph <- predict_fold(fit_h, mesh, held_rows, i)
    ps <- predict_fold(fit_s, mesh, held_rows, NULL)

    Yh <- Y_all[[sp]][held_rows, , drop = FALSE]
    for (r in seq_along(held_rows)) {
      yr <- Yh[r, ]
      if (all(is.na(yr))) next
      p_h <- plogis(ph$det[1] + ph$det[2] * Xday[held_rows[r], ] +
                      ph$det[3] * Xtod[held_rows[r], ])
      p_s <- plogis(ps$det[1] + ps$det[2] * Xday[held_rows[r], ] +
                      ps$det[3] * Xtod[held_rows[r], ])
      rows[[length(rows) + 1]] <- data.frame(
        fold = k, species = sp, row = held_rows[r],
        observed = as.integer(any(yr == 1, na.rm = TRUE)),
        log_h = log_score(ph$psi[r], p_h, yr),
        log_s = log_score(ps$psi[r], p_s, yr),
        prec_h = p_recorded(ph$psi[r], p_h, yr),
        prec_s = p_recorded(ps$psi[r], p_s, yr))
    }
  }
}

cv <- do.call(rbind, rows)
saveRDS(cv, file.path(OUTDIR, "site_cv_raw.rds"))

## If the species effect never reached the predictions, every species
## would carry the same predicted probability. Check before scoring.
sp_mean <- tapply(cv$prec_h, cv$species, mean, na.rm = TRUE)
if (diff(range(sp_mean, na.rm = TRUE)) < 0.01) {
  stop("hierarchical predictions barely differ across species (range ",
       signif(diff(range(sp_mean, na.rm = TRUE)), 3),
       "): the species effect is not reaching the prediction.")
}
if (mean(is.na(cv$prec_h)) > 0.01) {
  stop(sprintf("%.0f%% of hierarchical predictions are NA; check indexing.",
               100 * mean(is.na(cv$prec_h))))
}

## ---------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------
cv <- cv[!is.na(cv$log_h) & !is.na(cv$log_s), ]

overall <- data.frame(
  model = c("hierarchical", "independent"),
  log_score = c(mean(cv$log_h), mean(cv$log_s)),
  brier     = c(mean((cv$prec_h - cv$observed)^2),
                mean((cv$prec_s - cv$observed)^2)),
  auc       = c(auc(cv$prec_h, cv$observed), auc(cv$prec_s, cv$observed))
)
cat("\n", strrep("=", 70), "\n", sep = "")
print(overall, row.names = FALSE, digits = 4)

## Paired, since both models score the same held-out records.
d_log <- cv$log_h - cv$log_s
d_bri <- (cv$prec_h - cv$observed)^2 - (cv$prec_s - cv$observed)^2
cat(sprintf("\npaired difference, log score : %+.4f (se %.4f)\n",
            mean(d_log), sd(d_log) / sqrt(nrow(cv))))
cat(sprintf("paired difference, Brier     : %+.5f (se %.5f)\n",
            mean(d_bri), sd(d_bri) / sqrt(nrow(cv))))
cat("  negative Brier difference favours the hierarchical model;\n")
cat("  positive log-score difference favours it.\n")

## By species, since the claim concerns sparsely recorded ones.
by_sp <- do.call(rbind, lapply(split(cv, cv$species), function(d) data.frame(
  species = d$species[1], n = nrow(d), naive = mean(d$observed),
  d_log = mean(d$log_h - d$log_s),
  d_brier = mean((d$prec_h - d$observed)^2) - mean((d$prec_s - d$observed)^2)
)))
by_sp <- by_sp[order(by_sp$naive), ]
cat("\nby species, ordered by observed frequency:\n")
print(by_sp, row.names = FALSE, digits = 4)

cat("\nIf the hierarchical model helps where information is scarce, the\n")
cat("advantage should be largest at the top of this table and absent or\n")
cat("reversed at the bottom. A uniform difference would instead suggest the\n")
cat("shared field, not pooling, is doing the work.\n")

write.csv(overall, file.path(OUTDIR, "site_cv_overall.csv"), row.names = FALSE)
write.csv(by_sp,   file.path(OUTDIR, "site_cv_by_species.csv"), row.names = FALSE)
