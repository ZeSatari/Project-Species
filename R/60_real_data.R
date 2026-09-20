## =====================================================================
## 60_real_data.R -- Hubbard Brook analysis
## =====================================================================
## The previous real-data script contained twelve copy-pasted blocks. One
## of them -- CAWA -- never executed
##     data$y <- data$y[sp.names == "CAWA", , , ]
## so it operated on the full 4-D array. Two further rows of the
## published table (AMRE and BAWW) carried identical estimates and
## credible intervals despite different prevalences.
##
## Iterating over species makes both failure modes impossible: the
## species label is a loop variable, so it cannot be omitted or repeated.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)
data(hbefElev)

OUTDIR <- "results/real"
record_session(OUTDIR)

SP_NAMES <- dimnames(hbefTrends$y)[[1]]
N_SP     <- length(SP_NAMES)
N_SITE   <- dim(hbefTrends$y)[2]
N_YEAR   <- dim(hbefTrends$y)[3]
K        <- dim(hbefTrends$y)[4]
message(sprintf("hbefTrends: %d species, %d sites, %d years, %d visits.",
                N_SP, N_SITE, N_YEAR, K))

## ---------------------------------------------------------------------
## Stack one species: sites within years. ONE function, used for both
## the SSOMs and the MSOM, so the row ordering cannot diverge.
## ---------------------------------------------------------------------
stack_species <- function(sp) {
  stopifnot(sp %in% SP_NAMES)
  y <- hbefTrends$y[sp, , , ]                       # site x year x visit
  Y <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) y[, t, ]))
  stopifnot(nrow(Y) == N_SITE * N_YEAR, ncol(Y) == K)
  Y
}

stack_det_cov <- function(nm) {
  z <- hbefTrends$det.covs[[nm]]                    # site x year x visit
  M <- do.call(rbind, lapply(seq_len(N_YEAR), \(t) z[, t, ]))
  (M - mean(M, na.rm = TRUE)) / sd(M, na.rm = TRUE)
}

coords <- hbefTrends$coords / 1000                  # to km
elev   <- hbefTrends$occ.covs$elev
elev_s <- as.numeric(scale(elev))

Xocc_one <- data.frame(
  site        = rep(seq_len(N_SITE), times = N_YEAR),
  time        = rep(seq_len(N_YEAR), each  = N_SITE),
  x           = rep(coords[, 1], times = N_YEAR),
  y           = rep(coords[, 2], times = N_YEAR),
  elev_s      = rep(elev_s,      times = N_YEAR),
  elev_s2     = rep(elev_s^2,    times = N_YEAR),
  Int_occ     = 1
)
N_ROW_ONE <- nrow(Xocc_one)

Xdet_day <- stack_det_cov("day")
Xdet_tod <- stack_det_cov("tod")
stopifnot(nrow(Xdet_day) == N_ROW_ONE, nrow(Xdet_tod) == N_ROW_ONE)

## ---------------------------------------------------------------------
## Mesh and SPDE, built once.
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
## Mesh resolution.
##
## The original specification -- max.edge = c(0.1, 0.7), cutoff = 0.12 --
## produced 1724 nodes, so the latent field held 1724 * 9 = 15,516
## elements against 3,357 observation rows: 4.6 field elements per
## observation. Repeated fits of the IDENTICAL model then returned
## community intercepts spanning 0.183 on the logit scale, run times
## varying by a factor of six (1900 s to 13,100 s), and outright failures
## in 20-80% of attempts depending on the model and strategy. The
## optimiser was landing in different local optima.
##
## Coarsening to max.edge = c(0.35, 1.2), cutoff = 0.35 gives 283 nodes,
## hence 2,547 field elements and a ratio of 0.76 -- comparable to the
## simulation setting, which was stable. This is not a loss of quality:
## a mesh finer than the data can support promises a resolution that
## 373 sites do not contain, and buys it with irreproducibility.
##
## Verify stability with R/61_stability_check.R after any change here.
## ---------------------------------------------------------------------
bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)

message(sprintf("Mesh: %d nodes | field %d | obs rows %d | ratio %.2f",
                mesh$n, mesh$n * N_YEAR, N_ROW_ONE,
                mesh$n * N_YEAR / N_ROW_ONE))
if (mesh$n * N_YEAR > 2 * N_ROW_ONE) {
  warning("Latent field is more than twice the number of observations. ",
          "Expect unstable, non-reproducible fits. Coarsen the mesh.")
}

spde <- INLA::inla.spde2.pcmatern(mesh,
                                  prior.range = PRIORS$spde_real$range,
                                  prior.sigma = PRIORS$spde_real$sigma,
                                  constr = TRUE)

## =====================================================================
## PART A -- single-species framework: REVI with OVEN as observation-level
##           biotic information.
## =====================================================================
Y_revi <- stack_species("REVI")
Y_oven <- stack_species("OVEN")

stf_one <- make_st_field(mesh, spde, Xocc_one[, c("x", "y")],
                         time = Xocc_one$time, n_time = N_YEAR)

## NOTE the missing-data convention is now explicit and recorded.
## Report the proportion affected, and rerun with na_action = "keep"
## as a sensitivity check before making any claim about alpha2.
message(sprintf("Missing OVEN visit records: %.1f%%",
                100 * mean(is.na(Y_oven))))

X_M1 <- build_det_X(list(Xdet_day, Xdet_tod), K = K, na_action = "zero")
X_M2 <- build_det_X(list(Xdet_day, Xdet_tod, Y_oven), K = K,
                    na_action = "zero")

base_one <- c(as.list(Xocc_one), stf_one$idx)

f_one <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
  f(spatialfield, model = spde, A.local = stf_one$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

## Only ONE temporal process. The previous real-data models included
## f(time, model = "ar1") in addition to an AR(1)-grouped spatial field,
## i.e. two AR(1) processes competing for the same temporal signal, which
## was neither reported nor checked for identifiability. If a separate
## temporal main effect is wanted, justify it and show that the two
## components are separately estimable.

m1 <- fit_occupancy(f_one, c(base_one, list(Y = Y_revi, X = X_M1)), X_M1)
m2 <- fit_occupancy(f_one, c(base_one, list(Y = Y_revi, X = X_M2)), X_M2)

single_tab <- data.frame(
  model  = c("M1", "M2"),
  waic   = c(m1$waic$waic, m2$waic$waic),
  mlik   = c(m1$mlik[1, 1], m2$mlik[1, 1]),
  ulogcv = c(m1$ulogcv, m2$ulogcv)
)
write.csv(single_tab, file.path(OUTDIR, "single_species_comparison.csv"),
          row.names = FALSE)

## Detection coefficients. Index them by NAME, not by a hard-coded
## subscript: the manuscript reports alpha_3 in the application but
## defines alpha_2 in the Methods, with the empirical detection model
## never written down.
det_M2 <- det_draws(m2)
det_tab <- data.frame(
  coefficient = c("intercept", "day", "tod", "biotic_OVEN"),
  mean = colMeans(det_M2),
  sd   = apply(det_M2, 2, sd),
  lwr  = apply(det_M2, 2, quantile, 0.025),
  upr  = apply(det_M2, 2, quantile, 0.975)
)
write.csv(det_tab, file.path(OUTDIR, "detection_coefficients_M2.csv"),
          row.names = FALSE)
print(det_tab, digits = 3)

## =====================================================================
## PART A2 -- spatial confounding diagnostic
## =====================================================================
## Elevation is a smooth spatial surface, so it competes with the SPDE
## field for the same variation. The standard consequences are variance
## inflation and shifted point estimates for the environmental effects.
##
## Raising the issue in the Discussion without a number from this
## analysis invites the obvious question. Refitting without the field
## answers it cheaply.
##
## Reading the output:
##   means roughly unchanged, sd modestly larger  -> mild; note it
##   means shifted, or a sign change              -> serious; the
##                                                   interpretation of
##                                                   the elevation effect
##                                                   is affected
##   sd inflated several-fold                     -> report the inflation
##
## All three outcomes are reportable. Only silence is not.
## =====================================================================
f_nospat <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2

m_nospat <- fit_occupancy(f_nospat, c(base_one, list(Y = Y_revi, X = X_M1)),
                          X_M1, compute_cv = FALSE)

env_terms <- c("elev_s", "elev_s2")
conf <- data.frame(
  term          = env_terms,
  mean_with     = m1$summary.fixed[env_terms, "mean"],
  sd_with       = m1$summary.fixed[env_terms, "sd"],
  mean_without  = m_nospat$summary.fixed[env_terms, "mean"],
  sd_without    = m_nospat$summary.fixed[env_terms, "sd"],
  row.names     = NULL
)
## Shift expressed in posterior SD units of the full model: a raw
## difference of 0.02 means nothing if sd = 0.5 and a great deal if
## sd = 0.01.
conf$shift_in_sd <- abs(conf$mean_with - conf$mean_without) / conf$sd_with
conf$sd_inflation <- conf$sd_with / conf$sd_without
conf$sign_change  <- sign(conf$mean_with) != sign(conf$mean_without)

write.csv(conf, file.path(OUTDIR, "spatial_confounding.csv"), row.names = FALSE)
cat("\n-- spatial confounding: environmental effects with vs without the field --\n")
print(conf, digits = 3)
cat(sprintf("WAIC with field: %.1f | without: %.1f\n",
            m1$waic$waic, m_nospat$waic$waic))

## =====================================================================
## PART B -- hierarchical MSOM over all species.
## =====================================================================
Y_all <- do.call(rbind, lapply(SP_NAMES, stack_species))
Xocc_all <- do.call(rbind, lapply(seq_len(N_SP), function(i) {
  z <- Xocc_one; z$species_id <- i; z
}))
stopifnot(nrow(Y_all) == N_SP * N_ROW_ONE,
          nrow(Xocc_all) == nrow(Y_all))

## Detection covariates are survey-level and shared across species, but
## the design matrix must still be replicated explicitly. The previous
## script passed a 3357-row Xdet against a 40284-row Y_all and relied on
## recycling to line up.
X_msom <- build_det_X(
  list(do.call(rbind, replicate(N_SP, Xdet_day, simplify = FALSE)),
       do.call(rbind, replicate(N_SP, Xdet_tod, simplify = FALSE))),
  K = K
)
stopifnot(nrow(X_msom) == nrow(Y_all))

stf_all <- make_st_field(mesh, spde, Xocc_all[, c("x", "y")],
                         time = Xocc_all$time, n_time = N_YEAR)

dat_msom <- c(as.list(Xocc_all), list(Y = Y_all, X = X_msom), stf_all$idx)

f_msom <- inla.mdata(Y, X) ~ -1 + Int_occ + elev_s + elev_s2 +
  f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
  f(spatialfield, model = spde, A.local = stf_all$A,
    group = spatialfield.group,
    control.group = list(model = "ar1", hyper = PRIORS$ar1))

m_msom <- fit_occupancy(f_msom, dat_msom, X_msom)

## =====================================================================
## PART C -- independent SSOMs, one per species, via a loop.
## =====================================================================
## Species with very sparse detection histories can fail to converge:
## when a species is detected at almost no site-visit, psi and p are
## barely separable and the optimiser diverges. This is a substantive
## result, not merely an obstacle -- it is a direct argument for the
## hierarchical approach -- so failures are recorded rather than allowed
## to halt the run, and are reported alongside the results.
##
## First attempt uses the default starting strategy; on failure we retry
## with a Gaussian latent strategy, which is more forgiving, before
## giving up.
fit_ssom_safe <- function(sp) {
  Yi <- stack_species(sp)
  dat_i <- c(base_one, list(Y = Yi, X = X_M1))

  out <- tryCatch(
    fit_occupancy(f_one, dat_i, X_M1),
    error = function(e) e
  )
  if (!inherits(out, "error")) return(out)

  message("  ", sp, " failed; retrying with strategy = 'gaussian'")
  out2 <- tryCatch(
    fit_occupancy(f_one, dat_i, X_M1, strategy = "gaussian"),
    error = function(e) e
  )
  if (!inherits(out2, "error")) {
    attr(out2, "fallback") <- "gaussian"
    return(out2)
  }

  message("  ", sp, " did not converge under either strategy.")
  structure(list(failed = TRUE,
                 message = conditionMessage(out2),
                 naive_occ = naive_occupancy(Yi),
                 det_rate  = detection_rate(Yi)),
            class = "ssom_failed")
}

ssom_fits <- lapply(SP_NAMES, function(sp) {
  message("SSOM: ", sp)
  fit_ssom_safe(sp)
})
names(ssom_fits) <- SP_NAMES

ssom_ok <- !vapply(ssom_fits, inherits, logical(1), "ssom_failed")
if (any(!ssom_ok)) {
  message("\nSingle-species models that did not converge: ",
          paste(SP_NAMES[!ssom_ok], collapse = ", "))
  message("Report these: a species that cannot be fitted independently ",
          "but is estimable within the community model is evidence for ",
          "hierarchical pooling.")
}
saveRDS(data.frame(species = SP_NAMES, converged = ssom_ok),
        file.path(OUTDIR, "ssom_convergence.rds"))

## ---------------------------------------------------------------------
## Species table. Reports the IDENTIFIED quantity beta0 + gamma_i, naive
## occupancy (defined), and the like-for-like predictive difference.
## Normalised WAIC is deliberately absent: MSOM WAIC / 40284 and SSOM
## WAIC / 3357 are computed on different response sets and dividing by n
## does not make them comparable.
## ---------------------------------------------------------------------
sp_int <- species_intercepts(m_msom, N_SP)
sp_int$species_name <- SP_NAMES

species_tab <- do.call(rbind, lapply(seq_len(N_SP), function(i) {
  Yi  <- stack_species(SP_NAMES[i])
  ok  <- ssom_ok[i]
  data.frame(
    species      = SP_NAMES[i],
    ssom_converged = ok,
    naive_occ    = naive_occupancy(Yi),
    det_rate     = detection_rate(Yi),
    ssom_int     = if (ok) ssom_fits[[i]]$summary.fixed["Int_occ", "mean"]
                   else NA_real_,
    msom_int     = sp_int$mean[i],
    msom_lwr     = sp_int$lwr[i],
    msom_upr     = sp_int$upr[i]
  )
}))

## Shrinkage, computed rather than asserted.
b0 <- m_msom$summary.fixed["Int_occ", "mean"]
species_tab$dist_ssom <- abs(species_tab$ssom_int - b0)
species_tab$dist_msom <- abs(species_tab$msom_int - b0)
species_tab$shrinkage <- species_tab$dist_ssom - species_tab$dist_msom
species_tab$shrunk_toward_mean <- species_tab$shrinkage > 0

## Sanity check. The species-level intercept beta0 + gamma_i must track
## naive occupancy closely; anything else means the posterior components
## have been mis-indexed. Report this correlation -- it tells a reviewer
## the community model recovers something sensible.
r_naive <- cor(species_tab$naive_occ, species_tab$msom_int)
if (r_naive < 0.8) {
  warning(sprintf("Species intercepts correlate only %.2f with naive occupancy. Check indexing.",
                  r_naive))
}
message(sprintf("Correlation of beta0 + gamma_i with naive occupancy: %.3f",
                r_naive))

write.csv(species_tab, file.path(OUTDIR, "Table1_species.csv"),
          row.names = FALSE)
print(species_tab, digits = 3)

message("\nCommunity intercept: ",
        sprintf("%.3f (%.3f, %.3f)",
                b0,
                m_msom$summary.fixed["Int_occ", "0.025quant"],
                m_msom$summary.fixed["Int_occ", "0.975quant"]))
message("Species shrunk toward the community mean: ",
        sum(species_tab$shrunk_toward_mean), " / ", N_SP)

saveRDS(list(m1 = m1, m2 = m2, msom = m_msom, ssom = ssom_fits),
        file.path(OUTDIR, "fits.rds"))
