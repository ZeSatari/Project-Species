## =====================================================================
## 40_run_single_sim.R -- replicated single-species simulation
## =====================================================================
## Replaces the previous single-fit analysis. Every number that reaches
## Table S1 / S2 comes from this driver.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")
source("R/30_fit_occupancy.R")

## Set RUN_FULL <- FALSE before sourcing to load the configuration and
## the run_one_rep() definition WITHOUT launching the full study. Use
## that to time a single replicate first:
##     RUN_FULL <- FALSE
##     source("R/40_run_single_sim.R")
##     future::plan(future::sequential)
##     system.time(r1 <- run_one_rep(1, "independent"))
if (!exists("RUN_FULL")) RUN_FULL <- TRUE

OUTDIR <- "results/single"
record_session(OUTDIR)
check_crps_implementation()

N_REP     <- 200          # >= 200. With R = 1 nothing below is defined.
SEED_BASE <- 20260101L
SCENARIOS <- c("independent", "shared", "observer")

## ---------------------------------------------------------------------
## PROBLEM SIZE.
##
## Timing on the full-size problem (2000 sites, mesh 2217 nodes, nT = 5)
## was roughly 600 s per fit, i.e. ~33 h for 1200 fits on six workers.
## The settings below shrink the problem to make R = 200 feasible.
##
## This is the right trade-off. The number of sites and the mesh
## resolution affect the STATISTICAL POWER of the simulation study, not
## its validity: bias, RMSE and coverage remain well defined at any
## sample size. R = 1 makes them undefined. So shrink the problem, never
## the replication.
##
## Report the simulation dimensions in Appendix S1 exactly as set here.
## The approximation check (R/70_strategy_check.R) was run at FULL size,
## so its conclusion still covers the scale used in the application.
## ---------------------------------------------------------------------
SIM <- list(
  site_frac = 0.05,             # 500 of 10000 cells (was 0.20)
  max_edge  = c(25, 50),        # coarser mesh (was c(15, 30))
  strategy  = "simplified.laplace"
)
N_WORKERS <- 6                  # set by hand; see the note below

## Cross-validation on every replicate. Profiling on this problem size:
##   fit                86 s
##   inla.group.cv       6 s
##   posterior sampling 66 s   <- the real cost, now removed
## CV is cheap, so there is nothing to gain by skipping it. The former
## sampling cost is avoided by summarising from the stored marginals
## (summarise_marginal) instead of reconstructing the latent field.
CV_EVERY <- 1

DOM <- make_domain(max_edge = SIM$max_edge)
message("SIM config -- site_frac: ", SIM$site_frac,
        " | max.edge: ", paste(SIM$max_edge, collapse = ","),
        " | strategy: ", SIM$strategy)

## Time one fit before committing to the full run.
##   future::plan(future::sequential)
##   system.time(run_one_rep(1, "independent"))
## Target roughly 30-60 s. If it is much higher, reduce site_frac again
## or coarsen max_edge further.

## ---------------------------------------------------------------------
## One replicate: generate, fit M1 and M2, return tidy per-parameter rows
## plus the model-comparison row.
## ---------------------------------------------------------------------
run_one_rep <- function(r, scenario) {

  init_worker()

  sim <- simulate_single(DOM, scenario = scenario,
                         site_frac = SIM$site_frac, seed = SEED_BASE + r)
  d   <- sim$data
  K   <- sim$truth$K

  spde <- INLA::inla.spde2.pcmatern(
    DOM$mesh,
    prior.range = PRIORS$spde_sim$range,
    prior.sigma = PRIORS$spde_sim$sigma
  )
  stf <- make_st_field(DOM$mesh, spde, d[, c("x.loc", "y.loc")],
                       time = d$time, n_time = sim$truth$nT)

  Y  <- sim$Y_A
  gA <- as.matrix(d[, paste0("gA", seq_len(K))])
  yB <- as.matrix(d[, paste0("yB", seq_len(K))])

  X1 <- build_det_X(list(gA),     K = K)               # M1
  X2 <- build_det_X(list(gA, yB), K = K, na_action = "zero")  # M2

  base <- c(list(
    Int_occ = rep(1, nrow(d)),
    x_sq    = d$x_A^2            # precomputed: NEVER `x^2` inside a formula
  ), stf$idx)

  ## Identical occupancy component in both models, by construction.
  f <- inla.mdata(Y, X) ~ -1 + Int_occ + x_sq +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))

  do_cv <- (r %% CV_EVERY == 1L) || (CV_EVERY == 1L)
  m1 <- fit_occupancy(f, c(base, list(Y = Y, X = X1)), X1,
                      strategy = SIM$strategy, compute_cv = do_cv)
  m2 <- fit_occupancy(f, c(base, list(Y = Y, X = X2)), X2,
                      strategy = SIM$strategy, compute_cv = do_cv)

  ## ---- parameter recovery -------------------------------------------
  tr <- sim$truth
  collect <- function(fit, label) {
    fx <- fixed_marginals(fit, c("Int_occ", "x_sq"))
    dt <- det_marginals(fit)
    rows <- list(
      cbind(parameter = "beta0",  summarise_marginal(fx$Int_occ, tr$beta0)),
      cbind(parameter = "beta1",  summarise_marginal(fx$x_sq,    tr$beta1)),
      cbind(parameter = "alpha0", summarise_marginal(dt$alpha0,  tr$alpha0)),
      cbind(parameter = "alpha1", summarise_marginal(dt$alpha1,  tr$alpha1))
    )
    if (length(dt) >= 3 && !is.na(tr$alpha2)) {
      rows <- c(rows, list(
        cbind(parameter = "alpha2", summarise_marginal(dt$alpha2, tr$alpha2))
      ))
    }
    cbind(rep = r, scenario = scenario, model = label, do.call(rbind, rows))
  }

  params <- rbind(collect(m1, "M1"), collect(m2, "M2"))

  nn <- function(x) if (is.null(x)) NA_real_ else x
  fitstats <- data.frame(
    rep = r, scenario = scenario,
    waic_M1 = m1$waic$waic,   waic_M2 = m2$waic$waic,
    ## unname(): mlik[1,1] carries the name "log marginal-likelihood
    ## (integration)", which otherwise becomes the data.frame row name
    ## and leaks into the CSV.
    mlik_M1 = unname(m1$mlik[1, 1]),  mlik_M2 = unname(m2$mlik[1, 1]),
    ulogcv_M1 = nn(m1$ulogcv), ulogcv_M2 = nn(m2$ulogcv),
    secs_M1 = m1$elapsed,     secs_M2 = m2$elapsed,
    naive_occ = naive_occupancy(Y),
    det_rate  = detection_rate(Y),
    row.names = NULL
  )

  list(params = params, fitstats = fitstats)
}

## ---------------------------------------------------------------------
## Run. For a Slurm cluster swap the plan for future.batchtools:
##   future::plan(future.batchtools::batchtools_slurm,
##                resources = list(ncpus = 1, memory = 8000,
##                                 walltime = 3600))
## ---------------------------------------------------------------------
if (!RUN_FULL) {
  message("RUN_FULL = FALSE: configuration and run_one_rep() loaded; ",
          "the full study was NOT started.")
} else {

## Do NOT use detectCores() on a hybrid P-core/E-core CPU: workers landing
## on efficiency cores run several times slower, and future_lapply waits
## for the slowest. Each worker is also a full R process holding its own
## copy of the mesh and projector, so memory binds before cores do.
future::plan(future::multisession, workers = N_WORKERS)

all_params <- all_fits <- list()

for (sc in SCENARIOS) {
  message("=== scenario: ", sc, " ===")
  res <- future_lapply(
    seq_len(N_REP),
    \(r) tryCatch(run_one_rep(r, sc),
                  error = function(e) {
                    warning("rep ", r, " (", sc, ") failed: ",
                            conditionMessage(e))
                    NULL
                  }),
    future.seed = TRUE,

    future.globals = TRUE,

    future.packages = c("INLA", "fmesher", "sf", "terra", "dplyr", "MASS")
  )
  res <- Filter(Negate(is.null), res)
  message(sprintf("  %d / %d replicates completed.", length(res), N_REP))

  all_params[[sc]] <- do.call(rbind, lapply(res, `[[`, "params"))
  all_fits[[sc]]   <- do.call(rbind, lapply(res, `[[`, "fitstats"))
}

params_df <- do.call(rbind, all_params)
fits_df   <- do.call(rbind, all_fits)

saveRDS(params_df, file.path(OUTDIR, "per_replicate_params.rds"))
saveRDS(fits_df,   file.path(OUTDIR, "per_replicate_fitstats.rds"))

## ---------------------------------------------------------------------
## Table S2 replacement: bias / RMSE / coverage / CRPS by scenario+model.
## Read the `coverage` column first. If M2 shows coverage well below 0.95
## for beta0 under "shared" but near 0.95 under "independent", that is
## the paper's result: the method is valid exactly when Assumption 1
## holds. Say so plainly rather than describing it as a bias--variance
## trade-off.
## ---------------------------------------------------------------------
summary_tab <- params_df |>
  dplyr::group_by(scenario, model) |>
  dplyr::group_modify(~ aggregate_reps(.x)) |>
  dplyr::ungroup()

write.csv(summary_tab, file.path(OUTDIR, "TableS2_recovery.csv"),
          row.names = FALSE)
print(as.data.frame(summary_tab), digits = 3)

## Model comparison, with Monte Carlo error -- a mean WAIC difference is
## only interpretable alongside its variability across replicates.
cmp <- fits_df |>
  dplyr::group_by(scenario) |>
  dplyr::summarise(
    n_rep          = dplyr::n(),
    d_waic_mean    = mean(waic_M2 - waic_M1),
    d_waic_se      = sd(waic_M2 - waic_M1) / sqrt(dplyr::n()),
    n_cv           = sum(!is.na(ulogcv_M1)),
    d_ulogcv_mean  = mean(ulogcv_M2 - ulogcv_M1, na.rm = TRUE),
    d_ulogcv_se    = sd(ulogcv_M2 - ulogcv_M1, na.rm = TRUE) /
                       sqrt(sum(!is.na(ulogcv_M1))),
    mean_secs      = mean(secs_M1 + secs_M2),
    prop_M2_better = mean(waic_M2 < waic_M1),
    .groups = "drop"
  )
write.csv(cmp, file.path(OUTDIR, "TableS1_comparison.csv"), row.names = FALSE)
print(as.data.frame(cmp), digits = 3)

}  # end if (RUN_FULL)
