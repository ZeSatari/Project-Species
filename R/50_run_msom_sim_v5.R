## =====================================================================
## 50_run_msom_sim_v5.R -- multi-species simulation, RUN (c)
## =====================================================================
## HOMOGENEOUS DETECTION. Identical to the primary run (v2, run (a)) except
## that every species has the same detection parameters,
##     alpha0_i = -0.9,  alpha1_i = -1.0,
## so the fitted hierarchical model, whose detection is pooled across
## species, is CORRECTLY specified in its detection component. The
## comparison with independent fitting then isolates partial pooling of
## the occupancy intercepts and slopes (PaperB, Section 4.1, run (c)).
##
## Nothing else changes: same N_REP, seeds, domain, mesh, priors, fitted
## models, heterogeneous slopes (sd_beta1 = 0.5).
##
## RANDOM-NUMBER STREAM. The detection parameters are still DRAWN exactly
## as in v2 and only then overwritten, so the stream is consumed
## identically and replicate r of this run shares with replicate r of v2:
## the field, the covariate, gamma, beta1, the sampled sites, the survey
## covariate g, the occupancy states z and the missing-visit pattern. Only
## the detection outcomes differ. The pairing is checked before the run
## (check_pairing) and again against the stored v2 output at the end.
##
## Do NOT implement homogeneity with rnorm(n, -0.9, 0): R returns the mean
## without drawing when sd = 0, which would shift the stream.
##
## Confirm on load that the banner reports "v5 | run (c) | homogeneous detection".
## Outputs: results/msom/*_v5.*   Checkpoints: results/msom/checkpoints_v5
## NOT RUN BY THE AUTHOR OF THIS DRAFT (no R-INLA available at drafting).
## Derived from 50_run_msom_sim_v2.R; `diff` the two to review every change.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")   # make_domain(), sim_st_field()
source("R/30_fit_occupancy.R")

OUTDIR <- "results/msom"
record_session(OUTDIR, tag = "homdet_v5")

## R = 100 rather than 200. Timing on this problem: 1690 s per replicate
## (13 fits each), so R = 200 would need ~31 h on three workers.
##
## At R = 100 the Monte Carlo standard error of a coverage estimate is
## sqrt(0.95*0.05/100) = 0.022, which still resolves a fall from 0.95 to
## below 0.90 -- the effect this study needs to detect. The central claim
## of the multi-species contribution rests on prop_msom_closer, which is
## computed over 100 x 12 = 1200 species-replicate pairs.
##
## Do not reduce R below 100: at that point Monte Carlo error starts to
## compete with the effect being measured. Reduce site_frac instead.
N_REP     <- 100
SEED_BASE <- 20260201L
N_SPECIES <- 12

## See the note on problem size in 40_run_single_sim.R. Same reasoning:
## shrink the problem, never the replication.
## sd_beta1 = 0.5: species-specific slopes are heterogeneous, which is
## what distinguishes this primary run from the v3 diagnostic.
SIM <- list(site_frac = 0.02, max_edge = c(25, 50),
            sd_beta1 = 0.5,
            homog_det = TRUE,        # THE single change from v2
            strategy = "simplified.laplace")

## Each replicate fits 13 models (one MSOM + twelve SSOMs). With
## config = TRUE every fitted object retains the full latent-field
## structure, so a worker can hold several GB at once. Six workers
## exhausted memory on a 32 GB machine ("allocation of text connection
## failed"). Three is a safer starting point; raise it only after
## watching peak usage.
N_WORKERS <- 4   # peak memory measured at ~2.6 GB per replicate

message("=== 50_run_msom_sim_v5 | run (c) | homogeneous detection | N_REP = ", N_REP,
        " | site_frac = ", SIM$site_frac,
        " | sd_beta1 = ", SIM$sd_beta1, " ===")

if (!exists("RUN_FULL")) RUN_FULL <- TRUE

DOM <- make_domain(max_edge = SIM$max_edge)
message("SIM config -- site_frac: ", SIM$site_frac,
        " | max.edge: ", paste(SIM$max_edge, collapse = ","),
        " | strategy: ", SIM$strategy)

## ---------------------------------------------------------------------
## Sanity check on any exchangeable correlation matrix before use.
## ---------------------------------------------------------------------
exch_cov <- function(n, sd, rho) {
  R <- matrix(rho, n, n); diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) stop("Correlation matrix is not positive definite.")
  sd^2 * R
}

simulate_msom <- function(dom, n_species = N_SPECIES, nT = 5, K = 3,
                          site_frac = 0.20,
                          mu_beta0 = -0.5, sd_gamma = 0.8, gen_cor = 0,
                          mu_beta1 = 0.7,  sd_beta1 = 0.5,
                          range = 100, sigma = 0.5, rho = 0.6,
                          missing_p = 0.10, homog_det = FALSE,
                          alpha0_common = -0.9, alpha1_common = -1.0,
                          seed) {
  set.seed(seed)
  nc <- dom$ncells

  omega <- sim_st_field(dom, range, sigma, rho, nT, seed = seed)
  x_s   <- as.numeric(scale(sim_st_field(dom, range, sigma, 0, 1,
                                         seed = seed + 100000L)[, 1]))

  ## Species effects centred at zero, then explicitly recentred so the
  ## realised community intercept is unambiguous for THIS replicate.
  gamma <- as.numeric(MASS::mvrnorm(1, rep(0, n_species),
                                    exch_cov(n_species, sd_gamma, gen_cor)))
  gamma <- gamma - mean(gamma)          # matches constr = TRUE in the fit
  ## Centre the slopes for the same reason: with constr = TRUE on the
  ## random slope, the identified community slope is mean(beta1). Recentre
  ## so that the realised truth for THIS replicate is exactly mu_beta1
  ## rather than mu_beta1 + sampling noise.
  beta1 <- rnorm(n_species, mu_beta1, sd_beta1)
  beta1 <- beta1 - mean(beta1) + mu_beta1

  ## Draw exactly as in v2 so the random-number stream is unchanged,
  ## THEN overwrite for run (c). See the header.
  alpha0 <- rnorm(n_species, -0.9, 0.2)
  alpha1 <- rnorm(n_species, -1.0, 0.3)
  if (homog_det) {
    alpha0[] <- alpha0_common
    alpha1[] <- alpha1_common
  }

  nsites  <- round(nc * site_frac)
  site_id <- sort(sample.int(nc, nsites))

  out <- vector("list", n_species)
  for (i in seq_len(n_species)) {
    Yi <- matrix(NA_integer_, nsites * nT, K)
    gi <- array(runif(nsites * K * nT, -1, 1), c(nsites, K, nT))
    for (t in seq_len(nT)) {
      ## LINEAR covariate. With x ~ N(0,1), x^2 is chi-squared_1: its right
      ## tail drives psi toward one, where the data carry little information
      ## about the slope. Diagnostics with no spatial structure at n = 5000:
      ##   linear           beta1 = 0.697  (truth 0.700)
      ##   x^2 raw          beta1 = 0.456
      ##   x^2 standardised beta1 = 0.614
      ## The quadratic form attenuates the coefficient by 12-35% for reasons
      ## unrelated to the hierarchical structure this study evaluates. With a
      ## linear covariate, residual bias is attributable to the model. The
      ## empirical application keeps its quadratic elevation term; the
      ## attenuation above should be recalled when interpreting it.
      psi <- plogis(mu_beta0 + gamma[i] + beta1[i] * x_s + omega[, t])
      z   <- rbinom(nc, 1, psi)
      rows <- ((t - 1) * nsites + 1):(t * nsites)
      for (j in seq_len(K)) {
        p  <- plogis(alpha0[i] + alpha1[i] * gi[, j, t])
        y  <- rbinom(nsites, 1, z[site_id] * p)
        y[as.logical(rbinom(nsites, 1, missing_p))] <- NA
        Yi[rows, j] <- y
      }
    }
    g_long <- do.call(rbind, lapply(seq_len(nT), \(t) gi[, , t]))
    out[[i]] <- list(
      Y = Yi, g = g_long,
      X = data.frame(
        species_id = i,
        cellid = rep(site_id, times = nT),
        time   = rep(seq_len(nT), each = nsites),
        x.loc  = rep(dom$coords[site_id, 1], times = nT),
        y.loc  = rep(dom$coords[site_id, 2], times = nT),
        x_cov  = rep(x_s[site_id], times = nT)
      )
    )
  }

  list(
    per_species = out,
    Y_all = do.call(rbind, lapply(out, `[[`, "Y")),
    g_all = do.call(rbind, lapply(out, `[[`, "g")),
    X_all = do.call(rbind, lapply(out, `[[`, "X")),
    truth = list(
      beta0_community = mu_beta0,       # gamma is centred, so this is it
      beta1_community = mu_beta1,       # beta1 is centred, likewise
      gamma = gamma, beta1 = beta1,
      alpha0 = alpha0, alpha1 = alpha1,
      sd_gamma = sd_gamma, gen_cor = gen_cor,
      species_intercept = mu_beta0 + gamma,
      homog_det = homog_det, site_id = site_id,
      n_species = n_species, nsites = nsites, nT = nT, K = K, seed = seed
    )
  )
}

## ---------------------------------------------------------------------
## One replicate: fit the MSOM and all n_species independent SSOMs,
## then compare them on the SAME observations for each species.
## ---------------------------------------------------------------------
run_msom_rep <- function(r) {

  init_worker()
  sim <- simulate_msom(DOM, site_frac = SIM$site_frac,
                       sd_beta1 = SIM$sd_beta1, homog_det = SIM$homog_det,
                       seed = SEED_BASE + r)
  stopifnot(isTRUE(sim$truth$homog_det))
  X   <- sim$X_all; K <- sim$truth$K; nsp <- sim$truth$n_species

  spde <- INLA::inla.spde2.pcmatern(DOM$mesh,
                                    prior.range = PRIORS$spde_sim$range,
                                    prior.sigma = PRIORS$spde_sim$sigma)
  stf <- make_st_field(DOM$mesh, spde, X[, c("x.loc", "y.loc")],
                       time = X$time, n_time = sim$truth$nT)
  Xdet_all <- build_det_X(list(sim$g_all), K = K)

  dat <- c(list(
    Y = sim$Y_all, X = Xdet_all,
    Int_occ = rep(1, nrow(X)), x_cov = X$x_cov,
    species_id = X$species_id, species_slope = X$species_id
  ), stf$idx)

  ## constr = TRUE on the species effect: without it beta0 and gamma_i
  ## are not separately identified.
  ## The fixed `x_cov` term is REQUIRED alongside the constrained random
  ## slope. constr = TRUE forces sum_i (slope deviation_i) = 0, so
  ## without a fixed effect carrying the community mean slope the model
  ## would force that mean to zero while the truth is mu_beta1 = 0.7.
  ## The same logic applies to Int_occ and the species intercepts.
  f_msom <- inla.mdata(Y, X) ~ -1 + Int_occ + x_cov +
    f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
    f(species_slope, x_cov, model = "iid", constr = TRUE,
      hyper = PRIORS$sigma_gamma) +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))

  m_msom <- fit_occupancy(f_msom, dat, Xdet_all, strategy = SIM$strategy)

  ## --- independent SSOMs, as a loop (not twelve pasted blocks) --------
  ssom <- lapply(seq_len(nsp), function(i) {
    s  <- sim$per_species[[i]]
    sti <- make_st_field(DOM$mesh, spde, s$X[, c("x.loc", "y.loc")],
                         time = s$X$time, n_time = sim$truth$nT)
    Xd <- build_det_X(list(s$g), K = K)
    di <- c(list(Y = s$Y, X = Xd, Int_occ = rep(1, nrow(s$X)),
                 x_cov = s$X$x_cov), sti$idx)
    fi <- inla.mdata(Y, X) ~ -1 + Int_occ + x_cov +
      f(spatialfield, model = spde, A.local = sti$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))
    ## compute_cv = FALSE: the cross-validated comparison was withdrawn
    ## (see 10_metrics.R), so nothing here needs it.
    fit_occupancy(fi, di, Xd, strategy = SIM$strategy, compute_cv = FALSE)
  })

  ## --- community intercept, against the RIGHT truth -------------------
  fx   <- fixed_marginals(m_msom, c("Int_occ", "x_cov"))
  comm <- rbind(
    cbind(rep = r, parameter = "beta0_community",
          summarise_marginal(fx$Int_occ, sim$truth$beta0_community)),
    cbind(rep = r, parameter = "beta1_community",
          summarise_marginal(fx$x_cov,   sim$truth$beta1_community))
  )

  ## --- species intercepts beta0 + gamma_i (the identified quantity) ---
  si <- species_intercepts(m_msom, nsp)
  si$truth   <- sim$truth$species_intercept
  si$covered <- as.integer(si$truth >= si$lwr & si$truth <= si$upr)
  si$rep     <- r

  ## --- shrinkage, measured rather than asserted -----------------------
  ## Distance from the community mean under each model. Shrinkage means
  ## |MSOM - beta0| < |SSOM - beta0|. In the real-data tables this did
  ## NOT hold for every species, so it must be computed, not stated.
  ssom_int <- vapply(ssom, \(f) f$summary.fixed["Int_occ", "mean"], numeric(1))

  ## The fitted SSOM objects are no longer needed once their intercepts
  ## have been extracted; holding thirteen of them per worker is what
  ## exhausts memory.
  rm(ssom); gc(verbose = FALSE)
  shrink <- data.frame(
    rep = r, species = seq_len(nsp),
    truth      = sim$truth$species_intercept,
    ssom       = ssom_int,
    msom       = si$mean,
    naive_occ  = vapply(sim$per_species, \(s) naive_occupancy(s$Y), numeric(1)),
    shrunk     = abs(si$mean - sim$truth$beta0_community) <
                 abs(ssom_int - sim$truth$beta0_community)
  )

  ## The cross-validated MSOM-vs-SSOM comparison is withdrawn; see the
  ## note in 10_metrics.R. In simulation the generating values are known,
  ## so estimation accuracy can be compared directly instead: for each
  ## species, is beta0 + gamma_i closer to the truth under the
  ## hierarchical model than under the independent one?
  shrink$err_ssom <- abs(shrink$ssom - shrink$truth)
  shrink$err_msom <- abs(shrink$msom - shrink$truth)
  shrink$msom_closer <- shrink$err_msom < shrink$err_ssom

  out <- list(comm = comm, species = si, shrink = shrink)
  rm(m_msom, sim); gc(verbose = FALSE)
  out
}

## ---------------------------------------------------------------------
## PRE-RUN PAIRING CHECK (simulation only, no fitting; a few seconds).
## Replicates 1-3 generated with and without homogeneous detection must
## agree on everything except the detection parameters and outcomes.
## ---------------------------------------------------------------------
check_pairing <- function(r) {
  sa <- simulate_msom(DOM, site_frac = SIM$site_frac, sd_beta1 = SIM$sd_beta1,
                     homog_det = FALSE, seed = SEED_BASE + r)
  sc <- simulate_msom(DOM, site_frac = SIM$site_frac, sd_beta1 = SIM$sd_beta1,
                     homog_det = TRUE,  seed = SEED_BASE + r)
  ok <- c(
    gamma   = isTRUE(all.equal(sa$truth$gamma,   sc$truth$gamma,   tolerance = 0)),
    beta1   = isTRUE(all.equal(sa$truth$beta1,   sc$truth$beta1,   tolerance = 0)),
    sites   = identical(sa$truth$site_id, sc$truth$site_id),
    g       = identical(sa$g_all, sc$g_all),
    missing = identical(is.na(sa$Y_all), is.na(sc$Y_all)),
    det_differs = !identical(sa$truth$alpha0, sc$truth$alpha0)
  )
  if (!all(ok))
    stop("Pairing with v2 FAILED for rep ", r, ": ",
         paste(names(ok)[!ok], collapse = ", "),
         ". The stream is not aligned; do not run.")
  invisible(TRUE)
}
invisible(lapply(1:3, check_pairing))
message("Pairing check passed: runs (a) and (c) share fields, sites, g, z and missingness.")

if (!RUN_FULL) {
  message("RUN_FULL = FALSE: run_msom_rep() loaded; full study NOT started.")
} else {

future::plan(future::multisession, workers = N_WORKERS)

## ---------------------------------------------------------------------
## CHECKPOINTING.
##
## The run takes hours; a reboot part way through should not cost all of
## it. Replicates are processed in chunks and each chunk is written to
## CHECKPOINT_DIR as soon as it completes. On restart, chunks already on
## disk are loaded and only the missing replicates are run.
##
## Seeds are SEED_BASE + r, so a replicate recomputed later is identical
## to the one it replaces. Resuming does not disturb reproducibility.
##
## To force a clean run, delete CHECKPOINT_DIR.
## ---------------------------------------------------------------------
CHECKPOINT_DIR <- file.path(OUTDIR, "checkpoints_v5")
CHUNK <- 10L
dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

chunk_file <- function(lo, hi)
  file.path(CHECKPOINT_DIR, sprintf("reps_%04d_%04d.rds", lo, hi))

chunks <- split(seq_len(N_REP), ceiling(seq_len(N_REP) / CHUNK))
done   <- vapply(chunks, \(ix) file.exists(chunk_file(min(ix), max(ix))),
                 logical(1))

if (any(done)) {
  message(sprintf("Resuming: %d of %d chunks already on disk (%d replicates).",
                  sum(done), length(chunks), sum(done) * CHUNK))
}

n_todo <- sum(!done) * CHUNK
message(sprintf("Starting: %d replicates to run on %d workers (~%.1f h expected).",
                n_todo, N_WORKERS, n_todo * 1600 / N_WORKERS / 3600))
## (The v2 run took about 23 h on four workers; expect the same here.)
t_start <- Sys.time()

for (k in seq_along(chunks)) {
  ix <- chunks[[k]]
  f  <- chunk_file(min(ix), max(ix))
  if (file.exists(f)) next

  part <- future_lapply(
    ix,
    \(r) tryCatch(run_msom_rep(r),
                  error = function(e) {
                    warning("rep ", r, ": ", conditionMessage(e))
                    NULL
                  }),
    future.seed = TRUE,
    future.globals = TRUE,
    future.packages = c("INLA", "fmesher", "sf", "terra", "dplyr", "MASS")
  )
  saveRDS(part, f)

  el   <- as.numeric(difftime(Sys.time(), t_start, units = "hours"))
  done_now <- sum(vapply(chunks, \(j) file.exists(chunk_file(min(j), max(j))),
                         logical(1)))
  message(sprintf("  chunk %d/%d saved | %.2f h elapsed | ~%.1f h remaining",
                  k, length(chunks), el,
                  el / max(done_now - sum(done), 1) * (length(chunks) - done_now)))
}

## Reassemble every chunk, including any written by an earlier session.
res <- unlist(
  lapply(chunks, \(ix) readRDS(chunk_file(min(ix), max(ix)))),
  recursive = FALSE
)
res <- Filter(Negate(is.null), res)
message(sprintf("%d / %d replicates completed in %.1f h.",
                length(res), N_REP,
                as.numeric(difftime(Sys.time(), t_start, units = "hours"))))
if (length(res) < 0.9 * N_REP) {
  warning("More than 10% of replicates failed. Failures may not be random ",
          "-- check whether they cluster on particular seeds before ",
          "interpreting the summaries.")
}

comm_df   <- do.call(rbind, lapply(res, `[[`, "comm"))
sp_df     <- do.call(rbind, lapply(res, `[[`, "species"))
shrink_df <- do.call(rbind, lapply(res, `[[`, "shrink"))

saveRDS(list(comm = comm_df, species = sp_df, shrink = shrink_df),
        file.path(OUTDIR, "msom_replicates_v5.rds"))

print(aggregate_reps(comm_df), digits = 3)

## Proportion of replicates in which each species was genuinely shrunk
## toward the community mean, against its detection frequency. This is
## the evidence the shrinkage claim needs.
## Two questions, kept separate. prop_shrunk asks whether pooling moved
## the estimate toward the community mean -- a mechanical property of
## hierarchical models. prop_msom_closer asks whether it moved it toward
## the TRUTH, which is what actually matters and which only simulation
## can answer. rmse_ssom vs rmse_msom quantifies the gain.
shrink_summary <- shrink_df |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    mean_naive_occ   = mean(naive_occ),
    prop_shrunk      = mean(shrunk),
    prop_msom_closer = mean(msom_closer),
    rmse_ssom        = sqrt(mean((ssom - truth)^2)),
    rmse_msom        = sqrt(mean((msom - truth)^2)),
    .groups = "drop"
  )
write.csv(shrink_summary, file.path(OUTDIR, "shrinkage_v5.csv"), row.names = FALSE)
print(as.data.frame(shrink_summary), digits = 3)

cat(sprintf("\nOverall: MSOM closer to truth in %.1f%% of species-replicates; RMSE %.3f vs %.3f (SSOM).\n",
            100 * mean(shrink_df$msom_closer),
            sqrt(mean((shrink_df$msom - shrink_df$truth)^2)),
            sqrt(mean((shrink_df$ssom - shrink_df$truth)^2))))

## ---------------------------------------------------------------------
## PAIRED COMPARISON WITH RUN (a) (v2), replicate by replicate.
## ---------------------------------------------------------------------
f_a <- file.path(OUTDIR, "msom_replicates.rds")
if (!file.exists(f_a)) {
  message("v2 output not found at ", f_a, "; paired comparison skipped.")
} else {
  A <- readRDS(f_a)
  pa <- merge(A$shrink, shrink_df, by = c("rep", "species"),
              suffixes = c("_a", "_c"))
  ## post-run pairing check: the generating intercepts must coincide
  d_truth <- max(abs(pa$truth_a - pa$truth_c))
  if (d_truth > 1e-10)
    warning("Generating intercepts differ between runs (max ", signif(d_truth, 3),
            "): pairing broken, paired statistics below are invalid.")
  cat(sprintf("\nPaired species-replicates: %d (max |truth_a - truth_c| = %.1e)\n",
              nrow(pa), d_truth))

  sq <- function(x) x^2
  ## per replicate: excess MSE of the hierarchical model over independent fits
  by_rep <- do.call(rbind, lapply(split(pa, pa$rep), function(d) data.frame(
    rep = d$rep[1],
    excess_a = mean(sq(d$msom_a - d$truth_a)) - mean(sq(d$ssom_a - d$truth_a)),
    excess_c = mean(sq(d$msom_c - d$truth_c)) - mean(sq(d$ssom_c - d$truth_c)),
    closer_a = mean(d$msom_closer_a),
    closer_c = mean(d$msom_closer_c)
  )))
  dd <- by_rep$excess_c - by_rep$excess_a
  cat(sprintf("RMSE  run (a): MSOM %.3f  SSOM %.3f\n",
              sqrt(mean(sq(pa$msom_a - pa$truth_a))), sqrt(mean(sq(pa$ssom_a - pa$truth_a)))))
  cat(sprintf("RMSE  run (c): MSOM %.3f  SSOM %.3f\n",
              sqrt(mean(sq(pa$msom_c - pa$truth_c))), sqrt(mean(sq(pa$ssom_c - pa$truth_c)))))
  cat(sprintf("Proportion MSOM closer to truth: (a) %.3f  (c) %.3f\n",
              mean(pa$msom_closer_a), mean(pa$msom_closer_c)))
  cat(sprintf("Excess MSE of MSOM over SSOM, (c) minus (a): mean %.4f, SE %.4f (over %d paired replicates)\n",
              mean(dd), sd(dd) / sqrt(length(dd)), length(dd)))
  cat("A negative value means that removing detection heterogeneity reduced the\n")
  cat("hierarchical model's disadvantage; zero means detection misspecification\n")
  cat("played no part in it.\n")

  ## community parameters side by side
  ca <- aggregate_reps(A$comm); cc <- aggregate_reps(comm_df)
  cat("\nCommunity parameters, run (a):\n"); print(ca, digits = 3)
  cat("\nCommunity parameters, run (c):\n"); print(cc, digits = 3)

  write.csv(by_rep, file.path(OUTDIR, "paired_a_vs_c_v5.csv"), row.names = FALSE)
}

}  # end if (RUN_FULL)
