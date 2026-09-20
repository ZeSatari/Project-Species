## =====================================================================
## 50_run_msom_sim_v4.R -- heterogeneous prevalence
## =====================================================================
## ATTEMPTED RUN. Reported in Section 4.5 as an attempt that did not
## resolve the question, not as a result. Does not supersede v2 or v3.
## Confirm on load that the banner reports "v4 | heterogeneous
## prevalence".
##
## WHY THIS RUN EXISTS. Versions 2 and 3 simulated communities whose
## species had comparable prevalence, and found that pooling did NOT
## improve species-level accuracy (RMSE 0.247 vs 0.219 independent). The
## empirical community is nothing like that: prevalence runs from 2.3%
## (NAWA) to 73.3% (BTNW), and there pooling appeared to help. The
## manuscript therefore rested its central multi-species claim on an
## empirical observation with no simulation behind it. This run supplies
## the missing test.
##
## DESIGN, calibrated to hbefTrends rather than invented:
##   logit(naive occupancy) across the twelve species: mean -1.14, sd 1.77
##   logit(detection rate)                           : mean -2.10, sd 1.68
##   correlation between the two                     : 0.997
##
## That correlation is the point. Rare species at Hubbard Brook are also
## hard to detect; a simulation varying only occupancy would not
## reproduce the regime. Species effects are therefore drawn from a
## bivariate normal with correlation 0.9, and the detection intercept
## carries the same species-level deviation, scaled.
##
## EXPECT FAILURES. With species this sparse, independent fits may not
## converge -- NAWA did not, in the empirical data, under the initial
## mesh. Non-convergence is recorded rather than allowed to halt the run,
## and is itself part of the result.
##
## WHAT VARIES AND WHAT DOES NOT. Species-specific slopes are held at
## sd_beta1 = 0, as in v3, so that this run isolates heterogeneity in
## prevalence. The comparison of interest is v4 against v3: same slope
## structure, different prevalence spread.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")   # make_domain(), sim_st_field()
source("R/30_fit_occupancy.R")

OUTDIR <- "results/msom"
record_session(OUTDIR, tag = "heteroprev")

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
## sd_beta1 = 0, carried over from v3: slopes stay homogeneous so that
## this run isolates heterogeneity in PREVALENCE. Varying both at once
## would leave any difference from v2 unattributable between the two.
SIM <- list(site_frac = 0.02, max_edge = c(25, 50),
            sd_beta1 = 0,
            strategy = "simplified.laplace")

## Each replicate fits 13 models (one MSOM + twelve SSOMs). With
## config = TRUE every fitted object retains the full latent-field
## structure, so a worker can hold several GB at once. Six workers
## exhausted memory on a 32 GB machine ("allocation of text connection
## failed"). Three is a safer starting point; raise it only after
## watching peak usage.
N_WORKERS <- 4   # peak memory measured at ~2.6 GB per replicate

## ---------------------------------------------------------------------
## SPECIES-EFFECT PRIOR, WIDENED FOR THIS RUN ONLY.
##
## PRIORS$sigma_gamma states P(sigma_gamma > 1) = 0.01. That suits v2 and
## v3, where the generating value is 0.8. Here the community is built to
## match Hubbard Brook and the realised sd(gamma) is about 1.4, which the
## original prior places far in its tail: the model would be told, with
## 99% prior confidence, that the truth is impossible. A single test
## replicate under the original prior produced species intercepts that
## overshot the truth rather than shrinking toward the community mean
## (e.g. -5.55 against a generating -3.67), which is what prior-data
## conflict of this kind looks like.
##
## The prior is therefore widened to P(sigma_gamma > 3) = 0.01, which
## accommodates 1.4 without being uninformative.
##
## CONSEQUENCE FOR THE COMPARISON: v4 differs from v3 in prior as well as
## in prevalence. The manuscript must say so. The alternative -- rerunning
## v2 and v3 under the wider prior, another 44 hours -- was not taken.
## ---------------------------------------------------------------------
PRIORS_V4_GAMMA <- list(prec = list(prior = "pc.prec", param = c(3, 0.01)))

message("=== 50_run_msom_sim_v4 | heterogeneous prevalence | widened gamma prior | N_REP = ", N_REP,
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
                          mu_beta0 = -1.14, sd_gamma = 1.77, gen_cor = 0,
                          mu_beta1 = 0.7,  sd_beta1 = 0.5,
                          mu_alpha0 = -2.10, sd_alpha0 = 1.68,
                          cor_occ_det = 0.9,
                          range = 100, sigma = 0.5, rho = 0.6,
                          missing_p = 0.10, seed) {
  set.seed(seed)
  nc <- dom$ncells

  omega <- sim_st_field(dom, range, sigma, rho, nT, seed = seed)
  x_s   <- as.numeric(scale(sim_st_field(dom, range, sigma, 0, 1,
                                         seed = seed + 100000L)[, 1]))

  ## Species effects centred at zero, then explicitly recentred so the
  ## realised community intercept is unambiguous for THIS replicate.
  ## Occupancy and detection deviations drawn jointly. At Hubbard Brook
  ## the two correlate 0.997 on the logit scale: species that are rare
  ## are also hard to detect. Drawing them independently would produce a
  ## community that does not resemble the data this run is meant to
  ## represent.
  S <- matrix(c(sd_gamma^2,
                cor_occ_det * sd_gamma * sd_alpha0,
                cor_occ_det * sd_gamma * sd_alpha0,
                sd_alpha0^2), 2, 2)
  eff <- MASS::mvrnorm(n_species, mu = c(0, 0), Sigma = S)
  gamma      <- eff[, 1] - mean(eff[, 1])   # matches constr = TRUE in the fit
  alpha0_dev <- eff[, 2] - mean(eff[, 2])
  ## Centre the slopes for the same reason: with constr = TRUE on the
  ## random slope, the identified community slope is mean(beta1). Recentre
  ## so that the realised truth for THIS replicate is exactly mu_beta1
  ## rather than mu_beta1 + sampling noise.
  beta1 <- rnorm(n_species, mu_beta1, sd_beta1)
  beta1 <- beta1 - mean(beta1) + mu_beta1

  alpha0 <- mu_alpha0 + alpha0_dev
  alpha1 <- rnorm(n_species, -1.0, 0.3)

  ## Floor on detectability. With sd = 1.68 the lower tail can put a
  ## species' per-visit detection probability near zero, so that it is
  ## never recorded anywhere. Such a species carries no information for
  ## any model and would enter the comparison as a degenerate case rather
  ## than an informative one. The floor corresponds to the least
  ## detectable species at Hubbard Brook (NAWA, 0.009).
  alpha0 <- pmax(alpha0, qlogis(0.009))

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
      alpha0 = alpha0, alpha0_dev = alpha0_dev,
      beta0_community = mu_beta0,       # gamma is centred, so this is it
      beta1_community = mu_beta1,       # beta1 is centred, likewise
      gamma = gamma, beta1 = beta1,
      alpha0 = alpha0, alpha1 = alpha1,
      sd_gamma = sd_gamma, gen_cor = gen_cor,
      species_intercept = mu_beta0 + gamma,
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
  ## sd_beta1 = 0, carried over from v3: slopes stay homogeneous so that
  ## this run isolates heterogeneity in PREVALENCE. Letting both vary at
  ## once would leave any difference from v2 unattributable between the
  ## two sources.
  sim <- simulate_msom(DOM, site_frac = SIM$site_frac,
                       sd_beta1 = SIM$sd_beta1, seed = SEED_BASE + r)
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
  ## The random slope is retained even though the generating slopes are
  ## now homogeneous: the fitted model must be the one the manuscript
  ## describes. With sd_beta1 = 0 the term should shrink to near zero,
  ## which is itself worth checking in the output.
  ## The fixed `x_cov` term is REQUIRED alongside the constrained random
  ## slope. constr = TRUE forces sum_i (slope deviation_i) = 0, so
  ## without a fixed effect carrying the community mean slope the model
  ## would force that mean to zero while the truth is mu_beta1 = 0.7.
  ## The same logic applies to Int_occ and the species intercepts.
  f_msom <- inla.mdata(Y, X) ~ -1 + Int_occ + x_cov +
    f(species_id, model = "iid", constr = TRUE, hyper = PRIORS_V4_GAMMA) +
    f(species_slope, x_cov, model = "iid", constr = TRUE,
      hyper = PRIORS_V4_GAMMA) +
    f(spatialfield, model = spde, A.local = stf$A,
      group = spatialfield.group,
      control.group = list(model = "ar1", hyper = PRIORS$ar1))

  m_msom <- fit_occupancy(f_msom, dat, Xdet_all, strategy = SIM$strategy)

  ## --- independent SSOMs, as a loop (not twelve pasted blocks) --------
  ## With very sparse species an independent fit can fail outright. That is
  ## part of the result, not an obstacle: record it and carry on.
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
    tryCatch(
      fit_occupancy(fi, di, Xd, strategy = SIM$strategy, compute_cv = FALSE),
      error = function(e) NULL
    )
  })
  ssom_ok <- !vapply(ssom, is.null, logical(1))

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
  ssom_int <- vapply(seq_len(nsp), function(i)
    if (ssom_ok[i]) ssom[[i]]$summary.fixed["Int_occ", "mean"] else NA_real_,
    numeric(1))

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
                 abs(ssom_int - sim$truth$beta0_community),
    ssom_converged = ssom_ok
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
CHECKPOINT_DIR <- file.path(OUTDIR, "checkpoints_v4")
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
        file.path(OUTDIR, "msom_replicates_v4.rds"))

print(aggregate_reps(comm_df), digits = 3)

## Proportion of replicates in which each species was genuinely shrunk
## toward the community mean, against its detection frequency. This is
## the evidence the shrinkage claim needs.
## Two questions, kept separate. prop_shrunk asks whether pooling moved
## the estimate toward the community mean -- a mechanical property of
## hierarchical models. prop_msom_closer asks whether it moved it toward
## the TRUTH, which is what actually matters and which only simulation
## can answer. rmse_ssom vs rmse_msom quantifies the gain.
## Comparisons are restricted to replicates in which the independent fit
## converged; the failure rate is reported separately, since a species
## that cannot be fitted alone is the clearest case for pooling.
shrink_summary <- shrink_df |>
  dplyr::group_by(species) |>
  dplyr::summarise(
    mean_naive_occ   = mean(naive_occ),
    prop_ssom_failed = mean(!ssom_converged),
    prop_shrunk      = mean(shrunk[ssom_converged]),
    prop_msom_closer = mean(msom_closer[ssom_converged]),
    rmse_ssom        = sqrt(mean((ssom[ssom_converged] -
                                  truth[ssom_converged])^2)),
    rmse_msom_paired = sqrt(mean((msom[ssom_converged] -
                                  truth[ssom_converged])^2)),
    rmse_msom_all    = sqrt(mean((msom - truth)^2)),
    .groups = "drop"
  )
write.csv(shrink_summary, file.path(OUTDIR, "shrinkage_v4.csv"), row.names = FALSE)
print(as.data.frame(shrink_summary), digits = 3)

ok <- shrink_df$ssom_converged
cat(sprintf("\nIndependent fits failed in %.1f%% of species-replicates.\n",
            100 * mean(!ok)))
cat(sprintf("On the %d pairs where both converged: MSOM closer to truth in %.1f%%; RMSE %.3f vs %.3f (SSOM).\n",
            sum(ok), 100 * mean(shrink_df$msom_closer[ok]),
            sqrt(mean((shrink_df$msom[ok] - shrink_df$truth[ok])^2)),
            sqrt(mean((shrink_df$ssom[ok] - shrink_df$truth[ok])^2))))
cat(sprintf("MSOM over all %d species-replicates, including those the independent fit could not reach: RMSE %.3f\n",
            nrow(shrink_df),
            sqrt(mean((shrink_df$msom - shrink_df$truth)^2))))
cat("\nCompare with v2 (comparable prevalence): 45.3%% closer, RMSE 0.247 vs 0.219.\n")
cat("If pooling helps here and not there, the claim that its value depends on\n")
cat("the data configuration has simulation support rather than only empirical.\n")

}  # end if (RUN_FULL)
