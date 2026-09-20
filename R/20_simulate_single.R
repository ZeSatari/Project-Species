## =====================================================================
## 20_simulate_single.R -- data generation, single-species framework
## =====================================================================
## THREE SCENARIOS. This is the substantive change to the study design.
##
## The manuscript's Contribution 1 rests on the claim that y_B is
## observation-level information. But by construction
##       y_B | z_B ~ Bernoulli(z_B * p_B),
## so y_B always carries information about z_B. The claim is therefore
## conditional, and needs stating as such:
##
##   ASSUMPTION 1.  z_A independent of z_B, given omega(s,t) and x_s.
##
## The previous generator drove BOTH species with the same omega_shared,
## i.e. it sat squarely in the regime where Assumption 1 FAILS. That is
## the regime in which bias in beta is expected -- and indeed the
## reported occupancy MSE rose from 0.086 to 0.147. The simulation
## therefore never tested the hypothesis the paper advances.
##
##   scenario = "independent"  Assumption 1 HOLDS. Independent fields for
##                             z_A and z_B; the DGP contains a genuine
##                             alpha2 * y_B term. Expect: beta unbiased,
##                             alpha2 recovered, coverage near nominal.
##
##   scenario = "shared"       Assumption 1 FAILS. Both species driven by
##                             the same omega. This is the old design.
##                             Expect: bias in beta0, coverage below
##                             nominal. Report this honestly.
##
##   scenario = "confounded"   Assumption 1 genuinely FAILS. A SECOND
##                             latent field, not present in the fitted
##                             model, drives both z_A and z_B. The
##                             fitted omega(s,t) cannot absorb it, so
##                             residual dependence remains and y_B
##                             carries information about z_A.
##                             Expect: bias in beta0, coverage below
##                             nominal. This is the regime the earlier
##                             "shared" scenario was intended to
##                             represent but did not: sharing omega
##                             creates MARGINAL dependence while leaving
##                             z_A and z_B conditionally independent
##                             given the field, which is exactly what
##                             Assumption 1 requires. The R = 200 run
##                             confirmed this: beta0 coverage under
##                             "shared" (0.870) was indistinguishable
##                             from "independent" (0.905).
##
##   scenario = "observer"     The paper's actual ecological motivation:
##                             no alpha2 in the truth at all. Instead a
##                             shared visit-level nuisance u(s,j,t)
##                             (observer / weather / survey context)
##                             perturbs p_A and p_B simultaneously, and
##                             y_B acts as a proxy for it. z_A ind. z_B.
##                             There is no true alpha2 here; the question
##                             is whether including y_B reduces bias in
##                             alpha0, alpha1 and beta. This is the
##                             scenario that justifies the paper's story.
## =====================================================================

## ---------------------------------------------------------------------
## Build the spatial scaffolding once and reuse it across replicates.
## Only the random draws change between replicates.
## ---------------------------------------------------------------------
make_domain <- function(extent = 300, cell_size = 3,
                        max_edge = c(15, 30), offset = c(-0.1, -0.2)) {
  grd <- sf::st_make_grid(
    sf::st_as_sfc(sf::st_bbox(c(xmin = 0, ymin = 0,
                                xmax = extent, ymax = extent))),
    cellsize = c(cell_size, cell_size)
  ) |>
    sf::st_sf() |>
    dplyr::mutate(cellid = dplyr::row_number())

  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grd)))

  mesh <- fmesher::fm_mesh_2d(
    loc.domain = matrix(c(0, 0, extent, 0, extent, extent, 0, extent),
                        ncol = 2, byrow = TRUE),
    offset = offset, max.edge = max_edge
  )

  message(sprintf("Domain: %d cells, max.edge = c(%g, %g), mesh with %d nodes.",
                  nrow(grd), max_edge[1], max_edge[2], mesh$n))
  list(grid = grd, coords = coords, mesh = mesh,
       ncells = nrow(grd), extent = extent)
}

## ---------------------------------------------------------------------
## One stationary AR(1)-in-time Matern field on the grid.
## Marginal variance is sigma^2 at every t (innovation scaled by
## sqrt(1 - rho^2)). State this parameterisation in Methods: it is what
## control.group = list(model = "ar1") assumes, and it determines whether
## sigma is the innovation SD or the marginal SD.
## ---------------------------------------------------------------------
sim_st_field <- function(dom, range, sigma, rho, nT, seed) {
  spde <- INLA::inla.spde2.pcmatern(dom$mesh,
                                    prior.range = c(range, 0.5),
                                    prior.sigma = c(sigma, 0.5))
  Q <- INLA::inla.spde.precision(spde, theta = c(log(range), log(sigma)))
  A <- INLA::inla.spde.make.A(dom$mesh, loc = dom$coords)

  eps <- as.matrix(A %*% INLA::inla.qsample(n = nT, Q = Q, seed = seed))
  om  <- eps
  ## NOTE: `for (t in 2:nT)` runs BACKWARDS as 2:1 when nT == 1. This
  ## function is called with nT = 1 to draw the static covariates, so the
  ## sequence must be guarded.
  if (nT > 1L) {
    for (t in seq.int(2L, nT)) {
      om[, t] <- rho * om[, t - 1] + sqrt(1 - rho^2) * eps[, t]
    }
  }
  om
}

## ---------------------------------------------------------------------
## Main generator.
##
## Returns the data AND a `truth` list. Every quantity the manuscript
## compares an estimate against must come from `truth` -- never from a
## number typed into the text. (The MSOM generator previously drew gamma
## with mean -0.5 while adding mu_beta0 = -0.5 separately, so the true
## community intercept was -1.0, but the manuscript compared the estimate
## against -0.5.)
## ---------------------------------------------------------------------
simulate_single <- function(dom,
                            scenario = c("independent", "shared", "observer", "confounded"),
                            nT        = 5,
                            K         = 3,
                            site_frac = 0.20,
                            beta      = c(beta0 = -0.5, beta1 = 0.3),
                            alpha_A   = c(a0 = qlogis(0.30), a1 = 0.4),
                            alpha_B   = c(a0 = qlogis(0.25), a1 = 0.5),
                            alpha2    = 0.7,
                            eta_obs   = 1.0,   # used only by "observer"
                            sd_obs    = 0.8,   # SD of the shared nuisance
                            eta_conf  = 1.0,   # used only by "confounded"
                            range_conf = 25,   # short range: NOT absorbed
                            sigma_conf = 0.8,  #   by the fitted field
                            range     = 100,
                            sigma     = 0.5,
                            rho       = 0.65,
                            missing_p = 0.10,
                            seed) {

  scenario <- match.arg(scenario)
  set.seed(seed)
  nc <- dom$ncells

  ## --- latent fields --------------------------------------------------
  omega_A <- sim_st_field(dom, range, sigma, rho, nT, seed = seed)
  omega_B <- if (scenario == "shared") {
    omega_A            # marginal dependence; A1 still holds (see header)
  } else {
    sim_st_field(dom, range, sigma, rho, nT, seed = seed + 500000L)
  }

  ## Unmodelled confounder. Deliberately short-range relative to the
  ## fitted field, so the fitted omega(s,t) cannot absorb it and
  ## conditional independence of z_A and z_B genuinely fails.
  conf <- if (scenario == "confounded") {
    sim_st_field(dom, range_conf, sigma_conf, rho, nT,
                 seed = seed + 700000L)
  } else {
    matrix(0, dom$ncells, nT)
  }

  ## --- environmental covariates --------------------------------------
  ## Drawn from independent fields, and NOT reused as the random effect.
  ## (Previously x_s and omega_s were the same object, and x_s2 was drawn
  ##  with the same seed as the shared field, so the covariate was very
  ##  likely a time slice of the field itself.)
  x_A <- as.numeric(scale(sim_st_field(dom, range, sigma, 0, 1,
                                       seed = seed + 100000L)[, 1]))
  x_B <- as.numeric(scale(sim_st_field(dom, range, sigma, 0, 1,
                                       seed = seed + 200000L)[, 1]))
  ## Guard against accidentally re-using the field as the covariate. This
  ## is a WARNING, not an error: with range = 100 on a 300 x 300 domain
  ## there are few effective independent regions, so two genuinely
  ## independent realisations can show a large sample correlation by
  ## chance. Killing the replicate would bias the retained set.
  r_xf <- cor(x_A, omega_A[, 1])
  if (abs(r_xf) > 0.8) {
    warning(sprintf("seed %d: cor(x_A, omega_A[,1]) = %.2f -- check the seeds are distinct.",
                    seed, r_xf))
  }

  ## --- ecological state process --------------------------------------
  ## Quadratic in x, matching the fitted model. The fitted model MUST use
  ## I(x^2): in an R formula, `x^2` is the crossing operator and collapses
  ## to `x`, so the previous pipeline generated a quadratic and fitted a
  ## linear term.
  lin_A <- beta["beta0"] + beta["beta1"] * x_A^2
  lin_B <- -0.5 + 0.4 * x_B^2
  zA <- zB <- matrix(NA_integer_, nc, nT)
  for (t in seq_len(nT)) {
    zA[, t] <- rbinom(nc, 1, plogis(lin_A + omega_A[, t] + eta_conf * conf[, t]))
    zB[, t] <- rbinom(nc, 1, plogis(lin_B + omega_B[, t] + eta_conf * conf[, t]))
  }

  ## --- site selection -------------------------------------------------
  nsites  <- round(nc * site_frac)
  site_id <- sort(sample.int(nc, nsites))

  ## --- observation process -------------------------------------------
  gA <- array(runif(nsites * K * nT, -1, 1), c(nsites, K, nT))
  gB <- array(runif(nsites * K * nT, -1, 1), c(nsites, K, nT))
  u  <- if (scenario == "observer")
          array(rnorm(nsites * K * nT, 0, sd_obs), c(nsites, K, nT))
        else array(0, c(nsites, K, nT))

  YA <- YB <- matrix(NA_integer_, nsites * nT, K)
  for (t in seq_len(nT)) {
    rows <- ((t - 1) * nsites + 1):(t * nsites)
    for (j in seq_len(K)) {
      pB <- plogis(alpha_B["a0"] + alpha_B["a1"] * gB[, j, t] +
                     eta_obs * u[, j, t] * (scenario == "observer"))
      yB <- rbinom(nsites, 1, zB[site_id, t] * pB)

      ## Missingness in y_B. NOTE: the fitted model applies the
      ## "missing -> non-detection" convention. Here the TRUTH uses the
      ## complete y_B, so that convention is genuinely tested rather than
      ## assumed correct (previously the truth itself used the NA->0
      ## version, making the convention self-fulfilling).
      etaA <- alpha_A["a0"] + alpha_A["a1"] * gA[, j, t]
      etaA <- etaA + switch(scenario,
                            observer = eta_obs * u[, j, t],
                            alpha2 * yB)   # default: independent/shared/confounded
      yA <- rbinom(nsites, 1, zA[site_id, t] * plogis(etaA))

      yB[as.logical(rbinom(nsites, 1, missing_p))] <- NA
      yA[as.logical(rbinom(nsites, 1, missing_p))] <- NA
      YB[rows, j] <- yB
      YA[rows, j] <- yA
    }
  }

  ## --- assemble, with an explicit alignment contract ------------------
  dat <- data.frame(
    cellid = rep(site_id, times = nT),
    time   = rep(seq_len(nT), each = nsites),
    x.loc  = rep(dom$coords[site_id, 1], times = nT),
    y.loc  = rep(dom$coords[site_id, 2], times = nT),
    x_A    = rep(x_A[site_id], times = nT),
    x_B    = rep(x_B[site_id], times = nT)
  )
  stopifnot(nrow(dat) == nrow(YA), nrow(dat) == nrow(YB))
  ## covariate travels WITH the row; it is never re-attached by position

  colnames(YA) <- paste0("yA", seq_len(K))
  colnames(YB) <- paste0("yB", seq_len(K))
  gA_long <- do.call(rbind, lapply(seq_len(nT), \(t) gA[, , t]))
  gB_long <- do.call(rbind, lapply(seq_len(nT), \(t) gB[, , t]))
  colnames(gA_long) <- paste0("gA", seq_len(K))
  colnames(gB_long) <- paste0("gB", seq_len(K))

  list(
    data     = cbind(dat, YA, YB, gA_long, gB_long),
    Y_A      = YA, Y_B = YB, gA = gA_long,
    scenario = scenario,
    truth    = list(
      beta0  = unname(beta["beta0"]),
      beta1  = unname(beta["beta1"]),
      alpha0 = unname(alpha_A["a0"]),
      alpha1 = unname(alpha_A["a1"]),
      ## no true alpha2 exists under "observer": y_B is a proxy for u,
      ## not a term in the data-generating detection model
      alpha2 = if (scenario == "observer") NA_real_ else alpha2,
      range = range, sigma = sigma, rho = rho,
      nsites = nsites, nT = nT, K = K, seed = seed
    )
  )
}
