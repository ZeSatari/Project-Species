## =====================================================================
## 99_diag_confounding.R -- why is beta1_community biased by -46%?
## =====================================================================
## A single MSOM replicate returned beta1 = 0.380 against a generating
## value of 0.700, with a 95% interval of (0.310, 0.450) that excludes
## the truth, and pooling moved species intercepts AWAY from their true
## values in 9 of 12 species. All three are consistent with one cause:
## the environmental covariate and the latent field are not separable.
##
## In simulate_msom() both are Matern fields with range = 100 on a
## 300 x 300 domain -- roughly nine independent regions each. Two smooth
## surfaces of the same characteristic scale are nearly interchangeable, so
## the field absorbs the covariate effect and beta1 is attenuated.
##
## Three fits per configuration, one replicate. Runs in minutes because
## the SSOM loop is skipped entirely.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")
source("R/30_fit_occupancy.R")

SIM <- list(site_frac = 0.02, max_edge = c(25, 50),
            strategy = "simplified.laplace")
DOM <- make_domain(max_edge = SIM$max_edge)
SEED <- 20260201L + 1L

## ---------------------------------------------------------------------
## Generator with the covariate range exposed, so it can be varied
## independently of the field. Everything else matches simulate_msom().
## ---------------------------------------------------------------------
sim_msom_diag <- function(dom, range_x, n_species = 12, nT = 5, K = 3,
                          site_frac = 0.02,
                          mu_beta0 = -0.5, sd_gamma = 0.8,
                          mu_beta1 = 0.7,  sd_beta1 = 0.5,
                          range = 100, sigma = 0.5, rho = 0.6,
                          missing_p = 0.10, seed) {
  set.seed(seed)
  nc <- dom$ncells

  omega <- sim_st_field(dom, range, sigma, rho, nT, seed = seed)
  x_s   <- as.numeric(scale(sim_st_field(dom, range_x, sigma, 0, 1,
                                         seed = seed + 100000L)[, 1]))

  gamma <- as.numeric(MASS::mvrnorm(1, rep(0, n_species),
                                    diag(sd_gamma^2, n_species)))
  gamma <- gamma - mean(gamma)
  beta1 <- rnorm(n_species, mu_beta1, sd_beta1)
  beta1 <- beta1 - mean(beta1) + mu_beta1

  alpha0 <- rnorm(n_species, -0.9, 0.2)
  alpha1 <- rnorm(n_species, -1.0, 0.3)

  nsites  <- round(nc * site_frac)
  site_id <- sort(sample.int(nc, nsites))

  out <- vector("list", n_species)
  for (i in seq_len(n_species)) {
    Yi <- matrix(NA_integer_, nsites * nT, K)
    gi <- array(runif(nsites * K * nT, -1, 1), c(nsites, K, nT))
    for (t in seq_len(nT)) {
      psi <- plogis(mu_beta0 + gamma[i] + beta1[i] * x_s^2 + omega[, t])
      z   <- rbinom(nc, 1, psi)
      rows <- ((t - 1) * nsites + 1):(t * nsites)
      for (j in seq_len(K)) {
        y <- rbinom(nsites, 1, z[site_id] * plogis(alpha0[i] + alpha1[i] * gi[, j, t]))
        y[as.logical(rbinom(nsites, 1, missing_p))] <- NA
        Yi[rows, j] <- y
      }
    }
    out[[i]] <- list(
      Y = Yi,
      g = do.call(rbind, lapply(seq_len(nT), \(t) gi[, , t])),
      X = data.frame(species_id = i,
                     time  = rep(seq_len(nT), each = nsites),
                     x.loc = rep(dom$coords[site_id, 1], times = nT),
                     y.loc = rep(dom$coords[site_id, 2], times = nT),
                     x_sq  = rep(x_s[site_id]^2, times = nT))
    )
  }

  list(Y_all = do.call(rbind, lapply(out, `[[`, "Y")),
       g_all = do.call(rbind, lapply(out, `[[`, "g")),
       X_all = do.call(rbind, lapply(out, `[[`, "X")),
       x_full = x_s, omega1 = omega[, 1],
       truth = list(beta0 = mu_beta0, beta1 = mu_beta1, nT = nT, K = K,
                    n_species = n_species))
}

## ---------------------------------------------------------------------
## Fit, optionally without the spatial field. If beta1 is recovered when
## the field is dropped but attenuated when it is present, the covariate
## and the field are competing for the same variation.
## ---------------------------------------------------------------------
fit_diag <- function(sim, with_field = TRUE) {
  X <- sim$X_all; K <- sim$truth$K
  spde <- INLA::inla.spde2.pcmatern(DOM$mesh,
                                    prior.range = PRIORS$spde_sim$range,
                                    prior.sigma = PRIORS$spde_sim$sigma)
  stf <- make_st_field(DOM$mesh, spde, X[, c("x.loc", "y.loc")],
                       time = X$time, n_time = sim$truth$nT)
  Xd  <- build_det_X(list(sim$g_all), K = K)

  dat <- c(list(Y = sim$Y_all, X = Xd,
                Int_occ = rep(1, nrow(X)), x_sq = X$x_sq,
                species_id = X$species_id, species_slope = X$species_id),
           stf$idx)

  f <- if (with_field) {
    inla.mdata(Y, X) ~ -1 + Int_occ + x_sq +
      f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(species_slope, x_sq, model = "iid", constr = TRUE,
        hyper = PRIORS$sigma_gamma) +
      f(spatialfield, model = spde, A.local = stf$A,
        group = spatialfield.group,
        control.group = list(model = "ar1", hyper = PRIORS$ar1))
  } else {
    inla.mdata(Y, X) ~ -1 + Int_occ + x_sq +
      f(species_id, model = "iid", constr = TRUE, hyper = PRIORS$sigma_gamma) +
      f(species_slope, x_sq, model = "iid", constr = TRUE,
        hyper = PRIORS$sigma_gamma)
  }

  fit_occupancy(f, dat, Xd, strategy = SIM$strategy, compute_cv = FALSE)
}

report <- function(label, fit, truth) {
  s <- fit$summary.fixed[c("Int_occ", "x_sq"), c("mean", "0.025quant", "0.975quant")]
  data.frame(
    config = label,
    parameter = c("beta0", "beta1"),
    truth = c(truth$beta0, truth$beta1),
    est = s$mean, lwr = s$`0.025quant`, upr = s$`0.975quant`,
    covered = as.integer(c(truth$beta0, truth$beta1) >= s$`0.025quant` &
                         c(truth$beta0, truth$beta1) <= s$`0.975quant`),
    row.names = NULL
  )
}

## ---------------------------------------------------------------------
## Three configurations.
## ---------------------------------------------------------------------
res <- list()

## (a) as-is: covariate range 100, same as the field
s_a <- sim_msom_diag(DOM, range_x = 100, site_frac = SIM$site_frac, seed = SEED)
cat(sprintf("\n[a] range_x = 100 | cor(x, omega_1) = %.3f\n",
            cor(s_a$x_full, s_a$omega1)))
res$a1 <- report("a: range_x=100, with field",    fit_diag(s_a, TRUE),  s_a$truth)
res$a2 <- report("a: range_x=100, NO field",      fit_diag(s_a, FALSE), s_a$truth)

## (b) short-range covariate: structure the smooth field cannot mimic
s_b <- sim_msom_diag(DOM, range_x = 20, site_frac = SIM$site_frac, seed = SEED)
cat(sprintf("[b] range_x = 20  | cor(x, omega_1) = %.3f\n",
            cor(s_b$x_full, s_b$omega1)))
res$b1 <- report("b: range_x=20, with field",     fit_diag(s_b, TRUE),  s_b$truth)

## (c) independent covariate with no spatial structure at all
s_c <- sim_msom_diag(DOM, range_x = 5, site_frac = SIM$site_frac, seed = SEED)
cat(sprintf("[c] range_x = 5   | cor(x, omega_1) = %.3f\n\n",
            cor(s_c$x_full, s_c$omega1)))
res$c1 <- report("c: range_x=5, with field",      fit_diag(s_c, TRUE),  s_c$truth)

out <- do.call(rbind, res)
print(out, row.names = FALSE, digits = 3)

cat("\nReading this:\n")
cat("  a-with-field biased but a-no-field unbiased -> the covariate and the\n")
cat("    field are competing; the generator, not the estimator, is at fault.\n")
cat("  b and c recover beta1 -> confirmed: a covariate on the same spatial\n")
cat("    scale as the field is not separable from it.\n")
cat("  all three biased -> the cause lies elsewhere; do not change the\n")
cat("    generator until it is found.\n")

dir.create("results/msom", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "results/msom/confounding_diag.csv", row.names = FALSE)
