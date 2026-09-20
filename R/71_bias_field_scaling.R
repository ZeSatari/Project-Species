## =====================================================================
## 71_bias_field_scaling.R -- is the slope bias attenuation by the field?
## =====================================================================
## Appendix S1.1 reports a persistent negative bias in the environmental
## slope (about -0.024) whose source 68_bias_source.R and 69_bias_field.R
## did not find. Those two varied the FITTED model. This one varies the
## GENERATING field, which they did not.
##
## THE HYPOTHESIS. The model integrates over the latent field. Whatever
## part of omega the fitted field fails to recover acts as residual
## randomness in the linear predictor, and marginalising over a random
## effect under a logit link shrinks a slope toward zero. The sign fits:
## beta1 = 0.3 > 0 and the bias is negative. If this is the mechanism the
## bias must vanish when the generating field has no variance, and grow
## with it.
##
## THE PREDICTION. For a logistic model with a random effect of standard
## deviation sigma_r the usual approximation is
##
##     beta_marginal ~= beta_conditional * (1 + 0.346 * sigma_r^2)^(-1/2)
##
## with c = 16*sqrt(3)/(15*pi), so c^2 = 768/(225*pi^2) = 0.3458. Source:
## Zeger, Liang & Albert (1988), Biometrics 44(4):1049-1060.
##
## That result is derived for a random INTERCEPT. The general form carries
## c^2 z' G z, and omega(s,t) here is a spatio-temporal field, not an
## intercept. Treat the formula as an upper bound rather than a point
## prediction, and say so if it reaches the manuscript. Moreover sigma_r
## is not the generating sigma but the part the fitted field fails to
## absorb, so the observed attenuation should fall between zero and the
## line the formula draws. A result above that line, or one that does not
## scale with sigma at all, falsifies the hypothesis -- worth knowing
## just as much.
##
## WHAT IS HELD FIXED. simulate_single() generates the occupancy
## covariate as x_A^2, and the quadratic form has its own attenuation
## through the chi-squared geometry of x^2, already reported in Appendix
## S1.1. That form is held fixed here, so any scaling with sigma is
## attributable to the field rather than to the covariate.
##
## Everything below uses the project's own generator and its own fitting
## entry point. Nothing is reimplemented.
## =====================================================================

source("R/00_setup.R")
source("R/10_metrics.R")
source("R/20_simulate_single.R")
source("R/30_fit_occupancy.R")

suppressPackageStartupMessages(library(future.apply))

OUTDIR <- "results/bias_scaling"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
record_session(OUTDIR, tag = "bias_scaling")

## ---------------------------------------------------------------------
## Configuration. Deliberately smaller than the study's own simulation:
## this is a diagnostic, and the question is how the bias moves with
## sigma, not its value to three decimal places.
## ---------------------------------------------------------------------
SD_FIELD  <- c(0.001, 0.25, 0.50, 0.75, 1.00)   # 0.001 rather than 0:
                                                # sim_st_field takes log(sigma)
R_REP     <- 40
N_WORKERS <- 4
SEED_BASE <- 771000

SIM <- list(nT = 3, K = 3, site_frac = 0.06, range = 100, rho = 0.65)

DOM <- make_domain(extent = 300, cell_size = 3, max_edge = c(25, 50))
message(sprintf("=== 71_bias_field_scaling | sigma in {%s} | R = %d ===",
                paste(SD_FIELD, collapse = ", "), R_REP))

## Prior list name differs between pipelines; take whichever exists.
sp_prior <- if (!is.null(PRIORS$spde_sim)) PRIORS$spde_sim else PRIORS$spde_real
spde <- INLA::inla.spde2.pcmatern(DOM$mesh,
                                  prior.range = sp_prior$range,
                                  prior.sigma = sp_prior$sigma,
                                  constr = TRUE)

## ---------------------------------------------------------------------
## One replicate: generate at the requested sigma, fit the baseline model
## M1, return the slope.
## ---------------------------------------------------------------------
one_rep <- function(r, sd_field) {
  ## future::multisession workers are fresh R sessions: they inherit
  ## exported globals but not package attachments or inla.setOption().
  ## 00_setup.R provides init_worker() for exactly this, and it must be
  ## called at the top of every function passed to future_lapply().
  ## Omitting it here made all 200 parallel replicates fail while the
  ## sequential probe succeeded.
  init_worker()

  sim <- tryCatch(
    simulate_single(DOM, scenario = "independent",
                    nT = SIM$nT, K = SIM$K, site_frac = SIM$site_frac,
                    range = SIM$range, sigma = sd_field, rho = SIM$rho,
                    seed = SEED_BASE + r),
    error = function(e) e)
  if (inherits(sim, "error")) return(NULL)

  d <- sim$data
  ## The fitted covariate must be I(x^2): in a formula, x^2 is the
  ## crossing operator and collapses to x. The generator's own header
  ## records that the previous pipeline made exactly this mistake.
  d$xq      <- d$x_A^2
  d$Int_occ <- 1

  stf  <- make_st_field(DOM$mesh, spde, d[, c("x.loc", "y.loc")],
                        time = d$time, n_time = SIM$nT)
  ar1  <- list(rho = list(prior = "pc.cor0", param = PRIORS$ar1$rho$param))
  Xdet <- build_det_X(list(sim$gA), K = SIM$K)
  dat  <- c(as.list(d), stf$idx, list(Y = sim$Y_A, X = Xdet))

  f <- as.formula(paste(
    "inla.mdata(Y, X) ~ -1 + Int_occ + xq +",
    "f(spatialfield, model = spde, A.local = stf$A,",
    "  group = spatialfield.group,",
    "  control.group = list(model = 'ar1', hyper = ar1))"))
  environment(f) <- environment()

  fit <- tryCatch(fit_occupancy(f, dat, Xdet, compute_cv = FALSE),
                  error = function(e) e)
  if (inherits(fit, "error")) return(NULL)

  m  <- fit$summary.fixed["xq", ]
  hp <- fit$summary.hyperpar
  gr <- function(p) { i <- grep(p, rownames(hp))
                      if (length(i)) hp[i[1], "mean"] else NA_real_ }
  data.frame(rep = r, sd_field = sd_field, truth = sim$truth$beta1,
             est = m[["mean"]], sd = m[["sd"]],
             covered = as.integer(m[["0.025quant"]] <= sim$truth$beta1 &
                                  sim$truth$beta1 <= m[["0.975quant"]]),
             sigma_hat = gr("^Stdev for"), range_hat = gr("^Range for"))
}

## ---------------------------------------------------------------------
## Time one replicate before committing to the grid
## ---------------------------------------------------------------------
message("timing one replicate...")
t0    <- proc.time()[3]
probe <- one_rep(1, SD_FIELD[3])
el    <- proc.time()[3] - t0
if (is.null(probe)) stop("the probe replicate failed; fix that before running")
message(sprintf("one fit: %.0f s | full run ~%.1f h on %d workers",
                el, length(SD_FIELD) * R_REP * el / N_WORKERS / 3600, N_WORKERS))
print(probe, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------------
## Grid
## ---------------------------------------------------------------------
plan(multisession, workers = N_WORKERS)
res <- list()
for (k in seq_along(SD_FIELD)) {
  s <- SD_FIELD[k]
  message(sprintf("[%d/%d] sigma = %.3f", k, length(SD_FIELD), s))
  out <- future_lapply(seq_len(R_REP), \(r) one_rep(r, s), future.seed = TRUE)
  res[[k]] <- do.call(rbind, Filter(Negate(is.null), out))
  saveRDS(res, file.path(OUTDIR, "raw.rds"))
}
plan(sequential)
dat <- do.call(rbind, res)
if (is.null(dat) || !nrow(dat))
  stop("every replicate failed. Run one_rep(1, 0.5) at the console to see ",
       "the error; a sequential success with parallel failure usually means ",
       "the workers are missing something init_worker() provides.")
message(sprintf("%d of %d replicates returned",
                nrow(dat), length(SD_FIELD) * R_REP))

## ---------------------------------------------------------------------
## Summary against the prediction
## ---------------------------------------------------------------------
agg <- do.call(rbind, lapply(split(dat, dat$sd_field), function(d) {
  data.frame(sd_field = d$sd_field[1], n_ok = nrow(d),
             bias      = mean(d$est) - d$truth[1],
             mc_se     = sd(d$est) / sqrt(nrow(d)),
             coverage  = mean(d$covered),
             sigma_hat = mean(d$sigma_hat, na.rm = TRUE))
}))
B1 <- dat$truth[1]
agg$bias_if_none_absorbed <- B1 * ((1 + 0.346 * agg$sd_field^2)^(-0.5) - 1)
agg$frac_of_bound <- agg$bias / agg$bias_if_none_absorbed

print(agg, row.names = FALSE, digits = 3)
write.csv(agg, file.path(OUTDIR, "bias_scaling.csv"), row.names = FALSE)

cat("\n--- how to read this ---\n")
cat("sigma = 0.001 is the control. If the bias is not ~0 there, the field\n")
cat("is not the mechanism and this hypothesis is wrong.\n")
cat("If the bias grows with sigma and frac_of_bound stays in (0, 1), the\n")
cat("attenuation explanation holds and Appendix S1.1 can say so.\n")
cat("The study's own simulation used sigma = 0.5; compare that row against\n")
cat("the -0.024 reported there.\n")
cat("\nWatch coverage. If it stays near nominal as the bias grows, the bias\n")
cat("is small relative to posterior spread, which is the point the\n")
cat("manuscript should make whatever the mechanism turns out to be.\n")
