## =====================================================================
## 82_slim_outputs.R -- shippable artifacts from fits.rds
## =====================================================================
## results/real/fits.rds is ~2.9 GB and cannot be uploaded. It is large
## because it holds four things at once: the MSOM fit with
## config = TRUE (which retains a precision matrix and mean vector for
## every CCD configuration), twelve single-species fits, and M1 and M2.
##
## Nothing downstream needs all of that. This script extracts what is
## actually read and writes it to small files, so the repository can
## carry the derived results without the fitted objects.
##
##   What reads fits.rds now        What it needs
##   ---------------------------    ----------------------------------
##   90_verify_claims.R             ssom[[i]]$summary.fixed only
##   80_figure_species.R            joint posterior draws of the MSOM
##   81_figure_field.R              joint posterior draws of the MSOM
##
## The draws are the whole point of config = TRUE, but the figures use
## only a handful of latent components: the intercept, the elevation
## coefficients, the twelve species effects and the field at the target
## year. That is about 300 numbers per draw, against a latent vector of
## several thousand, so storing the draws costs a few megabytes.
##
## TRADE-OFF, stated plainly. After this, a reviewer can reproduce every
## reported number and both figures, but cannot derive a quantity nobody
## anticipated -- a different year, or a posterior summary not extracted
## here -- without rerunning 60_real_data.R. That seems the right side of
## the trade at 2.9 GB, but it is a loss and the README should say so.
## =====================================================================

source("R/00_setup.R")
suppressPackageStartupMessages(library(spOccupancy))
data(hbefTrends)

OUTDIR <- "results/real"
fits   <- readRDS(file.path(OUTDIR, "fits.rds"))

sz <- function(f) sprintf("%.1f MB", file.size(f) / 1e6)
message("fits.rds: ", sz(file.path(OUTDIR, "fits.rds")))

## ---------------------------------------------------------------------
## 1. Single-species fits -> summary tables only
## ---------------------------------------------------------------------
ssom_summaries <- lapply(fits$ssom, function(f) {
  if (inherits(f, "ssom_failed")) return(NULL)
  list(fixed = f$summary.fixed,
       hyper = f$summary.hyperpar,
       waic  = if (!is.null(f$waic)) f$waic$waic else NA_real_,
       mlik  = if (!is.null(f$mlik)) f$mlik[1, 1] else NA_real_)
})
names(ssom_summaries) <- names(fits$ssom)
saveRDS(ssom_summaries, file.path(OUTDIR, "ssom_summaries.rds"))

## Same for M1 and M2, which nothing downstream currently reads but
## which are cheap to keep and awkward to have discarded.
m_summaries <- lapply(list(m1 = fits$m1, m2 = fits$m2), function(f)
  list(fixed = f$summary.fixed,
       hyper = f$summary.hyperpar,
       waic  = f$waic$waic,
       mlik  = f$mlik[1, 1]))
saveRDS(m_summaries, file.path(OUTDIR, "m1_m2_summaries.rds"))

## ---------------------------------------------------------------------
## 2. MSOM -> the latent components the figures use, for every draw
##
## Drawn once here rather than separately in each figure script, so the
## two figures describe the same posterior sample. The seed is recorded.
## ---------------------------------------------------------------------
N_DRAWS   <- 1000
DRAW_SEED <- 20260913

## set.seed() does NOT control inla.posterior.sample(): the sampler uses
## INLA's own generator, and the seed must be passed to the function. The
## earlier figure scripts called set.seed() and sampled without it, so
## each run drew a different sample -- which is why two quantities in
## Section 3.2.2 that are computed from the draws moved when the figures
## were regenerated. A non-zero seed requires single-threaded operation,
## which 00_setup.R already sets (num.threads = "1:1").
stopifnot(identical(OPTS$num_threads, "1:1"))
set.seed(DRAW_SEED)

message("Drawing ", N_DRAWS, " joint posterior samples (seed ",
        DRAW_SEED, ")...")
samp <- INLA::inla.posterior.sample(n = N_DRAWS, result = fits$msom,
                                    seed = DRAW_SEED)
nm   <- rownames(samp[[1]]$latent)

N_SP   <- length(dimnames(hbefTrends$y)[[1]])
N_YEAR <- dim(hbefTrends$y)[3]

idx_of <- function(pattern, expected = 1L) {
  i <- grep(pattern, nm)
  if (length(i) != expected)
    stop(sprintf("idx_of('%s'): matched %d, expected %d",
                 pattern, length(i), expected))
  i
}

i_int  <- idx_of("^Int_occ:")
i_ele  <- idx_of("^elev_s:")
i_ele2 <- idx_of("^elev_s2:")
i_sp   <- vapply(seq_len(N_SP),
                 \(k) idx_of(sprintf("^species_id:%d$", k)), integer(1))

i_field_all <- grep("^spatialfield:", nm)
n_spde <- length(i_field_all) / N_YEAR
stopifnot(n_spde == floor(n_spde))

## All years are kept. The field is n_spde x N_YEAR per draw, which at
## 283 nodes and nine years is 2547 numbers -- still small, and it means
## a later question about a different year does not require the 2.9 GB
## object back.
draws <- list(
  seed    = DRAW_SEED,
  n_draws = N_DRAWS,
  n_spde  = n_spde,
  n_year  = N_YEAR,
  species = dimnames(hbefTrends$y)[[1]],
  Int_occ = vapply(samp, \(d) d$latent[i_int],  numeric(1)),
  elev_s  = vapply(samp, \(d) d$latent[i_ele],  numeric(1)),
  elev_s2 = vapply(samp, \(d) d$latent[i_ele2], numeric(1)),
  gamma   = vapply(samp, \(d) d$latent[i_sp],   numeric(N_SP)),
  field   = vapply(samp, \(d) d$latent[i_field_all],
                   numeric(length(i_field_all)))
)
saveRDS(draws, file.path(OUTDIR, "msom_draws.rds"))

## MSOM summaries too, for anything that wants them without sampling.
saveRDS(list(fixed = fits$msom$summary.fixed,
             hyper = fits$msom$summary.hyperpar,
             waic  = fits$msom$waic$waic,
             mlik  = fits$msom$mlik[1, 1]),
        file.path(OUTDIR, "msom_summaries.rds"))

## ---------------------------------------------------------------------
## 3. Report
## ---------------------------------------------------------------------
out <- c("ssom_summaries.rds", "m1_m2_summaries.rds",
         "msom_draws.rds", "msom_summaries.rds")
cat("\n--- written ---\n")
for (f in out) cat(sprintf("  %-24s %s\n", f, sz(file.path(OUTDIR, f))))
cat(sprintf("\ntotal replacing a %s object\n",
            sz(file.path(OUTDIR, "fits.rds"))))
cat("\nfits.rds itself is NOT deleted here. Check the sizes above, run\n")
cat("the figure scripts and 90_verify_claims.R against the slim files,\n")
cat("and only then decide whether to exclude it from the repository.\n")
