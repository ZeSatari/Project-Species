## =====================================================================
## 81_figure_field.R -- the shared spatio-temporal field in 2018
## =====================================================================
## Reads results/real/msom_draws.rds; does not refit or resample.
##
## Since the manuscripts were split this figure belongs to Paper B, where
## it is Figure 1; 84_figures_B.R draws the same figure from the fitted
## objects, while this script needs only the slim draws and so runs on
## what the repository ships.
##
## Why this figure exists. The species panels state that broad spatial
## gradients recur across the twelve species and that
## these are the shared field. On the probability scale that claim is
## not visible: with mean occupancy running from 0.028 to 0.837, the six
## sparsest species compress into the bottom few percent of a common
## 0-1 colour scale and their panels are uniformly dark. The claim is
## therefore asserted rather than shown.
##
## This figure shows omega(s, 2018) directly, which is the quantity the
## claim is about.
##
## Three choices, each of which could have gone the other way.
##
##   Projected onto the 373 survey locations, not interpolated onto a
##   grid. A grid surface looks better but between sites it is largely
##   prior, and Section 5.1 of Paper B already reports that the practical range
##   sits near the domain diameter. Showing only where there are data
##   keeps the figure from implying knowledge the data do not support.
##
##   Posterior mean AND posterior standard deviation, side by side. A
##   mean surface on its own would show structure without saying whether
##   it is determined. If the SD is comparable to the range of the mean,
##   the "recurring gradient" is not established and the manuscript
##   sentence should be weakened rather than illustrated.
##
##   A diverging scale centred at zero, with symmetric limits. The field
##   carries a sum-to-zero constraint (constr = TRUE in the SPDE), so
##   zero is meaningful and the two signs should be equally readable.
## =====================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(spOccupancy)
  library(patchwork)
})

data(hbefTrends)
OUTDIR <- "results/real"
FIGDIR <- "figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

N_SITE <- dim(hbefTrends$y)[2]
N_YEAR <- dim(hbefTrends$y)[3]
YEARS  <- dimnames(hbefTrends$y)[[3]]
stopifnot(length(YEARS) == N_YEAR, !is.null(YEARS))

TARGET_YEAR_INDEX <- N_YEAR
TARGET_YEAR       <- YEARS[TARGET_YEAR_INDEX]
stopifnot(TARGET_YEAR == "2018")          # the year named in the manuscript

coords <- hbefTrends$coords / 1000

## ---------------------------------------------------------------------
## Joint posterior draws. Same seed and draw count as 80_figure1.R so the
## two figures describe the same posterior sample rather than two
## independent ones, which would differ slightly and invite the question
## of why.
## ---------------------------------------------------------------------
## Read the pre-extracted draws (82_slim_outputs.R) rather than sampling
## here, so this figure and Figure 2 describe the same posterior sample
## and the 2.9 GB fitted object need not be shipped.
draws   <- readRDS(file.path(OUTDIR, "msom_draws.rds"))
N_DRAWS <- draws$n_draws
n_spde  <- draws$n_spde
stopifnot(draws$n_year == N_YEAR)
message("Loaded ", N_DRAWS, " posterior draws (seed ", draws$seed, ")")

i_field_t <- ((TARGET_YEAR_INDEX - 1) * n_spde + 1):(TARGET_YEAR_INDEX * n_spde)

bnd  <- fmesher::fm_nonconvex_hull(coords, convex = 0.3)
mesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = c(0.35, 1.2),
                            min.angle = 21, offset = c(0.05, 1.5),
                            cutoff = 0.35)
stopifnot(mesh$n == n_spde)
A_site <- INLA::inla.spde.make.A(mesh, loc = as.matrix(coords))

omega <- vapply(seq_len(N_DRAWS),
                \(d) as.numeric(A_site %*% draws$field[i_field_t, d]),
                numeric(N_SITE))

fld <- data.frame(x = coords[, 1], y = coords[, 2],
                  mean = rowMeans(omega),
                  sd   = apply(omega, 1, sd))
fld$snr <- abs(fld$mean) / fld$sd

saveRDS(fld, file.path(OUTDIR, "field_2018.rds"))
write.csv(fld, file.path(OUTDIR, "field_2018.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## Diagnostics. READ THESE BEFORE WRITING THE CAPTION -- they decide
## whether the manuscript sentence is supported or has to be weakened.
## ---------------------------------------------------------------------
cat("\n--- omega(s, 2018) at the 373 survey locations ---\n")
cat(sprintf("posterior mean : %+.3f to %+.3f  (range %.3f)\n",
            min(fld$mean), max(fld$mean), diff(range(fld$mean))))
cat(sprintf("posterior SD   :  %.3f to  %.3f  (median %.3f)\n",
            min(fld$sd), max(fld$sd), median(fld$sd)))
cat(sprintf("|mean|/SD      : median %.2f | %.1f%% of sites above 2\n",
            median(fld$snr), 100 * mean(fld$snr > 2)))
cat(sprintf("range of the mean field, in units of the median SD: %.1f\n",
            diff(range(fld$mean)) / median(fld$sd)))
cat("\nIf the range of the mean is not large relative to the SD, the\n")
cat("gradient is not established and the sentence in Section 5.6 of\n")
cat("Paper B should be weakened rather than illustrated by this figure.\n\n")

## ---------------------------------------------------------------------
## Panels
## ---------------------------------------------------------------------
base_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid  = element_blank(),
        axis.text   = element_text(size = 6),
        axis.title  = element_text(size = 7),
        legend.position   = "bottom",
        legend.key.width  = unit(1.4, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.title = element_text(size = 7))

lim <- max(abs(fld$mean))

p_mean <- ggplot(fld, aes(x, y, colour = mean)) +
  geom_point(size = 1.1) +
  scale_colour_distiller(name = expression(omega(s, 2018)),
                         palette = "RdBu", direction = -1,
                         limits = c(-lim, lim)) +
  coord_equal() +
  labs(subtitle = "(a) Posterior mean",
       x = "Easting (km)", y = "Northing (km)") +
  base_theme

p_sd <- ggplot(fld, aes(x, y, colour = sd)) +
  geom_point(size = 1.1) +
  ## Anchored at the observed range, not at zero. The SD spans only
  ## 0.218-0.333, so a scale starting at zero pushes every site into the
  ## top third of the ramp and the panel reads as uniform -- the same
  ## compression this figure exists to avoid in Figure 1.
  scale_colour_distiller(name = "Posterior SD",
                         palette = "YlOrRd", direction = 1,
                         limits = range(fld$sd)) +
  coord_equal() +
  labs(subtitle = "(b) Posterior standard deviation",
       x = "Easting (km)", y = "Northing (km)") +
  base_theme

fig <- p_mean | p_sd

f_pdf <- file.path(FIGDIR, "figure_field.pdf")
ggsave(f_pdf, fig, width = 7.5, height = 4.2)
message("Wrote ", f_pdf)
if (isTRUE(getOption("figures.png", FALSE)))   # nothing in the papers uses PNG
  ggsave(sub("\\.pdf$", ".png", f_pdf), fig, width = 7.5, height = 4.2, dpi = 300)
