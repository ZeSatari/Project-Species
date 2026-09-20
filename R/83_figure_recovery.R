## =====================================================================
## 83_figure_recovery.R -- bias and coverage, M1 against M2
## =====================================================================
## The one figure of Paper A. Table S2 already carries these numbers, but
## the paper's claim is a pattern across them: adding the biotic term to
## the detection component moves alpha_0 sharply and leaves everything
## else where it was. A reader has to hold twenty-eight cells in mind to
## see that in the table; the figure shows it at once.
##
## Two panels, sharing an x axis of parameters:
##   (a) bias, with a reference line at zero
##   (b) coverage, with a reference line at the nominal 0.95
## Colour distinguishes the two models; regimes are facets.
##
## Reads results/single/TableS2_recovery.csv. Refits nothing.
## =====================================================================

source("R/00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggh4x)     # per-facet y scales
})

## 00_setup.R attaches MASS, whose select() masks dplyr's, and stats::filter
## masks dplyr::filter. Every verb below is therefore qualified: the failure
## mode is a confusing 'unused argument' error rather than a wrong answer,
## but it is avoidable either way.

OUTDIR <- "results/single"
FIGDIR <- "figures"
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

rec <- read.csv(file.path(OUTDIR, "TableS2_recovery.csv"),
                stringsAsFactors = FALSE)

## The confounded regime was run separately and only under M2, so it is
## not in TableS2_recovery.csv. It belongs in the figure nonetheless: it
## is the regime in which Assumption 1 fails, and leaving it out would
## drop the one setting where the construction is under strain. It is
## drawn as a fourth column with a single series, and the caption says
## why there is no comparison there.
f_conf <- file.path(OUTDIR, "confounded_recovery.csv")
if (file.exists(f_conf)) {
  conf <- read.csv(f_conf, stringsAsFactors = FALSE)
  if (!"scenario" %in% names(conf)) conf$scenario <- "confounded"
  if (!"model"    %in% names(conf)) conf$model    <- "M2"
  keep <- intersect(names(rec), names(conf))
  rec  <- rbind(rec[, keep], conf[, keep])
  message("confounded regime added: ", nrow(conf), " rows")
} else {
  warning("confounded_recovery.csv not found; the figure will show three regimes")
}

## Print the structure once. The column names below are asserted rather
## than assumed; if your file differs the assertion fails loudly instead
## of the figure being silently wrong.
cat("columns:", paste(names(rec), collapse = ", "), "\n")
print(head(rec, 3))

NEEDED <- c("scenario", "model", "parameter", "bias", "mc_se",
            "coverage", "cov_se")
missing <- setdiff(NEEDED, names(rec))
if (length(missing))
  stop("TableS2_recovery.csv lacks: ", paste(missing, collapse = ", "))

## ---------------------------------------------------------------------
## The file is already one row per scenario x model x parameter, with the
## Monte Carlo standard error beside each quantity. Both are carried
## through, so the figure shows how well each point is itself determined
## rather than implying that 200 replicates pin it exactly.
## ---------------------------------------------------------------------
long <- rec |>
  dplyr::select(dplyr::all_of(NEEDED)) |>
  tidyr::pivot_longer(
    cols = c(bias, mc_se, coverage, cov_se),
    names_to = "col", values_to = "v") |>
  dplyr::mutate(quantity = ifelse(col %in% c("bias", "mc_se"), "Bias", "Coverage"),
         part     = ifelse(col %in% c("bias", "coverage"), "value", "se")) |>
  dplyr::select(-col) |>
  tidyr::pivot_wider(names_from = part, values_from = v) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::mutate(
    model    = factor(model, levels = c("M1", "M2")),
    quantity = factor(quantity, levels = c("Bias", "Coverage")),
    regime   = factor(tools::toTitleCase(scenario),
                      levels = c("Independent", "Shared", "Observer",
                                 "Confounded")),
    parameter = factor(parameter,
                       levels = c("beta0", "beta1", "alpha0", "alpha1", "alpha2"),
                       labels = c("beta[0]", "beta[1]",
                                  "alpha[0]", "alpha[1]", "alpha[2]"))
  ) |>
  dplyr::filter(!is.na(regime), !is.na(parameter))

stopifnot(nrow(long) > 0,
          all(levels(long$model) == c("M1", "M2")),
          all(c("Bias", "Coverage") %in% levels(long$quantity)))
cat("\nrows plotted:", nrow(long), "\n")

## Reference lines differ by panel, so they are supplied as data rather
## than as two geom_hline calls with hard-coded facet membership.
refs <- data.frame(quantity = factor(c("Bias", "Coverage"),
                                     levels = c("Bias", "Coverage")),
                   yint = c(0, 0.95))

theme_paper <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text       = element_text(size = 8),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        axis.title.x     = element_blank())

fig <- ggplot(long, aes(parameter, value,
                        colour = model, shape = model, group = model)) +
  geom_hline(data = refs, aes(yintercept = yint),
             linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  ## preserve = "single" keeps the M2 offset identical in the confounded
  ## facet, which has no M1 series; without it dodge centres the lone
  ## series and the x positions no longer line up across facets.
  geom_linerange(aes(ymin = value - se, ymax = value + se),
                 position = position_dodge(width = 0.45, preserve = "single"),
                 linewidth = 0.4, na.rm = TRUE) +
  geom_point(size = 1.9, position = position_dodge(width = 0.45, preserve = "single")) +
  facet_grid(quantity ~ regime, scales = "free_y", switch = "y") +
  scale_x_discrete(labels = scales::parse_format()) +
  ## The coverage panel exists to show distance from the nominal 0.95, so
  ## 0.95 must be a labelled break; the default breaks skip it and the
  ## reference line then has no readable position.
  ggh4x::facetted_pos_scales(y = list(
    quantity == "Bias"     ~ scale_y_continuous(),
    quantity == "Coverage" ~ scale_y_continuous(
      breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
      labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "0.95")))) +
  scale_colour_manual(values = c(M1 = "#B2182B", M2 = "#2166AC"),
                      labels = c(M1 = "M1 (baseline)",
                                 M2 = "M2 (detection-informed)")) +
  scale_shape_manual(values = c(M1 = 16, M2 = 17),
                     labels = c(M1 = "M1 (baseline)",
                                M2 = "M2 (detection-informed)")) +
  labs(y = NULL) +
  theme_paper +
  theme(strip.placement = "outside")

ggsave(file.path(FIGDIR, "figure_recovery.pdf"), fig, width = 6.5, height = 4.2)
ggsave(file.path(FIGDIR, "figure_recovery.png"), fig, width = 6.5, height = 4.2,
       dpi = 300)
message("Wrote ", file.path(FIGDIR, "figure_recovery.pdf"))

## ---------------------------------------------------------------------
## Check the figure against the claims the text makes about it
## ---------------------------------------------------------------------
cov <- long |> dplyr::filter(quantity == "Coverage", parameter == "alpha[0]")
cat("\nalpha_0 coverage  M1:",
    paste(sprintf("%.3f", cov$value[cov$model == "M1"]), collapse = " "),
    "| M2:",
    paste(sprintf("%.3f", cov$value[cov$model == "M2"]), collapse = " "), "\n")
d <- long |>
  dplyr::filter(quantity == "Bias", parameter %in% c("beta[0]", "beta[1]")) |>
  dplyr::select(regime, parameter, model, value) |>
  tidyr::pivot_wider(names_from = model, values_from = value) |>
  dplyr::mutate(diff = abs(M2 - M1))
cat("largest M1-M2 difference in occupancy bias:",
    sprintf("%.4f", max(d$diff, na.rm = TRUE)), "(text says 0.012)\n")
