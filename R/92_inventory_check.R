## =====================================================================
## 92_inventory_check.R -- is the repository complete and consistent?
## =====================================================================
## Run from the project root. It checks, and reports rather than fixes:
##
##   1. scripts: every script either paper needs is present, and nothing
##      is present that belongs to neither;
##   2. results: every file the two verification scripts read exists, so
##      that a reader gets 0 SKIPPED rather than a silent gap;
##   3. figures: the files the manuscripts include by name exist, and no
##      orphan figure is left in the folder;
##   4. manuscripts: found, free of \TBD markers, and their
##      \includegraphics targets resolve;
##   5. packaging: what is large enough to need excluding, and whether
##      .gitignore covers it.
##
## Nothing here refits or resamples; it takes seconds.
## =====================================================================

ok  <- function(...) cat("  OK   ", ..., "\n", sep = "")
bad <- function(...) cat("  MISS ", ..., "\n", sep = "")
hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")
problems <- 0L
note <- function(cond, msg) {
  if (cond) ok(msg) else { bad(msg); problems <<- problems + 1L }
}

## --- 1. scripts ------------------------------------------------------
shared <- c("00_setup.R", "10_metrics.R", "20_simulate_single.R",
            "30_fit_occupancy.R")
paperA <- c("40_run_single_sim.R", "60_real_data.R", "61_stability_check.R",
            "62_missing_sensitivity.R", "64_mesh_adequacy.R", "65_site_cv.R",
            "66_prior_sensitivity.R", "67_sigma_gamma_check.R",
            "68_bias_source.R", "69_bias_field.R", "70_strategy_check.R",
            "71_bias_field_scaling.R", "72_mcmc_benchmark.R",
            "73_prior_sensitivity_field.R", "82_slim_outputs.R",
            "83_figure_recovery.R", "90_verify_claims.R",
            "99_diag_confounding.R", "99_diag_group.R")
paperB <- c("50_run_msom_sim_v2.R", "50_run_msom_sim_v3.R",
            "50_run_msom_sim_v4.R", "50_run_msom_sim_v5.R",
            "51_paired_runs_a_c.R", "52_remaining_checks.R",
            "53_detection_heterogeneity_sim.R", "54_beta1_mechanism.R",
            "55_bias_and_coverage.R", "63_avg_occupancy.R",
            "74_strategy_check_msom.R", "75_field_prior_msom.R",
            "76_detection_heterogeneity.R", "77_detection_gap_check.R",
            "78_detection_overdispersion.R", "79_detection_ladder.R",
            "80_figure_species.R", "81_figure_field.R", "84_figures_B.R",
            "85_mcmc_msom.R", "90_verify_claims_B.R", "91_text_audit.R")
expected <- c(shared, paperA, paperB)

hdr("1. Scripts")
here <- list.files("R", pattern = "\\.R$")
for (f in setdiff(expected, here)) { bad("absent: R/", f); problems <- problems + 1L }
extra <- setdiff(here, c(expected, "92_inventory_check.R"))
for (f in extra) cat("  ?    unaccounted for: R/", f, "\n", sep = "")
note(length(setdiff(expected, here)) == 0,
     sprintf("%d of %d expected scripts present", sum(expected %in% here), length(expected)))

## --- 2. results the verification scripts read ------------------------
hdr("2. Results the verification scripts read")
needed <- c(
  "results/claim_verification.csv",
  "results/real/avg_occupancy.csv",
  "results/real/detection_heterogeneity.csv",
  "results/real/detection_gap_check.csv",
  "results/real/detection_overdispersion.csv",
  "results/real/remaining_checks.rds",
  "results/real/field_summary_2018.csv",
  "results/cv/site_cv_raw.rds",
  "results/msom/msom_replicates.rds",
  "results/msom/msom_replicates_v3.rds",
  "results/msom/msom_replicates_v5.rds",
  "results/msom/species_by_tercile_a_c.csv",
  "results/ladder/ladder_full_summary.csv",
  "results/ladder/ladder_cv_overall.csv",
  "results/prior_sens/field_prior_sensitivity.csv",
  "results/strategy/strategy_comparison.csv",
  "results/het_sim/het_sim.rds",
  "results/beta1/beta1_mechanism.csv",
  "results/beta1/bias_and_coverage.csv")
for (p in needed) if (file.exists(p)) ok(p) else { bad(p); problems <- problems + 1L }

## the MCMC chains, in whichever form they survive
mc <- list.files("results/mcmc", pattern = "tMsPGOcc.*\\.(rds|csv)$",
                 full.names = TRUE)
note(length(mc) > 0, sprintf("MCMC reference output present (%d file(s))", length(mc)))

## --- 3. figures ------------------------------------------------------
hdr("3. Figures")
figs <- c("figures/figure_field.pdf", "figures/figure_species_M3.pdf",
          "figures/figure_species_M1.pdf", "figures/figure_recovery.pdf")
for (p in figs) if (file.exists(p)) ok(p) else { bad(p); problems <- problems + 1L }
orphan <- setdiff(list.files("figures", full.names = TRUE), figs)
for (p in orphan) cat("  ?    not used by either manuscript: ", p, "\n", sep = "")

## --- 4. manuscripts --------------------------------------------------
hdr("4. Manuscripts")
tex <- list.files("..", pattern = "^Paper[AB]\\.tex$",
                  recursive = TRUE, full.names = TRUE)
if (!length(tex)) {
  bad("no PaperA.tex or PaperB.tex found under the parent directory")
  problems <- problems + 1L
}
for (p in tex) {
  tx <- readLines(p, warn = FALSE)
  ## count uses of the marker, not the line that defines the macro
  n_tbd <- sum(grepl("\\\\TBD\\{", tx) & !grepl("newcommand", tx))
  gs <- regmatches(tx, regexpr("figures/[A-Za-z0-9_]+\\.pdf", tx))
  miss <- gs[!file.exists(gs)]
  cat("  ", basename(p), " (", dirname(p), ")\n", sep = "")
  note(n_tbd == 0, sprintf("    %d TBD marker(s) in the text", n_tbd))
  if (!any(grepl("newcommand.*TBD", tx)) && n_tbd > 0)
    cat("  NOTE  \\TBD is used but not defined; the file will not compile\n")
  note(length(miss) == 0,
       sprintf("    %d of %d included figures resolve", length(gs) - length(miss),
               length(gs)))
}
if (length(tex) > length(unique(basename(tex))))
  cat("  ?    more than one copy of a manuscript exists; edits can diverge\n")

## --- 5. packaging ----------------------------------------------------
hdr("5. Packaging")
all_f <- list.files(c("R", "results", "figures"), recursive = TRUE,
                    full.names = TRUE)
big <- all_f[file.size(all_f) > 20e6]
cat(sprintf("  package size excluding files over 20 MB: %.1f MB\n",
            sum(file.size(setdiff(all_f, big))) / 1e6))
if (length(big)) {
  cat("  large files, which should not be uploaded:\n")
  for (p in big) cat(sprintf("     %6.0f MB  %s\n", file.size(p) / 1e6, p))
  if (file.exists(".gitignore")) {
    for (p in big) {
      covered <- suppressWarnings(
        system2("git", c("check-ignore", "-q", shQuote(p)))) == 0
      if (!covered) {
        cat("  NOTE  .gitignore may not cover ", p,
            " (git reports otherwise, or git is unavailable)\n", sep = "")
      }
    }
  } else {
    cat("  NOTE  no .gitignore in the project root\n")
  }
}

hdr("Summary")
if (problems == 0) cat("  nothing missing.\n") else
  cat("  ", problems, " item(s) missing; see above.\n", sep = "")
cat("  Entries marked ? are not errors: decide whether each belongs.\n")
