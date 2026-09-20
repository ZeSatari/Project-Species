## =====================================================================
## 51_paired_runs_a_c.R -- summaries of runs (a), (b), (c) from stored output
## =====================================================================
## PaperB, Section 4 and Appendix S2. No model is refitted.
##
## Inputs (results/msom/):
##   msom_replicates.rds      run (a), v2
##   msom_replicates_v3.rds   run (b), v3   (optional)
##   msom_replicates_v5.rds   run (c), v5
## Each is list(comm, species, shrink); species has one row per
## species-replicate with lwr, upr, truth, covered, rep.
##
## Produces:
##   1. paired comparison of the community parameters, (c) minus (a):
##      mean difference in error with its paired SE, McNemar's test for
##      coverage, and mean CRPS and interval width per run;
##   2. run (b): proportion of species-replicates closer to truth;
##   3. species-level RMSE, proportion closer and coverage of
##      beta0 + gamma_i by tercile of |beta0 + gamma_i - beta0|, the
##      generating distance from the community mean, for runs (a) and (c).
## Runs (a) and (c) are paired by replicate (identical generating
## intercepts); run (b) is not (see PaperB, Section 4.1).
##
## NOT RUN BY THE AUTHOR OF THIS DRAFT. Column names of the `species`
## element are checked on load; adjust SPECIES_KEY if yours differ.
## =====================================================================

OUTDIR <- "results/msom"
BETA0  <- -0.5
SPECIES_KEY <- "species"   # column identifying the species in `species`

rd <- function(f) {
  p <- file.path(OUTDIR, f)
  if (!file.exists(p)) return(NULL)
  readRDS(p)
}
A <- rd("msom_replicates.rds")
B <- rd("msom_replicates_v3.rds")
C <- rd("msom_replicates_v5.rds")
stopifnot(!is.null(A), !is.null(C))

## ---------------------------------------------------------------------
## 1. Community parameters, paired
## ---------------------------------------------------------------------
cm <- merge(A$comm, C$comm, by = c("rep", "parameter"), suffixes = c("_a", "_c"))
stopifnot(all(cm$truth_a == cm$truth_c))

comm_tab <- do.call(rbind, lapply(split(cm, cm$parameter), function(d) {
  err_a <- d$est_a - d$truth_a
  err_c <- d$est_c - d$truth_c
  dd    <- err_c - err_a
  ## McNemar on the 2x2 table of coverage indicators
  tb <- table(factor(d$covered_a, levels = 0:1), factor(d$covered_c, levels = 0:1))
  mc <- if (tb[1, 2] + tb[2, 1] > 0) mcnemar.test(tb, correct = TRUE)$p.value else NA_real_
  data.frame(
    parameter      = d$parameter[1],
    n_pairs        = nrow(d),
    bias_a         = mean(err_a),
    bias_c         = mean(err_c),
    diff_bias      = mean(dd),
    se_diff_paired = sd(dd) / sqrt(nrow(d)),
    cov_a          = mean(d$covered_a),
    cov_c          = mean(d$covered_c),
    discordant_a1_c0 = tb[2, 1],
    discordant_a0_c1 = tb[1, 2],
    mcnemar_p      = mc,
    crps_a         = mean(d$crps_a),
    crps_c         = mean(d$crps_c),
    width_a        = mean(d$ci_width_a),
    width_c        = mean(d$ci_width_c),
    post_sd_a      = mean(d$post_sd_a),
    post_sd_c      = mean(d$post_sd_c)
  )
}))
cat("\n=== 1. Community parameters, paired (c) minus (a) ===\n")
print(comm_tab, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------------
## 2. Run (b): proportion closer to truth (not paired with (a))
## ---------------------------------------------------------------------
if (!is.null(B)) {
  cat("\n=== 2. Run (b) ===\n")
  cat(sprintf("closer to truth: %.3f | RMSE hier %.3f | RMSE indep %.3f | n = %d\n",
              mean(B$shrink$msom_closer),
              sqrt(mean((B$shrink$msom - B$shrink$truth)^2)),
              sqrt(mean((B$shrink$ssom - B$shrink$truth)^2)),
              nrow(B$shrink)))
  cat("community (b):\n")
  print(aggregate(cbind(crps, ci_width, covered) ~ parameter, B$comm, mean), digits = 3)
} else {
  cat("\nRun (b) output not found; section 2 skipped.\n")
}

## ---------------------------------------------------------------------
## 3. By tercile of the generating distance from the community mean
## ---------------------------------------------------------------------
if (!SPECIES_KEY %in% names(A$species)) {
  message("Column '", SPECIES_KEY, "' not in A$species; columns are: ",
          paste(names(A$species), collapse = ", "),
          ". Falling back to row order within replicate.")
  add_key <- function(sp) {
    sp[[SPECIES_KEY]] <- ave(sp$rep, sp$rep, FUN = seq_along); sp
  }
  A$species <- add_key(A$species); C$species <- add_key(C$species)
}

## tercile breaks from run (a); identical truths in (c)
dist_a <- abs(A$shrink$truth - BETA0)
brks   <- quantile(dist_a, c(0, 1/3, 2/3, 1))
brks[1] <- -Inf; brks[4] <- Inf
lab <- c("near mean", "middle", "far from mean")

by_tercile <- function(R, run) {
  sh <- R$shrink
  sh$terc <- cut(abs(sh$truth - BETA0), brks, labels = lab)
  sh$species <- as.character(sh$species)
  R$species[[SPECIES_KEY]] <- as.character(R$species[[SPECIES_KEY]])
  sp <- merge(sh[, c("rep", "species", "terc")],
              R$species[, c("rep", SPECIES_KEY, "covered")],
              by.x = c("rep", "species"), by.y = c("rep", SPECIES_KEY))
  stopifnot(nrow(sp) == nrow(sh))
  out <- do.call(rbind, lapply(split(sh, sh$terc), function(d) data.frame(
    run        = run,
    tercile    = as.character(d$terc[1]),
    n          = nrow(d),
    rmse_indep = sqrt(mean((d$ssom - d$truth)^2)),
    rmse_hier  = sqrt(mean((d$msom - d$truth)^2)),
    closer     = mean(d$msom_closer)
  )))
  cov <- tapply(sp$covered, sp$terc, mean)
  out$coverage_hier <- as.numeric(cov[out$tercile])
  out
}
terc_tab <- rbind(by_tercile(A, "a"), by_tercile(C, "c"))
cat("\n=== 3. Species-level, by tercile of |generating intercept - community mean| ===\n")
cat(sprintf("breaks: %.3f, %.3f\n", brks[2], brks[3]))
print(terc_tab, row.names = FALSE, digits = 3)
cat(sprintf("\nOverall coverage of beta0 + gamma_i (hierarchical): (a) %.3f  (c) %.3f\n",
            mean(A$species$covered), mean(C$species$covered)))

write.csv(comm_tab, file.path(OUTDIR, "paired_community_a_c.csv"), row.names = FALSE)
write.csv(terc_tab, file.path(OUTDIR, "species_by_tercile_a_c.csv"), row.names = FALSE)
