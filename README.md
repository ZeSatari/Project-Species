# Spatio-temporal occupancy models in INLA-SPDE

Analysis code for two manuscripts that share a model-fitting layer:

* **Paper A** — *Observation-level biotic information in spatio-temporal
  occupancy models* (single focal species; `PaperA.tex`).
* **Paper B** — *Hierarchical multi-species spatio-temporal occupancy models
  in INLA-SPDE* (`PaperB.tex`).

They were one manuscript until the multi-species material outgrew it. The
common code is in `R/00_setup.R`, `R/10_metrics.R` and `R/30_fit_occupancy.R`;
everything else belongs to one paper or the other and is listed under the
corresponding heading below.

Each paper has its own verification script, which restates every numerical
claim in the text as a testable expression, recomputes it from the saved
output, and reports any disagreement. Run the relevant one after any rerun:

```r
source("R/90_verify_claims.R")     # Paper A: 69 claims
source("R/90_verify_claims_B.R")   # Paper B: 122 claims, plus a text audit
source("R/91_text_audit.R")        # both papers: are the verified numbers in the text?
source("R/92_inventory_check.R")   # is anything missing from the repository?
```

The Paper B script also greps the manuscript for each verified value, which
catches numbers that were superseded but left in the text, and counts the
`\TBD{}` markers, which must be zero before submission.

The repository is at
<https://anonymous.4open.science/r/Project-Species-D075>.

## For reviewers: what can be checked, and at what cost

The fitted model objects are outputs, not inputs, and are too large to ship
(`results/real/fits.rds` alone is 2.9 GB). Nothing in them cannot be
regenerated from the data and the code, since the INLA fits are
deterministic. The repository is therefore about 32 MB, and three levels of
checking are open, in increasing order of cost.

**1. Check every number in the manuscripts. Minutes, no fitting.** The
verification scripts recompute each claim from the stored summaries, which are
all in the repository, and report any disagreement. This is the level at which
a reader can confirm that the text matches the analysis.

```r
source("R/90_verify_claims.R")     # Paper A
source("R/90_verify_claims_B.R")   # Paper B, and a text audit
```

**2. Redraw the figures. Minutes, no fitting.** `81_figure_field.R` and
`80_figure_species.R` read the 1000 joint posterior draws in
`results/real/msom_draws.rds`, which is small enough to ship, and reproduce
Paper B's Figure 1 and Figure S1. Figure 2 needs a refit under M3 and belongs
to level 3.

**3. Rebuild everything from the data. Hours.** The empirical data are the
`hbefTrends` object in **spOccupancy**, so nothing has to be obtained from us.
`60_real_data.R` refits the models and writes `fits.rds` in about an hour,
after which every other script runs. The simulations take about a day each;
the run order and the timings are in the tables below.

Two things cannot be reproduced bit-for-bit, and the manuscripts say so: the
MCMC reference fits, because the sampler is stochastic (their chains and
convergence diagnostics are kept in `results/mcmc*/**_slim.rds`), and
quantities read off posterior draws as order statistics, which move by about
0.01 between runs of `inla.posterior.sample()` even with its seed fixed.

## Environment

| | |
|---|---|
| R | 4.5.3 |
| R-INLA | 25.10.19 (`https://inla.r-inla-download.org/R/stable`) |
| Packages | `fmesher`, `sf`, `terra`, `dplyr`, `MASS`, `future.apply`, `sn`, `scoringRules`, `spOccupancy`, `coda`, `ggplot2`, `patchwork` |
| Platform | Windows 11, 4 workers, Intel Core Ultra 5 225U, 16 GB |

`sn` is needed by `inla.posterior.sample()`; without it the fits run but no
draws can be taken. `scoringRules` only cross-checks the CRPS implementation.
`spOccupancy` supplies the empirical data and, for Paper B, the MCMC reference
fits. Windows note: `inla.exe` is blocked by Smart App Control on some
installations, with the message "An Application Control policy has blocked this
file"; INLA is then unusable until that setting is changed.

**Fix one INLA version and reproduce every number under it.** The `occupancy`
family is recent enough that behaviour may differ between releases, and results
from two versions must not be mixed in one table. `record_session()` writes the
version alongside each set of results; scripts that share an output directory
pass a `tag` so their session records do not overwrite one another.

## Shared layer

```
R/00_setup.R              packages, priors, options, shared helpers
R/10_metrics.R            bias, RMSE, coverage, CRPS
R/20_simulate_single.R    single-species generator, four regimes
R/30_fit_occupancy.R      the single fitting entry point
```

`00_setup.R` reads the detection-coefficient cap at run time from the
hyperparameter list the occupancy family declares, rather than taking it from
the documentation, and `build_det_X()` errors at the cap and warns within two
of it. The cap counts **coefficients**, not design columns: `build_det_X()`
lays out one column per coefficient and visit, so a three-coefficient detection
model has nine columns with `K = 3`.

---

# Paper A: observation-level biotic information

```
R/70_strategy_check.R     approximation adequacy        -> results/strategy/
R/40_run_single_sim.R     single-species simulation     -> results/single/
R/61_stability_check.R    empirical fit reproducibility
R/60_real_data.R          Hubbard Brook application     -> results/real/
R/82_slim_outputs.R       shippable artifacts from fits.rds
R/90_verify_claims.R      checks the claims
```

Diagnostics, each cited in Paper A:

```
R/62_missing_sensitivity.R  treatment of absent co-occurring records
R/64_mesh_adequacy.R        three mesh resolutions      -> results/mesh/
R/66_prior_sensitivity.R    species-effect prior, four specifications
R/67_sigma_gamma_check.R    reproducibility of sigma_gamma
R/68_bias_source.R          slope bias, likelihood in isolation
R/69_bias_field.R           slope bias, field and covariate form
R/71_bias_field_scaling.R   slope bias vs generating field sigma
R/99_diag_group.R           the `A.local` grouping convention
R/99_diag_confounding.R     covariate/field separability
```

Figure scripts are named for their content, not their number: the shared-field
figure precedes the species panels in the manuscript, so a file called
`figure1.pdf` that renders as Figure 2 would be a trap.

| Script | Fits | Time |
|---|---|---|
| `70_strategy_check` | 3 | ~45 min |
| `40_run_single_sim` | 1600 | ~10 h (4 workers) |
| `60_real_data` | 16 | ~1 h |

---

# Paper B: hierarchical multi-species structure

Scripts written for Paper B occupy 74-79, 84-85 and the 5x block, so that they
do not collide with Paper A's, which use 61-73, 82 and 83. Two numbers in
Paper A's range, `80_figure_species.R` and `81_figure_field.R`, draw Paper B's
figures: they predate the split and were left where they were rather than
renumbered, since the manuscripts name no scripts. `63_avg_occupancy.R`
is a deliberate replacement, noted below. Two of the 5x scripts,
`54_beta1_mechanism.R` and `55_bias_and_coverage.R`, are diagnostics rather
than parts of the reported analysis: they test explanations for the
community-slope displacement and for the species-level interval coverage, and
both are cheap because they strip the model to the mechanism under test.

```
R/50_run_msom_sim_v2.R    multi-species simulation, run (a)   -> results/msom/
R/50_run_msom_sim_v3.R    run (b), homogeneous slopes
R/50_run_msom_sim_v5.R    run (c), homogeneous detection
R/51_paired_runs_a_c.R    paired comparison of (a) and (c); tercile tables
R/52_remaining_checks.R   mesh, field constraint, detection coefficients,
                          predictive check, cluster CV standard errors
R/53_detection_heterogeneity_sim.R  what unmodelled heterogeneity does
R/54_beta1_mechanism.R    is the community slope an information-weighted mean?
R/55_bias_and_coverage.R  where the slope bias enters; are the species-level
                          intervals too narrow, or is it shrinkage?
R/63_avg_occupancy.R      average occupancy and P(recorded)   [see note]
R/74_strategy_check_msom.R    approximation strategies, hierarchical model
R/75_field_prior_msom.R       field-prior sensitivity, hierarchical model
R/76_detection_heterogeneity.R  conditional estimates of per-visit detection
R/77_detection_gap_check.R      fitted vs conditional detection
R/78_detection_overdispersion.R dispersion of the detection counts
R/79_detection_ladder.R         the M1-M2-M3 ladder, with cross-validation
R/85_mcmc_msom.R                MCMC references (stMsPGOcc, tMsPGOcc)
R/80_figure_species.R           Figure S1 (occupancy under M1) from the
                                slim draws, without needing fits.rds
R/81_figure_field.R             Figure 1 (the shared field), likewise
R/84_figures_B.R                all three figures, and the only route to
                                Figure 2, which needs a refit under M3
R/90_verify_claims_B.R          checks the claims and audits the text
R/91_text_audit.R               text audit for both papers
```

**Note on `63_avg_occupancy.R`.** The Paper B version replaces the earlier one.
It seeds `inla.posterior.sample()`, computes P(recorded) over the visits
actually made rather than over all three, samples the hierarchical model once
and reuses the draws across species, and reports P(recorded) beside average
occupancy. The earlier output is kept as `avg_occupancy_unseeded.csv`.

**`50_run_msom_sim_v4.R`** is the miscalibrated heterogeneous-prevalence
attempt. It is kept for the record, is described in Section 4.4 of Paper B, and
is not part of the reported results.

## Order of execution

Later scripts read what earlier ones store:

| Step | Script | Roughly |
|---|---|---|
| 1 | `60_real_data.R` (shared with Paper A) | ~1 h |
| 2 | `63_avg_occupancy.R` | ~10 min |
| 3 | `65_site_cv.R` (shared) | hours |
| 4 | `76_detection_heterogeneity.R` | seconds |
| 5 | `77_detection_gap_check.R` | minutes |
| 6 | `78_detection_overdispersion.R` | seconds |
| 7 | `79_detection_ladder.R` | ~3 h |
| 8 | `84_figures_B.R` | ~10 min |
| 9 | `50_run_msom_sim_v2/v3/v5.R` | ~23 h each |
| 10 | `51_paired_runs_a_c.R` | seconds |
| 11 | `74_strategy_check_msom.R`, `75_field_prior_msom.R` | ~15 min; ~40 min |
| 12 | `53_detection_heterogeneity_sim.R` | minutes |
| 12b | `54_beta1_mechanism.R`, `55_bias_and_coverage.R` | minutes each |
| 13 | `85_mcmc_msom.R` | ~75–110 min per engine |
| 14 | `52_remaining_checks.R` | minutes |
| 15 | `90_verify_claims_B.R` | minutes |
| 16 | `91_text_audit.R` | seconds |

`79_detection_ladder.R` has three stages, switched with `DO_FULL`, `DO_SENS`
and `DO_CV`; run the full-data fits first and inspect them before starting the
cross-validation, which is fifteen fits.

## Where each table and figure comes from

| Item | Script | Stored output |
|---|---|---|
| Table 1, priors | — (`PRIORS` in `00_setup.R`) | — |
| Table 2, community parameters | `50_run_msom_sim_v2/v3/v5.R` | `results/msom/msom_replicates{,_v3,_v5}.rds` |
| Table 3, species-level intercepts by run | `51_paired_runs_a_c.R` | `results/msom/msom_replicates*.rds` |
| Table 4, by tercile | `51_paired_runs_a_c.R` | `results/msom/species_by_tercile_a_c.csv` |
| Table 5, species intercepts | `60_real_data.R` | `results/real/fits.rds` |
| Table 6, detection diagnostic | `76_detection_heterogeneity.R` | `results/real/detection_heterogeneity.csv` |
| Figure 1, the shared field | `81_figure_field.R` or `84_figures_B.R` | `figures/figure_field.pdf`, `results/real/field_summary_2018.csv` |
| Table 7, site-level cross-validation | `65_site_cv.R`; errors from `52_remaining_checks.R` | `results/cv/site_cv_raw.rds` |
| Figure 2, occupancy under M3 | `84_figures_B.R` with `MODEL <- "M3"` | `figures/figure_species_M3.pdf` |
| Table 8, the ladder | `79_detection_ladder.R` | `results/ladder/ladder_full_summary.csv` |
| Table 9, cross-validation of the ladder | `79_detection_ladder.R` | `results/ladder/ladder_cv_overall.csv` |
| Table S1, average occupancy | `63_avg_occupancy.R` | `results/real/avg_occupancy.csv` |
| Table S2, hyperparameters | `52_remaining_checks.R` | `results/real/remaining_checks.rds` |
| Table S3, detection in the independent fits | `77_detection_gap_check.R` | `results/real/detection_gap_check.csv` |
| Table S4, predictive check | `52_remaining_checks.R` | `results/real/remaining_checks.rds` |
| Table S5, heterogeneity simulation | `53_detection_heterogeneity_sim.R` | `results/het_sim/het_sim.rds` |
| Table S6, overdispersion | `78_detection_overdispersion.R` | `results/real/detection_overdispersion.csv` |
| Figure S1, occupancy under M1 | `80_figure_species.R` or `84_figures_B.R` with `MODEL <- "M1"` | `figures/figure_species_M1.pdf` |
| Section 4.2, slope-bias diagnostics | `54_beta1_mechanism.R`, `55_bias_and_coverage.R` | `results/beta1/` |

## What Paper B's scripts need that the slim outputs do not provide

`82_slim_outputs.R` extracts what Paper A's figures and verification read.
Paper B's scripts call `inla.posterior.sample()` themselves and therefore need
the **full** `results/real/fits.rds`, fitted with `config = TRUE`, which is
about 2.9 GB and is not in the repository. Rerun `60_real_data.R` (about an
hour) before running steps 2, 8 or 14 above.

---

# Findings that shape the code

The first five apply to both papers.

**Mesh resolution follows the information the data can support, not the
geometry of the study area.** A mesh of 1724 nodes gave 15,516 field elements
against 3,357 observation rows; repeated fits of the identical model then
returned intercepts spanning 0.183 on the logit scale, run times between 1,900
and 13,100 seconds, and convergence failures in up to 80% of attempts.
Coarsening to 283 nodes reduced the spread to 0.001. `60_real_data.R` warns if
the ratio of field elements to rows exceeds 2. Because the occupancy
log-likelihood is not log-concave, verify that repeated fits reproduce rather
than assuming they do (`61_stability_check.R`).

**`int.strategy = "eb"` is not used.** Empirical Bayes fixes the
hyperparameters at their mode, so intervals understate uncertainty and coverage
is not interpretable. All reported results use `"ccd"`. The strategy-comparison
scripts vary the Laplace `strategy` only; the integration strategy is held
fixed.

**A quadratic covariate attenuates its own coefficient.** `x^2` for
`x ~ N(0,1)` is chi-squared with one degree of freedom, and its right tail
drives psi toward one where the data carry little information about the slope.
With no spatial structure at n = 5000 the estimator recovered 0.697 for a
linear covariate against a truth of 0.700, but 0.456 for raw `x^2` and 0.614
after standardisation. The multi-species simulation therefore uses a linear
covariate; the empirical application retains quadratic elevation.

**The "unused groups" warning is a false positive with `A.local`.** INLA counts
groups from the index vector, which is all `NA` by design. Checked on a
Gaussian analogue with a known AR(1) field: the specification used here recovers
rho = 0.794 for a generating value of 0.80, while the one that suppresses the
warning returns 0.000 (`R/99_diag_group.R`). Do not "fix" it.

**The detection model is capped at ten coefficients**, intercept included.
Raising the cap would not remove the underlying constraint: detection
coefficients enter as hyperparameters rather than latent field elements, and
integration over them becomes infeasible well before ten is reached (Belmont et
al. 2024). Species-specific detection intercepts for twelve species need
fourteen coefficients, which is why Paper B pools detection and then groups it.

The rest were established while preparing Paper B.

**`set.seed()` does not control `inla.posterior.sample()`, and `seed =` does
not make it exactly repeatable either.** The seed must be passed to the
function as `seed =`, which in turn requires `num.threads = "1:1"`; without
that, repeated calls differ enough to move a correlation reported to two
decimal places by 0.06. With it, two calls with the same seed in one session
still gave latent values differing by up to 0.0069. Summaries over 1000 draws
are stable to about 0.001 for means and 0.01 for extremes, so Paper B reports
the affected extremes as approximate. Anything derived from the draws should be
checked by rerunning, not assumed.

**Simulation replicates use L'Ecuyer-CMRG, because the workers do.**
`future_lapply(..., future.seed = TRUE)` sets that generator inside each
worker, and `set.seed(seed)` inside the replicate then seeds it. Reproducing a
replicate in the main session requires `RNGkind("L'Ecuyer-CMRG")` first;
the default Mersenne-Twister gives different data from the same seed. This cost
an afternoon: a replicate regenerated in the main session disagreed with the
stored one on every species intercept.

**`rnorm(n, mu, 0)` returns `mu` without drawing.** Run (b) of the
multi-species simulation generates homogeneous slopes that way, so its
random-number stream diverges from run (a) after the slopes: (a) and (b) share
their fields, covariate and species intercepts but not their sites, detection
parameters or data, and are compared as independent runs. Run (c) instead draws
the detection parameters exactly as (a) does and overwrites them, which leaves
the stream untouched and makes the two runs paired replicate by replicate.
`50_run_msom_sim_v5.R` checks the pairing before the run and again against the
stored output of (a) afterwards.

**The community-slope displacement needs two ingredients, not one.** A
hierarchical logistic model with shrunken random slopes recovers the community
slope without bias (`54_beta1_mechanism.R`), and Paper A reports that a
single-species occupancy model does too; adding the occupancy likelihood to
the hierarchical model brings back about 70% of the displacement
(`55_bias_and_coverage.R`, bias -0.054 against -0.076 in the full model), with
the spatial field accounting for the rest. An information-weighting
explanation was tested there and refuted. The species-level intervals were
checked in the same run: building them from joint posterior samples rather
than from marginal summaries changes neither width nor coverage, so the
under-coverage for species far from the community mean is shrinkage and not an
artefact of how the intervals are formed.

**Detection estimated from the full occupancy likelihood exceeds the estimate
from the recorded site-years alone**, for every species and in both INLA and
MCMC fits, by 0.016 to 0.067. With a free occupancy level the two are the same
estimator by construction, so a difference means the occupancy component is
constrained; simulation (`53_detection_heterogeneity_sim.R`) reproduces the
equality and not the difference. Unexplained, and recorded as such in Paper B.

**Spatial factor loadings are weakly identified on these data.** In the MCMC
reference (`stMsPGOcc`, one factor), chains disagreed on the sign of the
loading for 83% of species and R-hat reached 3; anchoring the fixed loading on
the most frequently recorded species instead of the alphabetically first did
not help. Dropping the field (`tMsPGOcc`) converges for detection but not for
the species-level occupancy coefficients.

# Defects corrected from the earlier pipeline

| Old behaviour | Consequence |
|---|---|
| `mvrnorm(mu = rep(-0.5, 12))` plus `mu_beta0 = -0.5` | true community intercept was −1.0, reported as −0.5 |
| Generated corr +0.5, fitted `Cmatrix = I + J` (corr −1/12) | generation and fitting misspecified in opposite directions |
| `scale_x_s^2` in a formula | `^` is the crossing operator; the quadratic term was never fitted |
| `as.vector(x_covariate)` assigned by position | lengths matched, ordering did not |
| `inla.group.cv(result = model_xs, ...)` | undefined object; a reported cross-validation value came from nowhere |
| CAWA block never subset the species | one table row came from a model that could not run |
| No `constr` on `f(species_id, "iid")` | beta0 and gamma_i not separately identified; CrI 5.3x too wide |
| WAIC/n compared across different response sets | the column tracked prevalence, not model quality |
| "MSE" = squared posterior SD of one fit | measures posterior spread, not distance from truth |
| R = 1 | bias, coverage and RMSE are undefined |
| P(recorded) multiplied over all three visits | overstated for the 18.6% of site-years with a missing visit |
| Ratios computed from rounded table entries | e.g. 0.028/0.023 = 1.22 where the exact ratio is 1.19 |

Two defects introduced during the rewrite were caught by testing:
`species_intercepts()` indexed posterior samples by position rather than name,
attaching species effects to the wrong species; and a `constr = TRUE` random
slope without a companion fixed effect forced the community slope to zero. Both
produced plausible output and neither raised an error.

# Outputs

```
results/strategy/   approximation comparison, timing, marginal overlays
results/single/     per-replicate parameters and fit statistics (Paper A)
results/msom/       community and species-level recovery, runs (a), (b), (c)
results/mesh/       three mesh resolutions
results/cv/         site-held-out cross-validation
results/prior/      species-effect prior sensitivity
results/prior_sens/ field-prior sensitivity, hierarchical model
results/bias/       slope-bias diagnostics
results/ladder/     the M1-M2-M3 ladder
results/het_sim/    detection-heterogeneity simulation
results/beta1/      slope-bias and interval-coverage diagnostics
results/mcmc/       MCMC reference fits and comparisons
results/real/       detection coefficients, species tables, average occupancy,
                    diagnostics, field summaries, verification tables
figures/            figure_field.pdf        Paper B, Figure 1
                    figure_species_M3.pdf   Paper B, Figure 2
                    figure_species_M1.pdf   Paper B, Figure S1
                    figure_recovery.pdf     Paper A, its only figure
```

Each results directory carries `sessionInfo*.txt` and
`priors_and_options*.rds`, suffixed by run tag where several scripts share a
directory.

**Fitted objects are not kept.** `results/real/fits.rds` (2.9 GB) is needed by
every script that samples the posterior and is regenerated by `60_real_data.R`
in about an hour; the ladder and strategy fits are regenerated by
`79_detection_ladder.R` and `74_strategy_check_msom.R` in about three hours and
one hour. The MCMC objects of both papers are instead kept in slimmed form
(`results/mcmc/*_slim.rds`, `results/mcmc_benchmark/*_slim.rds`), holding the
parameter chains and the convergence diagnostics but not the site-level arrays
`psi.samples` and `z.samples`, which is where the size is: a few megabytes
against one to two gigabytes each. They are kept rather than regenerated
because the samplers are stochastic and will not reproduce exactly. None of
this costs verifiability: both verification scripts read the summary tables
rather than the fitted objects, so every reported number survives the
slimming. With the fitted objects excluded, `results/` is about 30 MB.

**Two routes to the figures.** `81_figure_field.R` and `80_figure_species.R`
read the 1000 joint draws in `results/real/msom_draws.rds`, which is small
enough to ship, and reproduce Figure 1 and Figure S1 without the fitted
objects. `84_figures_B.R` reproduces the same two from `fits.rds` and is the
only route to Figure 2, since that needs a refit under M3. A reader who has
only the repository should use the first route. All write PDF. Nothing in either manuscript uses PNG, so the two older
scripts write one only under `options(figures.png = TRUE)`.

# Data

The empirical data are the `hbefTrends` object in `spOccupancy`, from the
Hubbard Brook Ecosystem Study. No data files are redistributed here.

# Not addressed

* **Species-specific loadings** on the shared field are discussed in both
  papers as an extension but not implemented; the MCMC reference suggests they
  are not well identified on nine years of data from twelve species.
* **The source of the slope bias** (−0.024 single-species, −0.050 and −0.076
  multi-species) was not found. `68_bias_source.R` and `69_bias_field.R`
  excluded the likelihood, the approximation, the covariate form and the field;
  `71_bias_field_scaling.R` excluded marginal attenuation induced by the field
  and showed that coverage for beta1 stays within 0.925–0.975 across generating
  field standard deviations from 0 to 1. Paper B adds that pooled detection is
  not the cause either: the bias was −0.076 with heterogeneous detection and
  −0.074 with homogeneous. An information-weighting conjecture is stated in
  Paper B and would need a rerun that stores the species-slope posteriors.
* **The redesigned heterogeneous-prevalence simulation** (run (d)) was not
  carried out. The attempt in `50_run_msom_sim_v4.R` was calibrated from the
  unconditional detection frequency rather than the conditional detection
  probability, which put its rarest species where occupancy is not identifiable
  at all; Section 4.4 of Paper B sets out what a correct design would be.
* **A grouped-detection simulation**, measuring how much of the benefit of
  pooling a coarse grouping recovers when detection is heterogeneous, is
  proposed in Paper B's future directions and not implemented.
