# BetaDanish 0.2.0

## Major new functionality

* **Bayesian inference**: `bayes_betadanish()` provides random-walk
  Metropolis sampling for the Exponentiated Danish submodel and the
  full four-parameter Beta-Danish model with vague Gamma priors.
* **Competing risks rewrite**: `fit_bd_competing()` now uses
  bound-constrained multi-start L-BFGS-B optimization. New
  `cif_compare()` overlays fitted cumulative incidence functions
  against the Aalen-Johansen estimator and reports Gray's test.
* **Structural properties**: closed-form Shannon entropy
  (`bd_entropy_shannon()`), order-statistic densities
  (`bd_order_stat_pdf()`), mean residual life, hazard-shape
  classification, and stress-strength reliability.
* **Diagnostics**: Cox-Snell residual plots for both AFT
  (`plot.bd_aft()`) and cure (`plot.bd_cure()`) fits.
* **Bootstrap confidence intervals** for AFT and cure models.
* **Finite-sample simulation-study runner** for Table 5.5 of the
  underlying thesis.

## Vignettes

Three new vignettes have been added:

* "Bayesian Estimation with BetaDanish"
* "Competing Risks with the Beta-Danish Distribution"
* "Cure Models with the Beta-Danish Distribution"

## Bug fixes

* `summary.bd_aft()` and `summary.bd_cure()` now apply the
  delta-method back-transform so that reported standard errors are
  on the natural parameter scale, not the log scale.
* `report_betadanish()` no longer prints NULL for AIC and BIC.
* `dbetadanish()` log-pdf is now numerically stable in the right
  tail.
* `qbetadanish()` clamps `p` to the unit interval.

## Infrastructure

* Continuous integration via GitHub Actions on four OS/R
  configurations: ubuntu-release, ubuntu-devel, macOS-release, and
  windows-release.
* Test coverage reporting via Codecov.
* Online package website built with pkgdown.
* All `Suggests` packages used via `requireNamespace()` guards at
  the call sites.

# BetaDanish 0.1.0

* First public release.
* Implements the four-parameter Beta-Danish distribution and its
  three-parameter Exponentiated Danish submodel for survival and
  reliability analysis.
* Maximum-likelihood estimation, goodness-of-fit, model comparison,
  and visualization.
* Built-in datasets: remission, carbon_fibres, transplant, aarset,
  leukemia, melanoma, brain_cancer.
