# BetaDanish 0.3.0

## Visualisation

* `bd_ttt_plot()` draws the scaled total time on test transform, a
  distribution-free read on hazard shape before any model is fitted. It
  returns the plotting coordinates and a suggested shape label, so it can
  be compared against `bd_hazard_shape()` from the fitted parameters.

* `bd_profile_plot()` draws a profile log-likelihood with the critical
  threshold and the resulting interval marked. Where the curve never
  returns below the threshold the plot says there is no finite upper
  bound, rather than implying one at the edge of the grid.

* `plot()` method for `bd_bayes` objects, giving trace and posterior
  density panels and restoring the graphical parameters on exit.

* The introduction and competing-risks vignettes cover the new
  functions.

## Competing risks and simulation

* `fit_bd_competing()` accepts covariates, either as a one-sided formula
  with `data` or as a matrix. Each cause carries its own coefficient
  vector acting on the time scale, `T_j = T_0j * exp(x'gamma_j)`, so a
  covariate may accelerate one cause and retard another. The intercept
  is dropped, since each cause already has its own scale parameter.

* `fit_bd_competing(submodel = TRUE)` fits Exponentiated Danish
  cause-specific kernels, and the fit now goes through the same
  degeneracy guard as the univariate model.

* `cif_betadanish()` gains an `x` argument for the covariate values at
  which to evaluate the cumulative incidence, defaulting to the
  reference subject. A covariate effect folds exactly into the scale
  parameter, since scaling time by lambda maps `k` to `k/lambda`.

* **Simulation studies**: `bd_simulation_study()` for the univariate
  model, `bd_simulation_cure()` and `bd_simulation_competing()`,
  alongside `simulate_bd_competing_data()`. Each reports bias, RMSE,
  the ratio of mean estimated standard error to empirical standard
  deviation, and Wald coverage. Replicates with no admissible optimum
  are counted rather than dropped.

## Estimation and inference

* `fit_betadanish(grouped = TRUE)` uses the grouped likelihood for times
  recorded on a coarse grid. An event at `t` contributes
  `log{F(t + delta/2) - F(t - delta/2)}` rather than `log f(t)`. The
  increment is inferred from the spacing of the times unless `delta` is
  supplied, and a warning now fires when the times look grid-recorded
  but the point-density likelihood is being used, which understates
  every standard error.

* `fit_betadanish(penalty = )` adds a ridge penalty on the
  log-parameter scale for the weakly identified regime, shrinking
  toward the unpenalised optimum unless `penalty_center` says
  otherwise.

  The reported `logLik` is the **unpenalised** log-likelihood evaluated
  at the penalised estimate, so AIC, BIC and likelihood ratio tests
  stay comparable with unpenalised fits. The objective actually
  maximised is kept separately as `penalised_logLik`.

* `bd_wald_ci()` builds Wald intervals on the log scale and
  exponentiates, so they cannot cross zero. The symmetric interval on
  the natural scale does exactly that for a weakly identified tail
  index.

* `bd_profile_ci()` inverts the likelihood ratio test for one
  parameter. When the profile stays inside the critical region to the
  top of the grid it returns `Inf` and says so: an unbounded interval
  is a result, and the honest report is a lower bound rather than a
  point estimate.

* `bd_identified_coef()` reports the four-parameter fit as `(ac, b, k)`.
  Near the lower tail `a` and `c` enter almost entirely through their
  product, so the composite is estimated precisely while the two
  separately are not. The log-likelihood is unchanged: this is a
  reparametrisation, not a different model.

## Theoretical properties

* **Moments**: `bd_moments()`, `bd_moment_summary()`,
  `bd_incomplete_moment()` and `bd_conditional_moment()`. Every integral
  is taken on the finite Beta(a, b) scale under the substitution
  `u = G(z)`, so there is no series truncation. The existence condition
  `b > r` is enforced: a moment that does not exist returns `Inf` rather
  than a large finite number from a truncated sum.

* **Ageing and inequality**: `bd_mrl()` and `bd_rmrl()` for mean residual
  life and mean inactivity time, `bd_mean_deviation()`, `bd_lorenz()`,
  `bd_bonferroni()` and `bd_pwm()`.

* **Entropies**: `bd_entropy_shannon()` now uses the closed form, with
  `method = "quadrature"` retained as an independent cross-check.
  `bd_entropy_renyi()` and `bd_entropy_tsallis()` are new.

  The quadrature `bd_entropy_shannon()` that shipped in 0.2.0 from
  `R/entropy.R` is superseded by the closed form and that file is
  removed. The signature changes from
  `(a, b, c, k, subdivisions, rel.tol)` to
  `(a, b, c, k, terms, method, rel.tol, subdivisions)`. Named calls are
  unaffected; a positional call passing `subdivisions` fifth would now
  set `terms`. The old behaviour is available as
  `method = "quadrature"`.

  The closed form's series terms decay like `i^-(b+1)`, so plain
  truncation is accurate for large `b` and not for small `b`. Against
  high-precision integration the untruncated-tail sum at `M = 2000` is
  out by about `1e-8` at `b = 3` but by about `0.11` at `b = 0.5`. An
  analytic Euler-Maclaurin tail is now added, restoring agreement to
  roughly eight digits throughout.

* **Reliability**: `bd_stress_strength()` for `R = P(X > Y)` with `X` the
  strength, evaluated on the Beta scale of the stress variable.

* **Shape**: `bd_hazard_shape()` classifies the hazard as increasing,
  decreasing, bathtub or upside-down bathtub using Glaser's
  `eta(t) = -f'(t)/f(t)`, formed analytically, and reports the mode.

* **Order statistics**: `bd_order_stat_cdf()` and
  `bd_order_stat_moments()`, alongside the existing
  `bd_order_stat_pdf()`.

* **Tail behaviour**: `bd_tail_index()` records that the survival
  function is regularly varying with index `-b`, that `E(Z^r)` is finite
  exactly when `b > r`, that the moment generating function does not
  exist, and that the distribution lies in the Frechet domain of
  attraction.

* **Exponentiated Danish API**: `ded()`, `ped()`, `qed()`, `red()`,
  `sed()` and `hed()` name the `a = 1` submodel directly instead of
  requiring `dbetadanish(x, 1, b, c, k)`.

## Datasets

* **`brain_cancer` has been removed.** The dataset of 500 brain cancer
  patients that shipped in 0.1.0 and 0.2.0 is no longer included, and
  the AFT example in the README now uses `melanoma`. The entry under
  0.1.0 below is left as written, since it records what that release
  actually contained.

* **`guinea_pig` added** (Bjerkedal 1960, n = 72): survival times in days
  of guinea pigs injected with virulent tubercle bacilli. This is the
  reference case for a well-identified four-parameter fit, with
  \eqn{\hat b = 3.64} (SE 1.20) sitting about 2.2 standard errors clear
  of the \eqn{b = 1} ridge. On `remission` and `carbon_fibres` the
  parent model sits on the flat \eqn{(a, c)} direction instead, so
  without this dataset no built-in example showed the full model
  behaving well.

* **`inst/extdata/` added** with four example CSVs -- complete,
  censored, covariate and competing-risks -- so that examples and
  vignettes can demonstrate the file-driven workflow.


## File-driven analysis

* **`bd_analyze_csv()`** reads a delimited or Excel file, fits the
  requested model, assembles tidy result tables, and optionally writes
  those tables and the diagnostic figures to a directory. It covers the
  univariate, AFT, cure and competing-risks paths. Nothing is written to
  disk unless `output_dir` is supplied, and each fit is attempted
  independently so one failure is recorded rather than fatal.

* **`bd_csv_template()`** writes a correctly shaped skeleton CSV in any
  of four layouts, so the expected column structure can be seen rather
  than read about.

* `print`, `summary` and `plot` methods for the new `bd_analysis` class.

* New vignette: "Analysing Your Own Data from a CSV File".
## Data input

* `read_survival_data()` gains `cause_col`, `sep`, `dec`, `encoding`,
  `na.strings`, `drop_na` and `quiet`, and will guess the time and
  status columns when they are not named. Guessing deliberately excludes
  names such as `censor`, whose polarity differs between conventions.

* The `bd_data_report` attribute now also records which columns were

* `read_survival_data(covar_cols = "all")` retains every column that is
  not part of the response. Columns that are dropped are now named in a
  message and recorded in the report, rather than disappearing silently.
  used and whether the times look recorded on a coarse grid.

## Correctness fixes

* **Degenerate optima are no longer reported as fits.** The Beta-Danish
  likelihood has ridges along which the objective can be driven
  arbitrarily high with no finite maximiser: as `k` and `c` grow with
  `c/k` fixed, the family collapses to a Frechet-type limit.
  `optim_multistart()` selected purely on the highest log-likelihood
  behind nothing but a finiteness test, so a runaway beat every honest
  optimum. A four-parameter fit on the `transplant` data returned a
  log-likelihood of `+5.2e77` on 91 observations, and the resulting
  likelihood ratio test reported a chi-squared statistic of `1.0e78`.

  `fit_betadanish()`, `fit_bd_aft()` and `fit_bd_cure()` now search over
  a deterministic start grid before any random starts, and accept only
  admissible optima: shape parameters bounded by 500, scale within
  `1e-8` to `1e8`, and a plausible per-observation log-likelihood. The
  thresholds and the start grid are those already used by the thesis
  analysis pipeline, so the two now agree.

* The fitted object records `starts_ok`, `starts_rejected` and
  `loglik_spread`. Warnings fire when starts are discarded as degenerate
  and when the accepted optima span more than two log-likelihood units,
  which is the signal that `n_starts` is too low for the surface.

* `plot.bd_aft()` and `plot.bd_cure()` now produce the Cox-Snell residual
  plot. Both looked up the shape parameters as `coefficients["b"]` and
  `coefficients["c"]`, but `fit_bd_aft()` and `fit_bd_cure()` store them as
  `log_b` and `log_c`. The lookup returned `NA`, the internal guard caught it,
  and the functions returned without drawing anything. The lookup is fixed and
  the values are exponentiated back to the natural scale.

* `dbetadanish()` is now accurate in the far right tail. The term
  \eqn{\log\{1 - G(t)\}} was formed as `log1p(-exp(log_G))`, which pins at
  about `log(1e-16)` once `G` rounds to one, flooring the log-density near
  -36.8 regardless of its true value. It is now formed as
  `log(-expm1(-c * log1p(1/(k*t))))`, which holds full relative precision for
  every `k*t > 0`.

* `pbetadanish()` computes the survival function through the beta mirror
  identity \eqn{1 - I_y(a,b) = I_{1-y}(b,a)}, so a probability near one is
  never subtracted from one.

* `qbetadanish()` obtains \eqn{1 - u} directly from
  \eqn{1 - q\beta(p; a, b) = q\beta(p; b, a)} with the tail flag reversed. The
  previous route computed `y^(1/c) / (k * (1 - y^(1/c)))`, which loses all
  significant digits as `u` approaches one. Round-trip accuracy now holds to
  `p = 1 - 1e-13`.

* `hbetadanish()` no longer substitutes `-700` for an exhausted log-survival.
  A positive density with zero survival gives `Inf`, an honest divergent
  hazard, and the indeterminate case gives `NaN`.

* The density, distribution, quantile, survival and hazard functions now
  recycle their parameters element-wise against `x`, `q` or `p`, following the
  usual convention for R distribution functions.

* `fit_betadanish()` records `npar`, `nobs`, `nevent`, `AIC` and `BIC` on the
  fitted object, and `report_betadanish()` computes information criteria
  through the `logLik` method. Both previously omitted AIC and BIC silently.

* `summary.bd_aft()` and `summary.bd_cure()` report shape parameters on the
  natural scale with delta-method standard errors and an exponentiated
  log-scale confidence interval, separately from the regression coefficients.
  No Wald test is reported for a shape parameter, since the implied null lies
  outside the parameter space.

* `read_survival_data()` selects columns by name rather than by position, so a
  covariate named `time` or `status` no longer collides with the response.
  Covariates keep their original type. The returned data frame carries a
  `bd_data_report` attribute.

* `extract_surv_data()` validates times and the event indicator up front and
  reports a clear error rather than failing inside the optimiser.

## New

* `sbetadanish()`, an explicit survival function.

* `fit_betadanish()` gains `check_identifiability`. It warns when the optimiser
  reports a poor code, when the information matrix is singular, when
  \eqn{\hat b} is within two standard errors of the \eqn{b = 1}
  non-identifiability ridge, and when the fitted correlation between
  \eqn{\hat a} and \eqn{\hat c} exceeds 0.95 in absolute value. The
  Identifiability section of `?fit_betadanish` explains each case.

## Documentation

* `?fit_betadanish` documents the \eqn{b = 1} ridge and the lower-tail
  \eqn{(a, c)} confounding.

* `?fit_bd_competing` documents the Tsiatis (1975) non-identifiability result
  and the direction of bias under latent dependence.

* `?BetaDanish` records that \eqn{S(t) \propto t^{-b}} in the upper tail and
  hence that \eqn{E(X^r)} is finite if and only if \eqn{b > r}.

* The `a = 1` submodel is called the Exponentiated Danish (ED) throughout, in
  line with the underlying thesis. "Complementary Exponentiated Danish (CED)"
  has been removed.

* Carbon fibre breaking stress is given in GPa. It was previously written
  "Gba", and the typo was whitelisted in `inst/WORDLIST`, which is why the
  spell check never caught it.

* The `melanoma` help page reported six columns and documented seven.

## Tests

* Regression tests pin the tail behaviour of the distribution functions. The
  survival function is regularly varying with index \eqn{-b}, so the slope of
  \eqn{\log S} against \eqn{\log t} must approach \eqn{-b}. The previous
  implementation cannot pass this test, which is why the tail defects went
  unnoticed.

* Smoke tests added for `fit_bd_aft()`, `fit_bd_cure()`, `plot.bd_aft()` and
  `plot.bd_cure()`, none of which were previously covered.

# BetaDanish 0.2.0

> **Changelog correction.** The 0.2.0 entry below has been rewritten to
> describe only what that release actually contained. As first published it
> also listed mean residual life, hazard-shape classification, stress-strength
> reliability, bootstrap confidence intervals for AFT and cure models, and a
> finite-sample simulation-study runner, none of which were implemented, and it
> reported four bug fixes that had not been applied. Those features are being
> added in the 0.3.0 development series and the fixes are recorded above.

## Major new functionality

* **Bayesian inference**: `bayes_betadanish()` provides random-walk Metropolis
  sampling for the Exponentiated Danish submodel and the full four-parameter
  Beta-Danish model with vague Gamma priors.
* **Competing risks rewrite**: `fit_bd_competing()` uses bound-constrained
  multi-start L-BFGS-B optimisation. `cif_compare()` overlays fitted cumulative
  incidence functions against the Aalen-Johansen estimator and reports Gray's
  test.
* **Structural properties**: Shannon entropy (`bd_entropy_shannon()`, by
  adaptive quadrature) and order-statistic densities (`bd_order_stat_pdf()`).
* **Diagnostics**: Cox-Snell residual plot methods for AFT and cure fits. These
  were shipped but non-functional; see the fix above.

## Vignettes

Three new vignettes were added:

* "Bayesian Estimation with BetaDanish"
* "Competing Risks with the Beta-Danish Distribution"
* "Cure Models with the Beta-Danish Distribution"

## Infrastructure

* Continuous integration via GitHub Actions on four OS/R configurations.
* Test coverage reporting via Codecov.
* Online package website built with pkgdown.
* All `Suggests` packages guarded with `requireNamespace()` at the call sites.

# BetaDanish 0.1.0

* First public release.
* Implements the four-parameter Beta-Danish distribution and its
  three-parameter Exponentiated Danish submodel for survival and reliability
  analysis.
* Maximum-likelihood estimation, goodness-of-fit, model comparison, and
  visualization.
* Built-in datasets: remission, carbon_fibres, transplant, aarset, leukemia,
  melanoma, brain_cancer.
