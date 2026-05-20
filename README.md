# BetaDanish: The Beta-Danish Distribution for Lifetime Data Analysis

<!-- badges: start -->
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.r-project.org/Licenses/GPL-3)
<!-- badges: end -->

The **BetaDanish** R package provides a comprehensive suite of tools for survival and reliability analysis using the highly flexible four-parameter Beta-Danish distribution and its three-parameter submodel.

Developed as part of doctoral research at the Department of Statistics, Allama Iqbal Open University (AIOU), Islamabad, this package addresses the limitations of classical distributions (like the Weibull or Gamma) by accommodating decreasing, increasing, unimodal, and bathtub-shaped hazard rates.


## Why use BetaDanish?

The Beta-Danish distribution is useful when classical lifetime models are too restrictive. 
It provides additional flexibility for modeling complex hazard-rate shapes commonly observed in survival, reliability, and biomedical data.

| Model | Typical hazard shape | Flexibility |
|---|---|---|
| Exponential | Constant | Low |
| Weibull | Increasing or decreasing | Moderate |
| Gamma | Flexible but limited | Moderate |
| Log-normal | Non-monotone | Moderate |
| Log-logistic | Non-monotone | Moderate |
| Beta-Danish | Increasing, decreasing, unimodal, bathtub-shaped | High |

## Features

* **Core Mathematics:** Numerically stable implementations of the density (`d`), distribution (`p`), quantile (`q`), random generation (`r`), and hazard (`h`) functions.
* **Robust Estimation:** Maximum Likelihood Estimation (MLE) engine supporting both complete and right-censored data via `survival::Surv` objects.
* **Model Comparison:** Automated benchmarking against standard lifetime distributions (Weibull, Log-Normal, Log-Logistic, Gamma, Exponential) with ranked Goodness-of-Fit tables.
* **Advanced Modules:** 
  * Accelerated Failure Time (AFT) regression.
  * Mixture and Promotion-Time (Non-Mixture) Cure Models.
  * Competing Risks analysis assuming independent latent failure times.
* **Publication-Ready Plots:** Built-in S3 methods for plotting survival, hazard, density, P-P, and Q-Q plots.

## Installation

You can install the development version of BetaDanish from GitHub using the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("bilal-aiou/BetaDanish")
```

## Quick Start

### 1. Basic Distribution Functions

```r
library(BetaDanish)

# Generate 100 random survival times
set.seed(2026)
sim_data <- rbetadanish(n = 100, a = 1.5, b = 2.0, c = 3.0, k = 0.5)

# Calculate the hazard rate at time t = 2
hbetadanish(x = 2, a = 1.5, b = 2.0, c = 3.0, k = 0.5)
```

### 2. Fitting a Model to Data

The package includes several built-in datasets, such as `remission` (bladder cancer remission times).

```r
# Load built-in dataset
data("remission", package = "BetaDanish")

# Fit the 4-parameter Beta-Danish model
fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission)

# View the summary (Estimates, Standard Errors, p-values, AIC/BIC)
summary(fit)

# Generate all diagnostic plots (Survival, Hazard, Density, PP, QQ)
plot(fit, type = "all")
```

### 3. Comparing Models

You can easily compare the 4-parameter full model against the 3-parameter submodel (where `a = 1`), or benchmark it against standard distributions.

```r
# Fit the 3-parameter submodel
fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission, submodel = TRUE)

# Likelihood Ratio Test
compare_models(full_model = fit, sub_model = fit_sub)

# Benchmark against Weibull, Gamma, Log-Normal, etc. (requires 'flexsurv')
compare_distributions(fit)
```

### 4. Advanced: Cure Models

For datasets with long-term survivors (e.g., the built-in `transplant` dataset), you can fit mixture or promotion-time cure models.

```r
data("transplant", package = "BetaDanish")

# Fit a mixture cure model (latency ~ 1, cure fraction ~ group)
cure_fit <- fit_bd_cure(
  formula_aft = survival::Surv(time, status) ~ 1, 
  formula_cure = ~ group, 
  data = transplant, 
  type = "mixture"
)

summary(cure_fit)
```

## Automated Pipeline

For a quick, end-to-end analysis of your own CSV data, use the automated reporting function:

```r
# analyze_betadanish("path/to/your_data.csv", time_col = "time", status_col = "status")
```


## Limitations

BetaDanish currently focuses primarily on complete and right-censored survival data. 
Users should check convergence, compare alternative models, and inspect diagnostic plots before drawing final conclusions. 
Covariate modeling is available through dedicated advanced functions such as `fit_bd_aft()`.

## License

This package is released under the GPL-3 License.
