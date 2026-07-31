## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      warning = FALSE, message = FALSE)
library(BetaDanish)

## -----------------------------------------------------------------------------
tmp <- tempfile(fileext = ".csv")
bd_csv_template(tmp, type = "covariate", n = 6)
read.csv(tmp)

## -----------------------------------------------------------------------------
f <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
dat <- read_survival_data(f)
head(dat)

## -----------------------------------------------------------------------------
dat_all <- read_survival_data(f, covar_cols = "all", quiet = TRUE)
names(dat_all)

## -----------------------------------------------------------------------------
str(attr(dat_all, "bd_data_report"))

## -----------------------------------------------------------------------------
res <- bd_analyze_csv(f, analysis = "univariate", model = "ED",
                      compare = FALSE, n_starts = 5, seed = 1, quiet = TRUE)
res

## -----------------------------------------------------------------------------
res$tables$estimates
res$tables$goodness_of_fit

## -----------------------------------------------------------------------------
both <- bd_analyze_csv(f, model = "both", compare = FALSE,
                       n_starts = 5, seed = 2, quiet = TRUE)
both$tables$information_criteria
both$tables$likelihood_ratio_test

## -----------------------------------------------------------------------------
out <- file.path(tempdir(), "bd_results")
saved <- bd_analyze_csv(f, model = "ED", compare = FALSE,
                        output_dir = out, n_starts = 5, seed = 3, quiet = TRUE)
basename(saved$files)

## ----eval = FALSE-------------------------------------------------------------
# g <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
# 
# # Accelerated failure time
# aft <- bd_analyze_csv(g, analysis = "aft",
#                       covariates = c("age", "thickness"))
# 
# # Mixture cure model, cure fraction depending on ulceration
# cure <- bd_analyze_csv(g, analysis = "cure", cure_formula = ~ ulcer)
# 
# # Competing risks
# h  <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
# cr <- bd_analyze_csv(h, analysis = "competing", cause_col = "cause")

## -----------------------------------------------------------------------------
res$failures

