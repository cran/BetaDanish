## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## -----------------------------------------------------------------------------
library(BetaDanish)
library(survival)
data('remission', package = 'BetaDanish')
head(remission)

## -----------------------------------------------------------------------------
fit <- fit_betadanish(Surv(time, status) ~ 1, data = remission, n_starts = 1)
summary(fit)

## -----------------------------------------------------------------------------
fit_sub <- fit_betadanish(Surv(time, status) ~ 1, data = remission, submodel = TRUE, n_starts = 1)
compare_models(fit, fit_sub)

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot(fit, type = 'survival')
plot(fit, type = 'hazard')

