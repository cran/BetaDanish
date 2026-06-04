## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  warning = FALSE,
  message = FALSE
)

## ----example------------------------------------------------------------------
library(BetaDanish)
set.seed(2026)
n <- 250
z <- stats::rbinom(n, 1, 0.5)
pi_susc <- stats::plogis(0.3 + 0.7 * z)
cured   <- stats::rbinom(n, 1, 1 - pi_susc) == 1
T_true  <- ifelse(cured, Inf,
                  rbetadanish(n, a = 1, b = 2, c = 1.5, k = 0.4))
C   <- stats::rexp(n, 0.04)
time   <- pmin(T_true, C)
status <- ifelse(T_true <= C, 1, 0)
dat <- data.frame(time = time, status = status, z = z)
cat("Sample size:", n, " Censoring rate:", round(mean(status == 0), 2), "\n")

## ----mix, error = TRUE--------------------------------------------------------
try({
fit_mix <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ z,
  data         = dat,
  type         = "mixture",
  n_starts     = 3
)
summary(fit_mix)
})

## ----prom, error = TRUE-------------------------------------------------------
try({
fit_prom <- fit_bd_cure(
  formula_aft  = survival::Surv(time, status) ~ 1,
  formula_cure = ~ z,
  data         = dat,
  type         = "promotion",
  n_starts     = 3
)
summary(fit_prom)
})

