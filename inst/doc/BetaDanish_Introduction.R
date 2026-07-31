## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----bd_intro_setup, include = FALSE------------------------------------------
# The blocks above this point are ```r display blocks, which knitr does
# not execute, so the package has not been attached yet.
library(BetaDanish)

## -----------------------------------------------------------------------------
p <- list(a = 1.5, b = 5, c = 2, k = 1)

bd_moment_summary(p$a, p$b, p$c, p$k)
bd_entropy_shannon(p$a, p$b, p$c, p$k)
bd_stress_strength(strength = p, stress = p)   # identical laws: one half

## -----------------------------------------------------------------------------
bd_moments(1:4, a = 1.5, b = 3, c = 2, k = 1)   # the fourth is Inf
bd_tail_index(a = 1.5, b = 3, c = 2, k = 1)$moment_condition

## ----fig.width = 6, fig.height = 4--------------------------------------------
data(guinea_pig)
ttt <- bd_ttt_plot(guinea_pig$time)
attr(ttt, "shape")

## ----eval = FALSE-------------------------------------------------------------
# fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission,
#                       submodel = TRUE)
# bd_wald_ci(fit)
# bd_profile_plot(bd_profile_ci(fit, "b"))
# 
# # For the four-parameter model, report the identified composite
# bd_identified_coef(fit_betadanish(survival::Surv(time, status) ~ 1,
#                                   data = remission))

## ----eval = FALSE-------------------------------------------------------------
# bd_simulation_study(n = c(50, 100, 200), n_sim = 500,
#                     truth = c(b = 3, c = 2, k = 0.5), submodel = TRUE)

