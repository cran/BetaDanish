
## The Patch 2c version of this test used a single seed and passed even with
## the guard absent, which made it worthless as a regression test. The observed
## +5.2e77 runaway appeared under the analyze-csv seeds, so those are swept
## here alongside the original.

test_that("the four-parameter fit never runs away on the transplant data", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  n <- nrow(transplant)

  for (s in c(11, 12, 14, 17, 99, 101, 202)) {
    set.seed(s)
    fit <- suppressWarnings(
      fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                     submodel = FALSE, n_starts = 3,
                     check_identifiability = FALSE))

    ## The bug gave logLik / n = 5.7e75. Any real fit is far below this.
    expect_lt(fit$logLik / n, 5, label = paste("seed", s, "loglik per obs"))
    expect_gt(fit$logLik, -1e4, label = paste("seed", s, "loglik floor"))

    shapes <- fit$coefficients[intersect(names(fit$coefficients),
                                         c("a", "b", "c"))]
    expect_true(all(shapes < 500), label = paste("seed", s, "shapes bounded"))
    expect_true(all(is.finite(fit$coefficients)),
                label = paste("seed", s, "finite estimates"))
  }
})

test_that("the guard is actually in force, not merely unexercised", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(11)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = FALSE, n_starts = 3,
                   check_identifiability = FALSE))

  ## These fields only exist once fit_models.R routes through the guarded
  ## multi-start, so they double as proof the patch applied.
  expect_true(is.numeric(fit$starts_ok))
  expect_gte(fit$starts_ok, 1L)
  expect_true(is.numeric(fit$starts_rejected))
  expect_true(is.numeric(fit$loglik_spread))
  expect_gte(fit$loglik_spread, 0)
  expect_true(is.numeric(fit$diagnostics$loglik_per_obs))
})

test_that("a likelihood ratio test on real data is no longer absurd", {
  skip_on_cran()
  data(transplant, package = "BetaDanish", envir = environment())
  set.seed(11)
  full <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = FALSE, n_starts = 3, check_identifiability = FALSE))
  sub <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = transplant,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))

  chisq <- 2 * (full$logLik - sub$logLik)
  expect_true(is.finite(chisq))
  expect_gt(chisq, -1)
  expect_lt(chisq, 100)
})
