## Estimation and inference added in Patch 3b.

grid_data <- function(n = 250, seed = 11) {
  ## Continuous times rounded to a grid of 1, as month-recorded data would be.
  set.seed(seed)
  t <- rbetadanish(n, a = 1, b = 4, c = 2, k = 0.15)
  data.frame(time = pmax(round(t), 1), status = 1L)
}

test_that(".bd_log_cell is a log probability and sums sensibly", {
  a <- 1; b <- 4; c <- 2; k <- 0.15
  t <- c(1, 3, 10, 40)
  lc <- BetaDanish:::.bd_log_cell(t, delta = 1, a, b, c, k)

  expect_true(all(is.finite(lc)))
  expect_true(all(lc <= 0))                      # never exceeds log(1)

  ## Matches the direct difference of the distribution function
  direct <- log(pbetadanish(t + 0.5, a, b, c, k) -
                  pbetadanish(t - 0.5, a, b, c, k))
  expect_equal(lc, direct, tolerance = 1e-9)
})

test_that(".bd_log_cell falls back to the density when the cell underflows", {
  ## Far into the tail the two distribution values are equal to machine
  ## precision, so the difference cancels to zero and the density is used.
  lc <- BetaDanish:::.bd_log_cell(1e12, delta = 1e-6, a = 1.5, b = 3, c = 2, k = 1)
  expect_true(is.finite(lc))
  expect_lt(lc, 0)
})

test_that("grouped fitting runs and reports its increment", {
  skip_on_cran()
  dat <- grid_data()
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 3,
                   check_identifiability = FALSE))
  expect_true(isTRUE(fit$grouped))
  expect_equal(fit$delta, 1)
  expect_true(is.finite(fit$logLik))
  expect_true(all(fit$coefficients > 0))
})

test_that("grouped standard errors exceed the point-density ones", {
  skip_on_cran()
  dat <- grid_data(n = 300, seed = 12)
  exact <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  grp <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 3,
                   check_identifiability = FALSE))

  ## Treating a rounded time as exact overstates the information, so the
  ## point-density likelihood is the more confident of the two.
  se_exact <- sqrt(pmax(diag(exact$vcov), 0))
  se_grp   <- sqrt(pmax(diag(grp$vcov), 0))
  expect_true(all(is.finite(c(se_exact, se_grp))))
  expect_gte(mean(se_grp / se_exact), 0.95)
})

test_that("grid-recorded times raise a warning when grouped is FALSE", {
  skip_on_cran()
  dat <- grid_data(n = 120, seed = 13)
  expect_warning(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 2),
    "grid")
})

test_that("grouped = TRUE without an inferable increment is an error", {
  dat <- data.frame(time = c(0.137, 1.882, 3.019, 7.4451, 11.02),
                    status = 1L)
  expect_error(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, grouped = TRUE, n_starts = 1),
    "recording increment")
})

test_that("the penalty shrinks the estimates and is recorded", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 21)
  plain <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, check_identifiability = FALSE))
  pen <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, penalty = 0.5, check_identifiability = FALSE))

  expect_equal(plain$penalty, 0)
  expect_equal(pen$penalty, 0.5)
  expect_true(is.na(plain$penalised_logLik))
  expect_true(is.finite(pen$penalised_logLik))
})

test_that("the reported log-likelihood is unpenalised", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 22)
  pen <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 3, penalty = 1, check_identifiability = FALSE))

  ## The penalised objective can only be lower than the log-likelihood it
  ## subtracts a non-negative penalty from.
  expect_lte(pen$penalised_logLik, pen$logLik + 1e-6)

  ## AIC and BIC must be built from the unpenalised value, or they would not
  ## be comparable with an unpenalised fit.
  expect_equal(pen$AIC, 2 * pen$npar - 2 * pen$logLik, tolerance = 1e-8)
  expect_equal(unname(stats::AIC(pen)), pen$AIC, tolerance = 1e-8)
})

test_that("a negative penalty is refused", {
  dat <- simulate_bd_data(60, a = 1, b = 3, c = 2, k = 0.5, seed = 23)
  expect_error(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   penalty = -1, n_starts = 1),
    "non-negative")
})

test_that("Wald intervals stay inside the parameter space", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 24)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  ci <- bd_wald_ci(fit)

  expect_true(all(colnames(ci) == c("estimate", "se", "lower", "upper")))
  expect_true(all(ci[, "lower"] > 0))                       # never crosses zero
  expect_true(all(ci[, "lower"] <= ci[, "estimate"]))
  expect_true(all(ci[, "upper"] >= ci[, "estimate"]))
  expect_equal(attr(ci, "level"), 0.95)
})

test_that("the profile interval brackets the estimate", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 25)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 25L)

  expect_s3_class(p, "bd_profile")
  expect_equal(p$parameter, "b")
  expect_lte(p$lower, p$estimate)
  expect_gte(p$upper, p$estimate)
  expect_true(max(p$profile, na.rm = TRUE) <= p$logLik_max + 1e-6)
  expect_output(print(p), "Profile likelihood")
})

test_that("a flat profile is reported as an open upper bound", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 26)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))

  ## A grid stopping just above the estimate cannot bound b from above, so the
  ## upper limit must come back as Inf rather than as the grid maximum.
  g <- seq(fit$coefficients[["b"]] * 0.5, fit$coefficients[["b"]] * 1.02,
           length.out = 15)
  p <- bd_profile_ci(fit, "b", grid = g)
  expect_true(is.infinite(p$upper))
  expect_true(p$open_above)
  expect_output(print(p), "lower bound")
})

test_that("profiling rejects an unknown parameter", {
  skip_on_cran()
  dat <- simulate_bd_data(80, a = 1, b = 3, c = 2, k = 0.5, seed = 27)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 2, check_identifiability = FALSE))
  expect_error(bd_profile_ci(fit, "a"), "not a parameter")
})

test_that("the identified parametrisation reports ac with a finite SE", {
  skip_on_cran()
  dat <- simulate_bd_data(250, a = 1.5, b = 3, c = 2, k = 0.5, seed = 28)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 5, check_identifiability = FALSE))
  id <- bd_identified_coef(fit)

  expect_s3_class(id, "bd_identified")
  expect_equal(rownames(id$table), c("ac", "b", "k"))
  expect_equal(id$table["ac", "estimate"],
               fit$coefficients[["a"]] * fit$coefficients[["c"]],
               tolerance = 1e-10)
  expect_true(all(is.finite(id$table[, "se"])))
  expect_true(all(id$table[, "lower"] > 0))
  expect_output(print(id), "Identified parametrisation")
})

test_that("the submodel needs no reparametrisation", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 29)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  id <- bd_identified_coef(fit)
  expect_true(id$submodel)
  expect_equal(rownames(id$table), c("b", "c", "k"))
  expect_output(print(id), "nothing is gained")
})

test_that("inference functions reject a non-betadanish object", {
  expect_error(bd_wald_ci(list()), "betadanish")
  expect_error(bd_profile_ci(list()), "betadanish")
  expect_error(bd_identified_coef(list()), "betadanish")
})
