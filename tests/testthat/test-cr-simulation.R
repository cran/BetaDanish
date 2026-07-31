## Competing risks with covariates, and the simulation runners.

test_that("the design matrix builder handles all three input forms", {
  d <- data.frame(a = rnorm(10), g = factor(rep(c("x", "y"), 5)))
  f <- BetaDanish:::.bd_cr_design(~ a + g, d, n = 10)
  expect_equal(nrow(f), 10L)
  expect_false("(Intercept)" %in% colnames(f))

  m <- BetaDanish:::.bd_cr_design(matrix(1:20, 10, 2), NULL, n = 10)
  expect_equal(dim(m), c(10L, 2L))
  expect_equal(colnames(m), c("x1", "x2"))

  expect_null(BetaDanish:::.bd_cr_design(NULL, NULL, n = 10))
})

test_that("the design matrix builder rejects bad input", {
  d <- data.frame(a = rnorm(5))
  expect_error(BetaDanish:::.bd_cr_design(~ a, NULL, n = 5), "'data' is required")
  expect_error(BetaDanish:::.bd_cr_design(~ a, d, n = 9), "row")
  expect_error(BetaDanish:::.bd_cr_design(y ~ a, d, n = 5), "one-sided")
  d2 <- data.frame(a = c(1, NA, 3, 4, 5))
  expect_error(BetaDanish:::.bd_cr_design(~ a, d2, n = 5), "missing values")
})

test_that("simulated competing risks data has the right shape", {
  d <- simulate_bd_competing_data(300, seed = 5)
  expect_named(d, c("time", "cause"))
  expect_true(all(d$time > 0))
  expect_true(all(d$cause %in% 0:2))
  expect_gt(sum(d$cause == 1), 10)
  expect_gt(sum(d$cause == 2), 10)
  expect_gt(sum(d$cause == 0), 0)

  dx <- simulate_bd_competing_data(200, gammas = c(0.5, -0.3), seed = 6)
  expect_true("x" %in% names(dx))
  expect_true(all(dx$x %in% c(0, 1)))
})

test_that("competing risks fits without covariates", {
  skip_on_cran()
  d <- simulate_bd_competing_data(300, seed = 7)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 2))

  expect_s3_class(fit, "bd_competing")
  expect_equal(dim(fit$coefficients), c(2L, 3L))
  expect_equal(colnames(fit$coefficients), c("b", "c", "k"))
  expect_true(all(fit$coefficients > 0))
  expect_true(is.finite(fit$logLik))
  expect_equal(fit$nobs, 300L)
})

test_that("competing risks fits with covariates and recovers their sign", {
  skip_on_cran()
  d <- simulate_bd_competing_data(600, gammas = c(0.8, -0.8),
                                  censor_rate = 0.1, seed = 8)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, covariates = ~ x, data = d,
                     submodel = TRUE, n_starts = 3))

  expect_true("gamma_x" %in% colnames(fit$coefficients))
  expect_equal(fit$covariates, "x")
  ## Cause 1 was accelerated and cause 2 retarded; the signs should differ.
  g <- fit$coefficients[, "gamma_x"]
  expect_gt(g[1], g[2])
})

test_that("competing risks validates its inputs", {
  expect_error(fit_bd_competing(1:5, c(1, 2)), "same length")
  expect_error(fit_bd_competing(c(1, 2, 3), c(1, 1, 1)), "at least two")
  expect_error(fit_bd_competing(c(-1, 2, 3), c(0, 1, 2)), "strictly positive")
})

test_that("the CIF is monotone, bounded, and covariate-aware", {
  skip_on_cran()
  d <- simulate_bd_competing_data(400, seed = 9)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 2))

  tv <- c(1, 5, 20, 100)
  f1 <- cif_betadanish(fit, tv, cause_idx = 1)
  f2 <- cif_betadanish(fit, tv, cause_idx = 2)

  expect_true(all(diff(f1) >= -1e-8))
  expect_true(all(f1 >= 0 & f1 <= 1))
  expect_equal(cif_betadanish(fit, 0, 1), 0)
  ## The competing CIFs cannot together exceed one.
  expect_true(all(f1 + f2 <= 1 + 1e-6))

  expect_error(cif_betadanish(fit, tv, cause_idx = 9), "not among")
})

test_that("a covariate shifts the CIF", {
  skip_on_cran()
  d <- simulate_bd_competing_data(500, gammas = c(1, 0), censor_rate = 0.1,
                                  seed = 10)
  fit <- suppressWarnings(
    fit_bd_competing(d$time, d$cause, covariates = ~ x, data = d,
                     submodel = TRUE, n_starts = 3))
  a <- cif_betadanish(fit, 10, cause_idx = 1, x = 0)
  b <- cif_betadanish(fit, 10, cause_idx = 1, x = 1)
  expect_true(is.finite(a) && is.finite(b))
  expect_false(isTRUE(all.equal(a, b)))
  expect_error(cif_betadanish(fit, 10, 1, x = c(0, 0)), "one value per covariate")
})

test_that("the simulation summary helper computes what it claims", {
  set.seed(1)
  est <- rnorm(500, mean = 2.1, sd = 0.4)
  se  <- rep(0.4, 500)
  s <- BetaDanish:::.bd_sim_summary(est, se, truth = 2)
  expect_equal(unname(s["truth"]), 2)
  expect_equal(unname(s["bias"]), mean(est) - 2, tolerance = 1e-12)
  expect_equal(unname(s["rmse"]), sqrt(mean((est - 2)^2)), tolerance = 1e-12)
  expect_gt(s[["coverage"]], 0.8)
  expect_lt(abs(s[["se_ratio"]] - 1), 0.2)
})

test_that("the univariate simulation study returns a tidy table", {
  skip_on_cran()
  s <- bd_simulation_study(n = 80, n_sim = 12,
                           truth = c(b = 3, c = 2, k = 0.5),
                           submodel = TRUE, n_starts = 2,
                           seed = 11, quiet = TRUE)
  expect_s3_class(s, "bd_simulation")
  expect_true(is.data.frame(s$results))
  expect_setequal(s$results$parameter, c("b", "c", "k"))
  expect_true(all(c("bias", "rmse", "se_ratio", "coverage", "n_fail") %in%
                    names(s$results)))
  expect_true(all(s$results$truth == c(3, 2, 0.5)))
  expect_output(print(s), "simulation study")
  expect_identical(summary(s), s$results)
})

test_that("the study rejects a truth vector missing a parameter", {
  expect_error(bd_simulation_study(n = 50, n_sim = 2, truth = c(b = 3, c = 2),
                                   submodel = TRUE, quiet = TRUE),
               "must name")
})

test_that("the cure data generator censors every cured subject", {
  set.seed(2)
  d <- BetaDanish:::.bd_sim_cure_data(2000, b = 2, c = 1.5, k = 0.5,
                                      cure_fraction = 0.4, censor_rate = 0.2)
  expect_named(d, c("time", "status"))
  expect_true(all(d$time > 0))
  expect_true(all(d$status %in% c(0L, 1L)))
  ## At least the cure fraction must end up censored.
  expect_gt(mean(d$status == 0), 0.35)
})

test_that("the competing-risks simulation runner returns per-cause rows", {
  skip_on_cran()
  s <- bd_simulation_competing(n = 250, n_sim = 4, n_starts = 2,
                               seed = 12, quiet = TRUE)
  expect_s3_class(s, "bd_simulation")
  expect_true(all(grepl("^c[12]:", s$results$parameter)))
  expect_equal(nrow(s$results), 6L)
  expect_output(print(s), "competing risks")
})
