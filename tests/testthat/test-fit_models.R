test_that("fit_betadanish works for complete data", {
  set.seed(42)
  sim_time <- rbetadanish(50, a = 1.5, b = 2.0, c = 1.5, k = 0.5)
  sim_status <- rep(1, 50)
  dat <- data.frame(time = sim_time, status = sim_status)

  # Fit full model
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, n_starts = 2, check_identifiability = FALSE)

  expect_s3_class(fit, "betadanish")
  expect_equal(length(fit$coefficients), 4)
  expect_true(all(fit$coefficients > 0))
  expect_true(is.numeric(fit$logLik))
})

test_that("fit_betadanish works for the 3-parameter submodel", {
  set.seed(42)
  sim_time <- rbetadanish(50, a = 1.0, b = 2.0, c = 1.5, k = 0.5)
  dat <- data.frame(time = sim_time, status = 1)

  fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, submodel = TRUE, n_starts = 2, check_identifiability = FALSE)

  expect_s3_class(fit_sub, "betadanish")
  expect_equal(length(fit_sub$coefficients), 3)
  expect_true(fit_sub$submodel)
})
