test_that("plot.betadanish runs without errors", {
  set.seed(123)
  dat <- data.frame(time = rbetadanish(30, 1.2, 1.5, 2, 0.5), status = 1)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, n_starts = 1, check_identifiability = FALSE)

  # Suppress plot output during testing
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_silent(plot(fit, type = "survival"))
  expect_silent(plot(fit, type = "hazard"))
  expect_silent(plot(fit, type = "all"))
})
