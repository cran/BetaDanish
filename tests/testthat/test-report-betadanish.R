test_that('advanced BetaDanish functions run on example data', {
  skip_on_cran()

  data('remission', package = 'BetaDanish')

  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission, n_starts = 1, check_identifiability = FALSE)
  expect_s3_class(fit, 'betadanish')

  rep <- report_betadanish(fit)
  expect_s3_class(rep, 'betadanish_report')
  expect_true(is.list(rep))
})
