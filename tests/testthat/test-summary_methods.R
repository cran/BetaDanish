test_that("S3 methods return correct structures", {
  set.seed(123)
  dat <- data.frame(time = rbetadanish(30, 1.2, 1.5, 2, 0.5), status = 1)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, n_starts = 1, check_identifiability = FALSE)

  # Summary
  sum_fit <- summary(fit)
  expect_s3_class(sum_fit, "summary.betadanish")
  expect_true("Pr(>|z|)" %in% colnames(sum_fit$coefficients))

  # logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 4)

  # vcov and coef
  expect_equal(dim(vcov(fit)), c(4, 4))
  expect_equal(length(coef(fit)), 4)
})
