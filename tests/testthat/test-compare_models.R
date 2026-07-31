test_that("gof_betadanish computes IC and KS correctly", {
  set.seed(123)
  dat <- data.frame(time = rbetadanish(30, 1.2, 1.5, 2, 0.5), status = 1)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, n_starts = 1, check_identifiability = FALSE)

  gof <- gof_betadanish(fit)
  expect_type(gof, "list")
  expect_true(all(c("AIC", "BIC", "HQIC", "AICC") %in% names(gof$IC)))
  expect_true(is.numeric(gof$KS_Statistic))
})

test_that("compare_models performs LRT correctly", {
  set.seed(123)
  dat <- data.frame(time = rbetadanish(40, 1.0, 1.5, 2, 0.5), status = 1)

  fit_full <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, submodel = FALSE, n_starts = 1, check_identifiability = FALSE)
  fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, submodel = TRUE, n_starts = 1, check_identifiability = FALSE)

  res <- compare_models(fit_full, fit_sub)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true("Pr(>Chisq)" %in% colnames(res))
})
