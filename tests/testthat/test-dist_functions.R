test_that("PDF integrates to 1 (Unity Check)", {
  # Full 4-parameter model
  int_full <- stats::integrate(
    f = function(x) dbetadanish(x, a = 1.5, b = 2.0, c = 3.0, k = 0.5),
    lower = 0, upper = Inf
  )
  expect_equal(int_full$value, 1.0, tolerance = 1e-4)

  # 3-parameter submodel (a = 1)
  int_sub <- stats::integrate(
    f = function(x) dbetadanish(x, a = 1.0, b = 2.5, c = 1.5, k = 1.0),
    lower = 0, upper = Inf
  )
  expect_equal(int_sub$value, 1.0, tolerance = 1e-4)
})

test_that("CDF matches numerical integration of PDF", {
  x_test <- 2.0
  a <- 2.0; b <- 1.5; c <- 2.0; k <- 0.5

  cdf_analytical <- pbetadanish(x_test, a, b, c, k)
  cdf_numerical <- stats::integrate(
    f = function(t) dbetadanish(t, a, b, c, k),
    lower = 0, upper = x_test
  )$value

  expect_equal(cdf_analytical, cdf_numerical, tolerance = 1e-5)
})

test_that("Quantile function perfectly inverts the CDF", {
  probs <- c(0.1, 0.5, 0.9)
  a <- 1.2; b <- 0.8; c <- 1.5; k <- 0.2

  quantiles <- qbetadanish(probs, a, b, c, k)
  cdf_back <- pbetadanish(quantiles, a, b, c, k)

  expect_equal(cdf_back, probs, tolerance = 1e-5)
})

test_that("Hazard function equals PDF / Survival", {
  x_test <- 1.5
  a <- 1.5; b <- 2.0; c <- 1.2; k <- 0.8

  haz_analytical <- hbetadanish(x_test, a, b, c, k)
  pdf_val <- dbetadanish(x_test, a, b, c, k)
  surv_val <- pbetadanish(x_test, a, b, c, k, lower.tail = FALSE)

  expect_equal(haz_analytical, pdf_val / surv_val, tolerance = 1e-5)
})

test_that("Invalid parameters return NaN and warnings", {
  expect_warning(res <- dbetadanish(1, a = -1, b = 2, c = 3, k = 0.5))
  expect_true(is.nan(res))
})
