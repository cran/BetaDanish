test_that("PDF integrates to one", {
  pars <- list(a = 1.5, b = 2.5, c = 2, k = 1)
  f <- function(t) dbetadanish(t, pars$a, pars$b, pars$c, pars$k)
  expect_equal(stats::integrate(f, 0, Inf, rel.tol = 1e-9)$value, 1,
               tolerance = 1e-6)
})

test_that("CDF matches numerical integration of the PDF", {
  a <- 1.2; b <- 3; c <- 2; k <- 0.8
  for (q in c(0.25, 1, 4)) {
    num <- stats::integrate(function(t) dbetadanish(t, a, b, c, k),
                            0, q, rel.tol = 1e-10)$value
    expect_equal(pbetadanish(q, a, b, c, k), num, tolerance = 1e-7)
  }
})

test_that("survival and CDF are complementary in the well-conditioned range", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.1, 0.5, 1, 2, 5)
  expect_equal(pbetadanish(t, a, b, c, k) +
                 sbetadanish(t, a, b, c, k),
               rep(1, length(t)), tolerance = 1e-12)
})

test_that("quantile function inverts the CDF, including far into the tail", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  p <- c(1e-10, 1e-4, 0.25, 0.5, 0.75, 1 - 1e-4, 1 - 1e-9, 1 - 1e-13)
  q <- qbetadanish(p, a, b, c, k)
  expect_true(all(is.finite(q)))
  expect_true(all(diff(q) > 0))
  expect_equal(pbetadanish(q, a, b, c, k), p, tolerance = 1e-8)
})

test_that("quantile boundaries are exact", {
  expect_equal(qbetadanish(0, 1.5, 3, 2, 1), 0)
  expect_identical(qbetadanish(1, 1.5, 3, 2, 1), Inf)
})

test_that("hazard equals density over survival", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.2, 1, 3)
  expect_equal(hbetadanish(t, a, b, c, k),
               dbetadanish(t, a, b, c, k) / sbetadanish(t, a, b, c, k),
               tolerance = 1e-12)
})

test_that("parameters recycle element-wise", {
  k <- c(0.5, 1, 2)
  got <- dbetadanish(c(1, 1, 1), a = 1, b = 2, c = 1.5, k = k)
  want <- vapply(k, function(ki) dbetadanish(1, 1, 2, 1.5, ki), numeric(1))
  expect_equal(got, want)
  expect_length(dbetadanish(1, a = 1, b = 2, c = 1.5, k = k), 3L)
})

test_that("invalid parameters give NaN and a warning", {
  expect_warning(res <- dbetadanish(1:3, a = -1, b = 2, c = 1, k = 1))
  expect_true(all(is.nan(res)))
  expect_true(all(is.nan(pbetadanish(1:3, a = 0, b = 2, c = 1, k = 1))))
  expect_true(all(is.nan(qbetadanish(0.5, a = 1, b = -2, c = 1, k = 1))))
})

test_that("boundary and missing values behave", {
  expect_equal(dbetadanish(c(-1, 0), 1.5, 3, 2, 1), c(0, 0))
  expect_equal(pbetadanish(c(-1, 0), 1.5, 3, 2, 1), c(0, 0))
  expect_equal(pbetadanish(Inf, 1.5, 3, 2, 1), 1)
  expect_equal(sbetadanish(Inf, 1.5, 3, 2, 1), 0)
  expect_true(is.na(dbetadanish(NA_real_, 1.5, 3, 2, 1)))
  expect_length(dbetadanish(numeric(0), 1.5, 3, 2, 1), 0L)
})

test_that("rbetadanish returns the requested length and respects the seed", {
  set.seed(11); x1 <- rbetadanish(50, 1.5, 3, 2, 1)
  set.seed(11); x2 <- rbetadanish(50, 1.5, 3, 2, 1)
  expect_length(x1, 50L)
  expect_identical(x1, x2)
  expect_true(all(x1 > 0))
})
