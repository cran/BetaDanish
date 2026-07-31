P <- list(a = 1.5, b = 5, c = 2, k = 1)     # b = 5 so all four moments exist

test_that("raw moments agree with direct integration of z^r f(z)", {
  for (r in 1:3) {
    direct <- stats::integrate(
      function(z) z^r * dbetadanish(z, P$a, P$b, P$c, P$k),
      0, Inf, rel.tol = 1e-10)$value
    expect_equal(bd_moments(r, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
  }
})

test_that("the zeroth moment is one", {
  expect_equal(bd_moments(0, P$a, P$b, P$c, P$k), 1, tolerance = 1e-8)
})

test_that("moments respect the b > r existence condition", {
  ## The survival tail has index -b, so E(Z^r) is finite iff b > r.
  expect_true(is.finite(bd_moments(2, a = 1.5, b = 2.5, c = 2, k = 1)))
  expect_identical(bd_moments(3, a = 1.5, b = 2.5, c = 2, k = 1), Inf)
  expect_identical(bd_moments(4, a = 1.5, b = 3,   c = 2, k = 1), Inf)
  expect_identical(bd_moments(3, a = 1.5, b = 3,   c = 2, k = 1), Inf)  # b = r
})

test_that("bd_moment_summary reports NA rather than nonsense past the boundary", {
  s <- bd_moment_summary(P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(s)))
  expect_gt(s[["variance"]], 0)
  expect_equal(s[["sd"]], sqrt(s[["variance"]]))

  s2 <- bd_moment_summary(a = 1.5, b = 2.5, c = 2, k = 1)
  expect_true(is.finite(s2[["mean"]]))
  expect_true(is.finite(s2[["variance"]]))
  expect_true(is.na(s2[["skewness"]]))
  expect_true(is.na(s2[["kurtosis"]]))
})

test_that("incomplete moments split the complete moment exactly", {
  for (t in c(0.2, 1, 5)) {
    lo <- bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = TRUE)
    up <- bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = FALSE)
    expect_equal(lo + up, bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  }
})

test_that("incomplete moment boundaries behave", {
  expect_equal(bd_incomplete_moment(1, 0, P$a, P$b, P$c, P$k), 0)
  expect_equal(bd_incomplete_moment(1, Inf, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  expect_true(is.na(bd_incomplete_moment(1, NA_real_, P$a, P$b, P$c, P$k)))
})

test_that("conditional moments match the incomplete/survival ratio", {
  t <- 1.5
  expect_equal(
    bd_conditional_moment(1, t, P$a, P$b, P$c, P$k, upper = TRUE),
    bd_incomplete_moment(1, t, P$a, P$b, P$c, P$k, lower = FALSE) /
      sbetadanish(t, P$a, P$b, P$c, P$k),
    tolerance = 1e-9)
})

test_that("MRL matches its definition and is zero-limit consistent", {
  t <- c(0.5, 1, 3)
  m <- bd_mrl(t, P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(m)))
  expect_true(all(m > 0))

  ## m(0+) = E(Z), by definition
  expect_equal(bd_mrl(1e-8, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-4)

  ## E(Z - t | Z > t) computed the long way
  tt <- 1
  direct <- stats::integrate(
    function(z) (z - tt) * dbetadanish(z, P$a, P$b, P$c, P$k),
    tt, Inf, rel.tol = 1e-10)$value / sbetadanish(tt, P$a, P$b, P$c, P$k)
  expect_equal(bd_mrl(tt, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("reversed MRL matches its definition", {
  tt <- 2
  direct <- tt - stats::integrate(
    function(z) z * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, tt, rel.tol = 1e-10)$value / pbetadanish(tt, P$a, P$b, P$c, P$k)
  expect_equal(bd_rmrl(tt, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("MRL requires a finite mean", {
  expect_warning(m <- bd_mrl(1, a = 1.5, b = 0.8, c = 2, k = 1), "b > 1")
  expect_true(is.na(m))
  expect_warning(bd_rmrl(1, a = 1.5, b = 0.8, c = 2, k = 1), "b > 1")
})

test_that("mean deviation about the mean matches E|Z - mu|", {
  mu <- bd_moments(1, P$a, P$b, P$c, P$k)
  direct <- stats::integrate(
    function(z) abs(z - mu) * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_mean_deviation(P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("mean deviation about the median matches E|Z - M|", {
  med <- qbetadanish(0.5, P$a, P$b, P$c, P$k)
  direct <- stats::integrate(
    function(z) abs(z - med) * dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_mean_deviation(P$a, P$b, P$c, P$k, about = "median"),
               direct, tolerance = 1e-6)
})

test_that("the Lorenz curve has the right endpoints and shape", {
  p <- c(0, 0.25, 0.5, 0.75, 1)
  L <- bd_lorenz(p, P$a, P$b, P$c, P$k)
  expect_equal(L[1], 0, tolerance = 1e-8)
  expect_equal(L[5], 1, tolerance = 1e-6)
  expect_true(all(diff(L) > 0))          # increasing
  expect_true(all(L[2:4] < p[2:4]))      # below the equality line
})

test_that("Bonferroni is Lorenz divided by p", {
  p <- c(0.25, 0.5, 0.75)
  expect_equal(bd_bonferroni(p, P$a, P$b, P$c, P$k),
               bd_lorenz(p, P$a, P$b, P$c, P$k) / p, tolerance = 1e-10)
})

test_that("PWM with zero weights reduces to the raw moment", {
  expect_equal(bd_pwm(1, 0, 0, P$a, P$b, P$c, P$k),
               bd_moments(1, P$a, P$b, P$c, P$k), tolerance = 1e-7)
  expect_equal(bd_pwm(0, 0, 0, P$a, P$b, P$c, P$k), 1, tolerance = 1e-8)
})

test_that("PWM matches direct integration", {
  direct <- stats::integrate(
    function(z) z * pbetadanish(z, P$a, P$b, P$c, P$k) *
      dbetadanish(z, P$a, P$b, P$c, P$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_pwm(1, 1, 0, P$a, P$b, P$c, P$k), direct, tolerance = 1e-6)
})

test_that("invalid parameters are rejected", {
  expect_error(bd_moments(1, a = -1, b = 3, c = 2, k = 1), "strictly positive")
  expect_error(bd_moments(-1, P$a, P$b, P$c, P$k), "non-negative")
})
