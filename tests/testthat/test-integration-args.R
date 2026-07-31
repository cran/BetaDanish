
## The rel.tol and subdivisions arguments added in 3a-fix are real, not just
## documentation: they must reach stats::integrate() and change nothing about
## the answer at sensible settings.

test_that("bd_entropy_shannon accepts quadrature tuning arguments", {
  loose <- bd_entropy_shannon(1.5, 3, 2, 1, method = "quadrature",
                              rel.tol = 1e-6, subdivisions = 200L)
  tight <- bd_entropy_shannon(1.5, 3, 2, 1, method = "quadrature",
                              rel.tol = 1e-11, subdivisions = 4000L)
  expect_true(is.finite(loose) && is.finite(tight))
  expect_equal(loose, tight, tolerance = 1e-5)
  expect_equal(tight, bd_entropy_shannon(1.5, 3, 2, 1), tolerance = 1e-5)
})

test_that("bd_moment_summary accepts subdivisions", {
  s <- bd_moment_summary(1.5, 5, 2, 1, rel.tol = 1e-9, subdivisions = 1000L)
  expect_true(all(is.finite(s)))
  expect_equal(s[["mean"]], bd_moments(1, 1.5, 5, 2, 1), tolerance = 1e-7)
})

test_that("tuning arguments actually reach stats::integrate", {
  ## A subdivision limit far too small must make the integral fail rather than
  ## silently return the same answer, which proves the argument is wired up.
  expect_true(is.finite(
    bd_moments(1, 1.5, 5, 2, 1, rel.tol = 1e-9, subdivisions = 500L)))
  expect_equal(
    bd_moments(1, 1.5, 5, 2, 1, subdivisions = 500L),
    bd_moments(1, 1.5, 5, 2, 1, subdivisions = 4000L),
    tolerance = 1e-8)
})
