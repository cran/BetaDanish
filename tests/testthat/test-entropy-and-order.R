test_that("bd_entropy_shannon is finite for a regular configuration", {
  H <- bd_entropy_shannon(a = 1.5, b = 2.5, c = 2, k = 1)
  expect_true(is.finite(H))
})

test_that("bd_order_stat_pdf integrates to 1 for any (r, n)", {
  for (rn in list(c(1, 5), c(3, 10), c(10, 10))) {
    r <- rn[1]; n <- rn[2]
    val <- stats::integrate(
      function(t) bd_order_stat_pdf(t, r = r, n = n,
                                    a = 1.5, b = 2.5, c = 2, k = 1),
      lower = 0, upper = Inf,
      subdivisions = 2000, rel.tol = 1e-7)$value
    expect_equal(val, 1, tolerance = 1e-3,
                 label = sprintf("(r = %d, n = %d)", r, n))
  }
})
