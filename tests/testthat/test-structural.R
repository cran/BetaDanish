P <- list(a = 1.5, b = 3, c = 2, k = 1)

test_that("closed-form Shannon entropy matches quadrature", {
  for (p in list(list(a = 1.5, b = 3, c = 2, k = 1),
                 list(a = 1,   b = 2.5, c = 3, k = 1.2),
                 list(a = 2,   b = 0.8, c = 1.5, k = 0.5))) {
    closed <- bd_entropy_shannon(p$a, p$b, p$c, p$k)
    quad   <- bd_entropy_shannon(p$a, p$b, p$c, p$k, method = "quadrature")
    expect_equal(closed, quad, tolerance = 1e-5,
                 label = sprintf("b = %.1f", p$b))
  }
})

test_that("the series tail correction matters at small b", {
  ## Terms decay like i^-(b+1). Without the analytic tail the truncated sum is
  ## badly short at small b; with it, closed and quadrature agree.
  small <- list(a = 1, b = 0.5, c = 2, k = 1)
  closed <- bd_entropy_shannon(small$a, small$b, small$c, small$k, terms = 2000L)
  quad   <- bd_entropy_shannon(small$a, small$b, small$c, small$k,
                               method = "quadrature")
  expect_equal(closed, quad, tolerance = 1e-3)

  ## And the correction is not merely cosmetic: the raw sum is far off.
  i <- seq_len(2000L)
  raw <- lbeta(small$a, small$b) - log(small$c * small$k) -
    (small$a - 1 / small$c) * (digamma(small$a) - digamma(small$a + small$b)) -
    (small$b - 1) * (digamma(small$b) - digamma(small$a + small$b)) +
    2 * sum(exp(lbeta(small$a + i / small$c, small$b) -
                  lbeta(small$a, small$b)) / i)
  expect_gt(abs(raw - quad), 0.05)
})

test_that("Renyi and Tsallis are finite and ordered sensibly", {
  r2 <- bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 2)
  r3 <- bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 3)
  sh <- bd_entropy_shannon(P$a, P$b, P$c, P$k)
  expect_true(all(is.finite(c(r2, r3, sh))))
  ## Renyi entropy is non-increasing in its order, and Shannon is the q -> 1 limit
  expect_lt(r3, r2)
  expect_lt(r2, sh)

  t2 <- bd_entropy_tsallis(P$a, P$b, P$c, P$k, order = 2)
  expect_true(is.finite(t2))
})

test_that("order one is refused for Renyi and Tsallis", {
  expect_error(bd_entropy_renyi(P$a, P$b, P$c, P$k, order = 1), "Shannon")
  expect_error(bd_entropy_tsallis(P$a, P$b, P$c, P$k, order = 1), "Shannon")
})

test_that("stress-strength is one half for identical distributions", {
  expect_equal(bd_stress_strength(P, P), 0.5, tolerance = 1e-6)
})

test_that("stress-strength respects direction", {
  strong <- list(a = 1.5, b = 3, c = 2, k = 0.5)   # larger scale => larger X
  weak   <- list(a = 1.5, b = 3, c = 2, k = 2)
  R  <- bd_stress_strength(strength = strong, stress = weak)
  Rr <- bd_stress_strength(strength = weak,   stress = strong)
  expect_gt(R, 0.5)
  expect_equal(R + Rr, 1, tolerance = 1e-6)
})

test_that("stress-strength matches direct integration", {
  x <- list(a = 1.2, b = 4, c = 2, k = 0.8)
  y <- list(a = 1.5, b = 3, c = 2, k = 1)
  direct <- stats::integrate(
    function(t) sbetadanish(t, x$a, x$b, x$c, x$k) *
      dbetadanish(t, y$a, y$b, y$c, y$k),
    0, Inf, rel.tol = 1e-10)$value
  expect_equal(bd_stress_strength(x, y), direct, tolerance = 1e-6)
})

test_that("stress-strength validates its arguments", {
  expect_error(bd_stress_strength(list(a = 1, b = 2), P), "named values")
  expect_error(bd_stress_strength(list(a = -1, b = 2, c = 1, k = 1), P),
               "strictly positive")
})

test_that("hazard shape classification returns a known label", {
  s <- bd_hazard_shape(P$a, P$b, P$c, P$k)
  expect_s3_class(s, "bd_shape")
  expect_true(s$shape %in% c("increasing", "decreasing", "bathtub",
                             "upside-down bathtub", "indeterminate"))
  expect_length(s$time, 400L)
  expect_true(all(is.finite(s$eta[is.finite(s$hazard)])))
  expect_output(print(s), "hazard shape")
})

test_that("Glaser eta agrees with a numerical derivative of the log-density", {
  t <- c(0.3, 1, 2.5)
  h <- 1e-5
  numeric_eta <- -(dbetadanish(t + h, P$a, P$b, P$c, P$k, log = TRUE) -
                   dbetadanish(t - h, P$a, P$b, P$c, P$k, log = TRUE)) / (2 * h)
  expect_equal(BetaDanish:::.bd_glaser_eta(t, P$a, P$b, P$c, P$k),
               numeric_eta, tolerance = 1e-5)
})

test_that("the classifier recognises each shape", {
  cl <- BetaDanish:::.bd_classify
  expect_equal(cl(1:10), "increasing")
  expect_equal(cl(10:1), "decreasing")
  expect_equal(cl(c(5, 3, 1, 2, 4, 6)), "bathtub")
  expect_equal(cl(c(1, 3, 6, 4, 2)), "upside-down bathtub")
})

test_that("order statistic CDF is a proper distribution function", {
  t <- c(0.2, 1, 5, 20)
  F3 <- bd_order_stat_cdf(t, i = 3, n = 5, P$a, P$b, P$c, P$k)
  expect_true(all(diff(F3) >= 0))
  expect_true(all(F3 >= 0 & F3 <= 1))

  ## The maximum is F^n and the minimum is 1 - (1 - F)^n
  Fx <- pbetadanish(t, P$a, P$b, P$c, P$k)
  expect_equal(bd_order_stat_cdf(t, 5, 5, P$a, P$b, P$c, P$k), Fx^5,
               tolerance = 1e-10)
  expect_equal(bd_order_stat_cdf(t, 1, 5, P$a, P$b, P$c, P$k), 1 - (1 - Fx)^5,
               tolerance = 1e-10)
})

test_that("order statistic moments are ordered and respect existence", {
  m <- vapply(1:5, function(i)
    bd_order_stat_moments(1, i, 5, P$a, P$b, P$c, P$k), numeric(1))
  expect_true(all(is.finite(m)))
  expect_true(all(diff(m) > 0))          # E(Z_(1)) < ... < E(Z_(n))

  ## The maximum of a sample needs b > r, exactly as the parent does.
  expect_identical(bd_order_stat_moments(4, 5, 5, a = 1.5, b = 3, c = 2, k = 1),
                   Inf)
  ## The minimum is far better behaved: b*n > r suffices.
  expect_true(is.finite(
    bd_order_stat_moments(4, 1, 5, a = 1.5, b = 3, c = 2, k = 1)))
})

test_that("order statistic index is validated", {
  expect_error(bd_order_stat_cdf(1, i = 6, n = 5, P$a, P$b, P$c, P$k), "between")
  expect_error(bd_order_stat_moments(1, 0, 5, P$a, P$b, P$c, P$k), "between")
})

test_that("the tail index is b, with the consequences recorded", {
  ti <- bd_tail_index(P$a, P$b, P$c, P$k)
  expect_equal(ti$tail_index, P$b)
  expect_equal(ti$survival_exponent, -P$b)
  expect_false(ti$mgf_exists)
  expect_equal(ti$domain_of_attraction, "Frechet")
})

test_that("the ED API matches the a = 1 parent exactly", {
  b <- 3; c <- 2; k <- 1; x <- c(0.5, 1, 4)
  expect_equal(ded(x, b, c, k), dbetadanish(x, 1, b, c, k))
  expect_equal(ped(x, b, c, k), pbetadanish(x, 1, b, c, k))
  expect_equal(sed(x, b, c, k), sbetadanish(x, 1, b, c, k))
  expect_equal(hed(x, b, c, k), hbetadanish(x, 1, b, c, k))
  p <- c(0.1, 0.5, 0.9)
  expect_equal(qed(p, b, c, k), qbetadanish(p, 1, b, c, k))
  expect_equal(ded(x, b, c, k, log = TRUE), dbetadanish(x, 1, b, c, k, log = TRUE))
  expect_equal(ped(x, b, c, k, lower.tail = FALSE),
               pbetadanish(x, 1, b, c, k, lower.tail = FALSE))
})

test_that("the ED closed form F(t) = 1 - (1 - G)^b holds", {
  b <- 3; c <- 2; k <- 1; t <- c(0.5, 1, 4, 50)
  G <- (k * t / (1 + k * t))^c
  expect_equal(ped(t, b, c, k), 1 - (1 - G)^b, tolerance = 1e-9)
})

test_that("red respects the seed and returns positive values", {
  set.seed(3); x1 <- red(20, 3, 2, 1)
  set.seed(3); x2 <- red(20, 3, 2, 1)
  expect_identical(x1, x2)
  expect_length(x1, 20L)
  expect_true(all(x1 > 0))
})
