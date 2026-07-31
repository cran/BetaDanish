## Upper-tail regression tests.
##
## These pin the defects fixed in 0.2.0.9000. The Beta-Danish survival function
## is regularly varying with index -b, so on a log-log scale it must approach a
## straight line of slope -b. An implementation that floors log(1 - G) or that
## forms the survival by subtracting a near-one probability from one cannot
## satisfy this, which is exactly how the earlier defects escaped notice.

test_that("the survival tail is regularly varying with index -b", {
  a <- 1.5; c <- 2; k <- 1
  for (b in c(1.5, 3, 5)) {
    t1 <- 1e8; t2 <- 1e12
    ls1 <- sbetadanish(t1, a, b, c, k, log = TRUE)
    ls2 <- sbetadanish(t2, a, b, c, k, log = TRUE)
    expect_true(is.finite(ls1) && is.finite(ls2))
    slope <- (ls2 - ls1) / (log(t2) - log(t1))
    expect_equal(slope, -b, tolerance = 1e-4)
  }
})

test_that("log-survival stays finite far beyond double-precision saturation", {
  ls <- sbetadanish(c(1e14, 1e16, 1e18), a = 1.5, b = 3, c = 2, k = 1,
                    log = TRUE)
  expect_true(all(is.finite(ls)))
  expect_true(all(diff(ls) < 0))
})

test_that("the log-density tail has index -(b+1)", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t1 <- 1e8; t2 <- 1e12
  ld1 <- dbetadanish(t1, a, b, c, k, log = TRUE)
  ld2 <- dbetadanish(t2, a, b, c, k, log = TRUE)
  expect_true(is.finite(ld1) && is.finite(ld2))
  slope <- (ld2 - ld1) / (log(t2) - log(t1))
  expect_equal(slope, -(b + 1), tolerance = 1e-4)
})

test_that("density and survival match their analytic asymptotes", {
  ## As t -> Inf, with 1 - G(t) ~ c/(kt):
  ##   log f(t) -> b log(c/k) - (b+1) log t - log B(a,b)
  ##   log S(t) -> b log(c/k) -  b    log t - log b - log B(a,b)
  ## Both were verified against a 60-digit computation. At t = 1e18 with
  ## (a,b,c,k) = (1.5,3,2,1) the density asymptote is -161.8253135; the
  ## previous implementation, which floored log(1 - G) at about log(1e-16),
  ## returned -154.0013 instead, so this test separates the two decisively.
  a <- 1.5; b <- 3; c <- 2; k <- 1; t <- 1e18

  expect_equal(dbetadanish(t, a, b, c, k, log = TRUE),
               b * log(c / k) - (b + 1) * log(t) - lbeta(a, b),
               tolerance = 1e-6)

  expect_equal(sbetadanish(t, a, b, c, k, log = TRUE),
               b * log(c / k) - b * log(t) - log(b) - lbeta(a, b),
               tolerance = 1e-6)
})

test_that("the hazard is asymptotically 1/t and does not saturate", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(1e8, 1e10, 1e12)
  h <- hbetadanish(t, a, b, c, k)
  expect_true(all(is.finite(h)))
  expect_equal(unname(h * t), rep(b, length(t)), tolerance = 1e-3)
})

test_that("quantile round-trip survives the extreme upper tail", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  p <- 1 - 10^-(4:13)
  q <- qbetadanish(p, a, b, c, k)
  expect_true(all(is.finite(q)))
  expect_true(all(diff(q) > 0))
  expect_equal(pbetadanish(q, a, b, c, k), p, tolerance = 1e-9)
})

test_that("upper-tail probabilities agree with the lower-tail complement", {
  ## Only where the complement is itself well conditioned; the point of the
  ## mirror identity is that the survival stays right after this range.
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.5, 1, 5, 20)
  expect_equal(sbetadanish(t, a, b, c, k),
               1 - pbetadanish(t, a, b, c, k), tolerance = 1e-11)
})

test_that("log.p is consistent with the probability scale", {
  a <- 1.5; b <- 3; c <- 2; k <- 1
  t <- c(0.5, 2, 10)
  expect_equal(pbetadanish(t, a, b, c, k, log.p = TRUE),
               log(pbetadanish(t, a, b, c, k)), tolerance = 1e-12)
  expect_equal(sbetadanish(t, a, b, c, k, log = TRUE),
               log(sbetadanish(t, a, b, c, k)), tolerance = 1e-12)
})
