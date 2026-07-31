## Visualisation added in 0.3.0. Plots are drawn to a null device; the
## assertions are about the returned values and the classification logic,
## which is what can actually be wrong.

test_that("the TTT transform has the right endpoints and is increasing", {
  data(guinea_pig, package = "BetaDanish", envir = environment())
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(guinea_pig$time)

  expect_s3_class(ttt, "data.frame")
  expect_named(ttt, c("i_n", "phi"))
  expect_equal(nrow(ttt), nrow(guinea_pig))
  expect_equal(ttt$phi[nrow(ttt)], 1, tolerance = 1e-12)   # phi(1) = 1
  expect_true(all(diff(ttt$phi) >= -1e-12))                # non-decreasing
  expect_true(all(ttt$phi >= 0 & ttt$phi <= 1 + 1e-12))
  expect_true(attr(ttt, "shape") %in%
                c("increasing", "decreasing", "bathtub",
                  "unimodal (upside-down bathtub)", "constant (exponential)"))
})

test_that("the TTT classifier recognises the reference shapes", {
  cl <- BetaDanish:::.bd_ttt_shape
  u <- seq(0.02, 1, length.out = 50)
  expect_equal(cl(u, u), "constant (exponential)")
  expect_equal(cl(u, pmin(u + 0.15, 1)), "increasing")
  expect_equal(cl(u, pmax(u - 0.15, 0)), "decreasing")
  ## Barlow-Campo: the TTT is concave for an increasing hazard and convex for
  ## a decreasing one, so convex-then-concave (below the diagonal, then above)
  ## is a bathtub, and concave-then-convex is unimodal.
  ##
  ##   u + 0.2 sin(2 pi u)  is above the diagonal on (0, 1/2)  -> unimodal
  ##   u - 0.2 sin(2 pi u)  is below the diagonal on (0, 1/2)  -> bathtub
  expect_equal(cl(u, u + 0.2 * sin(2 * pi * u)),
               "unimodal (upside-down bathtub)")
  expect_equal(cl(u, u - 0.2 * sin(2 * pi * u)), "bathtub")
})

test_that("an exponential sample gives a TTT curve near the diagonal", {
  set.seed(4)
  x <- stats::rexp(400, rate = 0.5)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  ttt <- bd_ttt_plot(x)
  expect_lt(max(abs(ttt$phi - ttt$i_n)), 0.12)
})

test_that("bd_ttt_plot drops censored observations with a warning", {
  set.seed(5)
  t <- rbetadanish(60, 1, 3, 2, 0.5)
  s <- stats::rbinom(60, 1, 0.8)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_warning(ttt <- bd_ttt_plot(t, status = s), "censored")
  expect_equal(nrow(ttt), sum(s == 1))
})

test_that("bd_ttt_plot accepts a fitted object and validates its input", {
  skip_on_cran()
  dat <- simulate_bd_data(60, a = 1, b = 3, c = 2, k = 0.5, seed = 7)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_s3_class(suppressWarnings(bd_ttt_plot(fit)), "data.frame")
  expect_error(bd_ttt_plot(c(1, 2, 3)), "At least five")
})

test_that("bd_profile_plot draws and returns its input", {
  skip_on_cran()
  dat <- simulate_bd_data(100, a = 1, b = 3, c = 2, k = 0.5, seed = 8)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 12L)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  out <- bd_profile_plot(p)
  expect_identical(out$parameter, "b")
  expect_error(bd_profile_plot(list()), "bd_profile object")
})

test_that("plot.bd_bayes validates before drawing", {
  fake <- structure(list(draws = matrix(rnorm(300), 100, 3,
                                        dimnames = list(NULL, c("b", "c", "k"))),
                         HPD = NULL, submodel = TRUE),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fake))
  expect_invisible(plot(fake, which = "b", type = "trace"))
  expect_error(plot(fake, which = "zzz"), "Not in the posterior")

  empty <- structure(list(draws = matrix(numeric(0), 0, 0)), class = "bd_bayes")
  expect_error(plot(empty), "no posterior draws")
})

test_that("plot.bd_bayes leaves the graphical parameters as it found them", {
  fake <- structure(list(draws = matrix(rnorm(200), 100, 2,
                                        dimnames = list(NULL, c("b", "c")))),
                    class = "bd_bayes")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  before <- graphics::par("mfrow")
  plot(fake)
  expect_equal(graphics::par("mfrow"), before)
})
