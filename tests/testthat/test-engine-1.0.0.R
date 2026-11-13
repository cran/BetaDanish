## Tests for the estimation, goodness-of-fit, interval and Bayesian routines
## brought into the package from the final thesis code in version 1.0.0.

test_that("censored goodness of fit equals the classical statistics without censoring", {
  set.seed(4242)
  x  <- stats::rexp(60)
  Ff <- function(t) stats::pexp(t, 1.3)
  us <- sort(Ff(x)); n <- length(us); i <- seq_len(n)
  classical <- c(D  = max(pmax(i / n - us, us - (i - 1) / n)),
                 W2 = 1 / (12 * n) + sum((us - (2 * i - 1) / (2 * n))^2),
                 A2 = -n - sum((2 * i - 1) * (log(us) + log(1 - rev(us)))) / n)
  got <- BetaDanish:::.bd_gof_pit(x, rep(1L, n), Ff)
  expect_equal(got, classical, tolerance = 1e-6)
})

test_that("censored goodness of fit is finite and ordered sensibly", {
  set.seed(7)
  t0 <- rbetadanish(200, 1, 3, 2, 0.5)
  cc <- stats::runif(200, 0, 8)
  time <- pmin(t0, cc); status <- as.integer(t0 <= cc)
  good <- BetaDanish:::.bd_gof_pit(time, status, function(t) ped(t, 3, 2, 0.5))
  bad  <- BetaDanish:::.bd_gof_pit(time, status, function(t) ped(t, 3, 2, 5))
  expect_true(all(is.finite(good)))
  expect_true(all(bad > good))
})

test_that("the analytic score matches a numerical gradient", {
  x <- c(0.05, 0.2, 0.5, 0.9, 1.4, 2.2, 3.7, 6.5, 12, 40)
  num_grad <- function(f, eta, h = 1e-5)
    vapply(seq_along(eta), function(j) {
      e1 <- eta; e2 <- eta; e1[j] <- e1[j] + h; e2[j] <- e2[j] - h
      (f(e1) - f(e2)) / (2 * h)
    }, numeric(1))
  ll4 <- function(eta) sum(dbetadanish(x, exp(eta[1]), exp(eta[2]),
                                       exp(eta[3]), exp(eta[4]), log = TRUE))
  ll3 <- function(eta) sum(ded(x, exp(eta[1]), exp(eta[2]), exp(eta[3]),
                               log = TRUE))
  for (th in list(c(1.5, 2.5, 2, 1), c(0.3, 6, 1.2, 0.4), c(4, 0.7, 0.8, 2.5))) {
    g <- BetaDanish:::.bd_score_eta(x, th, submodel = FALSE)
    expect_equal(g, num_grad(ll4, log(th)), tolerance = 1e-5)
  }
  for (th in list(c(2, 1.5, 0.8), c(0.6, 3, 0.3), c(7, 0.8, 2))) {
    g <- BetaDanish:::.bd_score_eta(x, th, submodel = TRUE)
    expect_equal(g, num_grad(ll3, log(th)), tolerance = 1e-5)
  }
})

test_that("pbetadanish stays finite where G underflows", {
  a <- 1.5; b <- 3; cc <- 2; k <- 1
  q_at <- function(logG) 1 / expm1(-logG / cc) / k   # time with that log G
  ## Below the switch the leading term is used; it must match the closed form.
  lg <- -701
  lF <- pbetadanish(q_at(lg), a, b, cc, k, log.p = TRUE)
  expect_true(is.finite(lF))
  expect_equal(lF, a * lg - log(a) - lbeta(a, b), tolerance = 1e-10)
  ## Just above the switch pbeta is used; the two branches agree.
  lg2 <- -699
  expect_equal(pbetadanish(q_at(lg2), a, b, cc, k, log.p = TRUE),
               a * lg2 - log(a) - lbeta(a, b), tolerance = 1e-8)
  ## The upper tail is one to double precision and its log is -F.
  expect_equal(pbetadanish(q_at(lg), a, b, cc, k, lower.tail = FALSE), 1)
  expect_equal(pbetadanish(q_at(lg), a, b, cc, k, lower.tail = FALSE, log.p = TRUE),
               -exp(lF), tolerance = 1e-10)
})

test_that("pbetadanish returns NaN, not an error, for an undefined log G", {
  ## c and k both infinite make log G = -Inf * 0; an optimiser can propose this.
  v <- suppressWarnings(pbetadanish(c(1, 2), a = 1, b = 2, c = Inf, k = Inf))
  expect_length(v, 2L)
  expect_true(all(is.nan(v)))
  v2 <- suppressWarnings(pbetadanish(c(1, 2), 1, 2, Inf, Inf,
                                     lower.tail = FALSE, log.p = TRUE))
  expect_true(all(is.nan(v2)))
})

test_that("competing risks seeds its joint fit from admissible cause-wise starts", {
  skip_on_cran()
  skip_if_not_installed("cmprsk")
  set.seed(2026)
  n  <- 400
  T1 <- rbetadanish(n, a = 1.2, b = 1.5, c = 1.0, k = 0.4)
  T2 <- rbetadanish(n, a = 1.0, b = 2.0, c = 1.0, k = 0.2)
  C  <- stats::rexp(n, 0.05)
  time  <- pmin(T1, T2, C)
  cause <- ifelse(time == C, 0L, ifelse(T1 <= T2, 1L, 2L))
  fit <- suppressWarnings(fit_bd_competing(time = time, cause = cause))
  expect_s3_class(fit, "bd_competing")
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("qbetadanish rejects malformed tail flags", {
  expect_error(qbetadanish(0.5, 1, 2, 2, 1, lower.tail = c(TRUE, FALSE)),
               "single TRUE or FALSE")
  expect_error(qbetadanish(0.5, 1, 2, 2, 1, log.p = NA), "single TRUE or FALSE")
})

test_that("a degenerate standard error is NA, never zero", {
  V <- diag(c(4, 0, -1, NA, 1e-30))
  expect_equal(BetaDanish:::.bd_se_from_vcov(V, letters[1:5]),
               c(a = 2, b = NA, c = NA, d = NA, e = NA))
})

test_that("the four-parameter fit is never below its Stage-1 ED fit", {
  skip_on_cran()
  set.seed(31)
  dat <- data.frame(time = rbetadanish(250, 1, 2.2, 1.6, 0.7), status = 1L)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   n_starts = 2, check_identifiability = FALSE))
  expect_true(is.finite(fit$submodel_logLik))
  expect_gte(fit$logLik, fit$submodel_logLik - 1e-6)
  expect_true(is.logical(fit$at_boundary))
  expect_true(is.character(fit$flag))
  expect_true(all(c("lower", "upper") %in% names(fit$bounds)))
})

test_that("a fit on the edge of the box reports no AIC, SE or goodness of fit", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 41)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, box = c(0.5, 1.5),
                   check_identifiability = FALSE))
  expect_true(fit$at_boundary)
  expect_true(nzchar(fit$edge_par))
  expect_true(is.na(fit$AIC) && is.na(fit$BIC))
  expect_true(all(is.na(fit$vcov)))
  g <- gof_betadanish(fit)
  expect_true(all(is.na(g$IC)))
  expect_true(is.na(g$KS_Statistic))
  expect_output(print(summary(fit)), "supremum")
  ## Other diagnostics may also fire on so narrow a box; only this one matters.
  w <- testthat::capture_warnings(BetaDanish:::.bd_warn_diagnostics(fit$diagnostics))
  expect_true(any(grepl("edge of the parameter box", w)))
})

test_that("gof_betadanish reports three statistics on one scale", {
  skip_on_cran()
  set.seed(12)
  x <- rbetadanish(120, 1, 3, 2, 0.5)
  dat <- data.frame(time = x, status = 1L)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 1, check_identifiability = FALSE))
  g <- gof_betadanish(fit)
  cf <- fit$coefficients
  us <- sort(ped(x, cf[["b"]], cf[["c"]], cf[["k"]]))
  n <- length(us); i <- seq_len(n)
  expect_equal(g$KS_Statistic, max(pmax(i / n - us, us - (i - 1) / n)),
               tolerance = 1e-8)
  expect_equal(g$CvM_Statistic,
               1 / (12 * n) + sum((us - (2 * i - 1) / (2 * n))^2),
               tolerance = 1e-8)
  expect_true(is.finite(g$AD_Statistic))
  expect_false(g$censored)
})

test_that("the default profile is extended until it closes", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 25)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  p <- bd_profile_ci(fit, "b", n_grid = 15L)
  expect_true(isTRUE(p$extended))
  expect_lte(p$lower, p$estimate)
  expect_gte(p$upper, p$estimate)
  if (is.finite(p$upper)) expect_true(p$upper_closed)
  expect_true(is.logical(p$lower_closed))
  expect_true(max(p$profile, na.rm = TRUE) <= p$logLik_max + 1e-6)
  expect_output(print(p), "endpoints")
})

test_that("the profile is measured from the highest log-likelihood reached", {
  skip_on_cran()
  dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 25)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3, check_identifiability = FALSE))
  ## Pretend the fit stopped 0.01 short of its maximum.
  short <- fit; short$logLik <- fit$logLik - 0.01
  expect_warning(p <- bd_profile_ci(short, "b", n_grid = 15L), "stopped short")
  expect_gte(p$logLik_max, max(p$profile, na.rm = TRUE))
  expect_equal(p$logLik_fit, short$logLik)
})

test_that("the KL separation is zero at a = 1 and positive otherwise", {
  skip_on_cran()
  k0 <- bd_kl_submodel(a = 1, b = 2, c = 2, k = 1)
  expect_lt(k0$kl, 1e-6)
  expect_equal(k0$mass, 1, tolerance = 1e-6)
  k3 <- bd_kl_submodel(a = 3, b = 2, c = 2, k = 1)
  expect_gt(k3$kl, k0$kl)
  expect_equal(k3$n_aic, 6 / k3$kl, tolerance = 1e-10)
  expect_named(k3$ed, c("b", "c", "k"))
  expect_error(bd_kl_submodel(a = -1, b = 2, c = 2, k = 1), "single positive")
})

test_that("the four-parameter cure kernel is nested over the ED kernel", {
  skip_on_cran()
  set.seed(21)
  n <- 250
  grp <- rep(0:1, length.out = n)
  cured <- stats::rbinom(n, 1, ifelse(grp == 1, 0.4, 0.2)) == 1
  t_u <- rbetadanish(n, a = 1, b = 3, c = 1.5, k = 0.3)
  cens <- stats::runif(n, 0, 30)
  dat <- data.frame(time = ifelse(cured, cens, pmin(t_u, cens)),
                    status = as.integer(!cured & t_u <= cens), group = grp)
  fit <- suppressWarnings(
    fit_bd_cure(survival::Surv(time, status) ~ 1, ~ group, data = dat,
                type = "mixture", n_starts = 2, submodel = FALSE))
  expect_s3_class(fit, "bd_cure")
  expect_true("log_a" %in% names(fit$coefficients))
  expect_false(fit$submodel)
  expect_gte(fit$logLik, fit$logLik_submodel - 1e-6)
  expect_equal(fit$npar, length(fit$coefficients))
  expect_true(is.logical(fit$separated))
  expect_output(print(summary(fit)), "LR test of a = 1")
})

test_that("the adaptive sampler returns chains, a summary and diagnostics", {
  skip_on_cran()
  dat <- simulate_bd_data(80, a = 1, b = 3, c = 2, k = 0.5, seed = 6)
  fit <- suppressWarnings(
    bayes_betadanish(dat$time, dat$status, submodel = TRUE,
                     burnin = 300, mcmc = 600, n_chains = 2, seed = 1))
  expect_s3_class(fit, "bd_bayes")
  expect_identical(fit$method, "adaptive")
  expect_length(fit$chains, 2L)
  expect_equal(nrow(fit$chains[[1]]), 600L)
  expect_true(all(c("b", "c", "k") %in% fit$summary$parameter))
  expect_true(all(c("Rhat", "ESS") %in% names(fit$diagnostics)))
  expect_true(is.logical(fit$converged))
  expect_true(all(fit$acceptance > 0 & fit$acceptance < 1))
  expect_equal(dim(fit$HPD), c(3L, 2L))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_invisible(plot(fit))
  expect_output(print(fit), "Adaptive Metropolis")
})
