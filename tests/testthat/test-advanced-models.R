## Smoke tests for the AFT and cure paths, which had no coverage at all. The
## Cox-Snell plot methods are the specific reason: they were broken from the
## start and nothing exercised them.

make_aft_data <- function(n = 200, seed = 7) {
  set.seed(seed)
  xcov  <- stats::rnorm(n)
  k_i   <- exp(-0.5 - 0.3 * xcov)
  t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k_i)
  cens  <- stats::rexp(n, rate = 0.05)
  data.frame(time   = pmin(t_sim, cens),
             status = as.integer(t_sim <= cens),
             x      = xcov,
             group  = rep(c(0, 1), length.out = n))
}

test_that("fit_bd_aft returns a usable object", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)

  expect_s3_class(fit, "bd_aft")
  expect_true(all(c("log_b", "log_c") %in% names(fit$coefficients)))
  expect_true(is.finite(fit$logLik))
  expect_output(print(fit), "AFT")
})

test_that("summary.bd_aft reports shapes on the natural scale", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)
  s   <- summary(fit)

  expect_s3_class(s, "summary.bd_aft")
  expect_true(all(c("b", "c") %in% rownames(s$tables$shape)))
  expect_true(all(s$tables$shape[, "Estimate"] > 0))
  expect_true(all(s$tables$shape[, "Lower 95%"] > 0))
  expect_output(print(s), "natural scale")
})

test_that("plot.bd_aft draws instead of warning about NA coefficients", {
  skip_on_cran()
  dat <- make_aft_data()
  fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 3)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  ## The defect was a warning about non-finite coefficients followed by an
  ## early return. Target that specifically rather than demanding total
  ## silence, which survfit could break for unrelated reasons.
  msg <- tryCatch({ plot(fit); character(0) },
                  warning = function(w) conditionMessage(w))
  expect_false(any(grepl("non-finite coefficients", msg)))
  expect_false(any(grepl("cannot be computed", msg)))
})

test_that("the natural-scale shape helper reads log_b and log_c", {
  fake <- list(coefficients = c(log_b = log(2), log_c = log(1.5),
                                `delta_(Intercept)` = -0.5))
  got <- BetaDanish:::.bd_shape_natural(fake)
  expect_equal(got$b, 2)
  expect_equal(got$c, 1.5)
  expect_equal(got$delta, -0.5)
})

test_that("fit_bd_cure runs for both formulations", {
  skip_on_cran()
  dat <- make_aft_data(n = 250, seed = 21)
  for (ty in c("mixture", "promotion")) {
    fit <- fit_bd_cure(survival::Surv(time, status) ~ 1,
                       formula_cure = ~ group, data = dat,
                       type = ty, n_starts = 3)
    expect_s3_class(fit, "bd_cure")
    expect_true(is.finite(fit$logLik))
    expect_identical(fit$type, ty)
  }
})

test_that("fit_betadanish stores information criteria and diagnostics", {
  skip_on_cran()
  dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 3)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                   submodel = TRUE, n_starts = 3))

  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))
  expect_equal(fit$nobs, 150L)
  expect_equal(fit$npar, 3L)
  expect_true(is.list(fit$diagnostics))
  expect_equal(unname(stats::AIC(fit)), fit$AIC, tolerance = 1e-8)

  rep <- report_betadanish(fit)
  expect_output(print(rep), "AIC")
})

test_that("read_survival_data selects by name and reports on the data", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(
    data.frame(t = c(1, 2, 3, NA), ev = c(1, 0, 1, 1), grp = c("a","b","a","b")),
    tmp, row.names = FALSE)

  dat <- suppressWarnings(
    read_survival_data(tmp, time_col = "t", status_col = "ev",
                       covar_cols = "grp"))

  expect_named(dat, c("time", "status", "grp"))
  expect_equal(nrow(dat), 3L)
  expect_type(dat$grp, "character")
  expect_equal(attr(dat, "bd_data_report")$rows_dropped, 1L)
})

test_that("a covariate named time is renamed rather than colliding", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(dur = c(1, 2), ev = c(1, 1), time = c(9, 9)),
                   tmp, row.names = FALSE)

  expect_warning(
    dat <- read_survival_data(tmp, time_col = "dur", status_col = "ev",
                              covar_cols = "time"),
    "collide")
  expect_named(dat, c("time", "status", "time_cov"))
  expect_equal(dat$time, c(1, 2))
})
