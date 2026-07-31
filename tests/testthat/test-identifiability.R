## The identifiability diagnostics are tested directly rather than through a
## fitted model, so the test does not depend on where the optimiser lands.

mk_diag <- function(...) {
  base <- list(converged = TRUE, convergence_code = 1L, vcov_singular = FALSE,
               near_b_ridge = FALSE, b_distance_se = 5, ac_correlation = 0.2,
               starts_ok = 12L, starts_rejected = 0L, loglik_spread = 0.1,
               loglik_per_obs = -2.5)
  utils::modifyList(base, list(...))
}

test_that("the b = 1 ridge warning fires when b-hat sits close to one", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(near_b_ridge = TRUE,
                                              b_distance_se = 0.6)),
    "non-identifiability ridge")
})

test_that("a singular information matrix is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(vcov_singular = TRUE)),
    "singular")
})

test_that("(a, c) confounding is reported above the 0.95 threshold", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = 0.99)),
    "only the product")
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(ac_correlation = -0.98)),
    "only the product")
})

test_that("a poor convergence code is reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(converged = FALSE,
                                              convergence_code = 4L)),
    "code 4")
})

test_that("a well-behaved fit produces no diagnostic warnings", {
  expect_silent(BetaDanish:::.bd_warn_diagnostics(mk_diag()))
})

test_that("check_identifiability = FALSE suppresses the warnings", {
  skip_on_cran()
  dat <- simulate_bd_data(120, a = 1, b = 1.1, c = 2, k = 0.5, seed = 99)
  fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                        n_starts = 3, check_identifiability = FALSE)
  expect_s3_class(fit, "betadanish")
  expect_true(is.list(fit$diagnostics))
})

test_that("zero-length input gives zero-length output for every function", {
  z <- numeric(0)
  expect_length(dbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(pbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(qbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(sbetadanish(z, 1.5, 3, 2, 1), 0L)
  expect_length(hbetadanish(z, 1.5, 3, 2, 1), 0L)
  ## A zero-length parameter also collapses the result, as in stats::dnorm.
  expect_length(dbetadanish(1, a = z, b = 3, c = 2, k = 1), 0L)
})

test_that("degenerate starts and local-optima spread are reported", {
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(starts_rejected = 4L)),
    "degenerate ridge")
  expect_warning(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(loglik_spread = 7.5)),
    "several local maxima")
  expect_silent(
    BetaDanish:::.bd_warn_diagnostics(mk_diag(starts_rejected = 0L,
                                              loglik_spread = 0.4)))
})

test_that("an incomplete diagnostics list degrades quietly", {
  ## Objects saved before a later diagnostic existed carry no starts_rejected
  ## or loglik_spread. Every field is normalised before it is tested, so an
  ## absent one is simply not reported.
  legacy <- list(converged = TRUE, convergence_code = 1L,
                 vcov_singular = FALSE, near_b_ridge = FALSE,
                 b_distance_se = 5, ac_correlation = 0.2)
  expect_silent(BetaDanish:::.bd_warn_diagnostics(legacy))

  ## Only the convergence flag survives.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(list(converged = TRUE)))

  ## Fields present but of the wrong shape or type.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(
    list(converged = TRUE, ac_correlation = NULL,
         starts_rejected = character(0), loglik_spread = "n/a")))
  expect_silent(BetaDanish:::.bd_warn_diagnostics(
    list(converged = TRUE, ac_correlation = c(0.1, 0.2))))

  ## An empty list has no convergence flag, so it warns -- but it must warn,
  ## not error.
  expect_warning(BetaDanish:::.bd_warn_diagnostics(list()), "optimiser")

  ## A non-list is ignored outright.
  expect_silent(BetaDanish:::.bd_warn_diagnostics(NULL))
})
