## Regression tests for the degenerate-ridge guard.
##
## The motivating failure: a four-parameter fit on the transplant data returned
## a log-likelihood of +5.2e77 on 91 observations, and the unguarded
## multi-start selected it as the maximum.

fake_fit <- function(ll, ...) {
  list(maximum = ll, estimate = c(...))
}

test_that("the acceptance predicate rejects exploded shapes", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)

  ok <- fake_fit(-250, log_a = log(1.5), log_b = log(3),
                 log_c = log(2), log_k = log(0.5))
  expect_true(acc(ok))

  ## The master script's recorded case: c = 296209.74
  blown <- fake_fit(-934.68, log_a = log(1), log_b = log(0.82),
                    log_c = log(296209.74), log_k = log(0.5))
  expect_false(acc(blown))

  ## Exactly at the limit is admissible; past it is not.
  expect_true(acc(fake_fit(-250, log_a = log(1), log_b = log(499),
                           log_c = log(2), log_k = log(0.5))))
  expect_false(acc(fake_fit(-250, log_a = log(1), log_b = log(501),
                            log_c = log(2), log_k = log(0.5))))
})

test_that("the acceptance predicate rejects an implausible log-likelihood", {
  acc <- BetaDanish:::.bd_make_accept(n = 91)
  ## The observed runaway. Shapes alone would not have caught it.
  expect_false(acc(fake_fit(5.210644e77, log_a = log(2), log_b = log(3),
                            log_c = log(2), log_k = log(0.1))))
  expect_false(acc(fake_fit(Inf, log_a = 0, log_b = 0, log_c = 0, log_k = 0)))
  expect_false(acc(fake_fit(NaN, log_a = 0, log_b = 0, log_c = 0, log_k = 0)))
  ## A large but attainable positive log-likelihood stays admissible.
  expect_true(acc(fake_fit(91 * 3, log_a = 0, log_b = log(2),
                           log_c = 0, log_k = 0)))
})

test_that("the acceptance predicate rejects an exploded or vanished scale", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)
  expect_false(acc(fake_fit(-250, log_a = 0, log_b = log(2),
                            log_c = 0, log_k = log(1e9))))
  expect_false(acc(fake_fit(-250, log_a = 0, log_b = log(2),
                            log_c = 0, log_k = log(1e-9))))
})

test_that("non-finite estimates are inadmissible", {
  acc <- BetaDanish:::.bd_make_accept(n = 100)
  expect_false(acc(fake_fit(-250, log_a = NaN, log_b = 0,
                            log_c = 0, log_k = 0)))
  expect_false(acc(NULL))
})

test_that("the deterministic start grid is well formed", {
  for (sub in c(TRUE, FALSE)) {
    g <- BetaDanish:::.bd_default_starts(sub, k_base = 0.02)
    expect_gt(length(g), 10L)
    for (s in g) {
      expect_true(all(is.finite(s)))
      expect_true(all(grepl("^log_", names(s))))
      expect_equal(length(s), if (sub) 3L else 4L)
      ## Every start is itself admissible, so the grid can never seed the ridge.
      expect_true(all(exp(s[names(s) != "log_k"]) < 500))
    }
  }
})

