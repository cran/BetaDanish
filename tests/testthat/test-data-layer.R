test_that("brain_cancer is gone", {
  items <- utils::data(package = "BetaDanish")$results[, "Item"]
  expect_false("brain_cancer" %in% items)
})

test_that("guinea_pig has the documented shape", {
  data(guinea_pig, package = "BetaDanish", envir = environment())
  expect_s3_class(guinea_pig, "data.frame")
  expect_equal(nrow(guinea_pig), 72L)
  expect_named(guinea_pig, c("time", "status"))
  expect_true(all(guinea_pig$status == 1))
  expect_true(all(guinea_pig$time > 0))
  expect_equal(min(guinea_pig$time), 12)
  expect_equal(max(guinea_pig$time), 376)
})

test_that("guinea_pig gives a well-identified four-parameter fit", {
  skip_on_cran()
  data(guinea_pig, package = "BetaDanish", envir = environment())
  set.seed(1)
  fit <- suppressWarnings(
    fit_betadanish(survival::Surv(time, status) ~ 1, data = guinea_pig,
                   n_starts = 12, check_identifiability = FALSE))

  ## The thesis reports b-hat = 3.64 (SE 1.20), i.e. 2.2 SEs clear of the
  ## b = 1 ridge. Tolerances are loose: the point is that the optimum is
  ## interior and b is separated from one, not exact reproduction.
  expect_true(all(fit$coefficients > 0))
  expect_gt(fit$coefficients[["b"]], 1.5)
  expect_true(is.finite(fit$AIC))
})

test_that("the example CSVs ship and read back", {
  for (f in c("complete_sample.csv", "censored_sample.csv",
              "covariate_sample.csv", "competing_sample.csv")) {
    p <- system.file("extdata", f, package = "BetaDanish")
    expect_true(nzchar(p), info = f)
    expect_gt(nrow(utils::read.csv(p)), 20L)
  }
})

test_that("column names are guessed when not supplied", {
  p   <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, quiet = TRUE)
  ## Without covar_cols only the response is retained.
  expect_named(dat, c("time", "status"))
  rep <- attr(dat, "bd_data_report")
  expect_equal(rep$time_col, "time")
  expect_equal(rep$status_col, "status")
  expect_true("group" %in% rep$dropped_columns)
})

test_that("covar_cols = 'all' retains every non-response column", {
  p   <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, covar_cols = "all", quiet = TRUE)
  expect_true("group" %in% names(dat))
  expect_length(attr(dat, "bd_data_report")$dropped_columns, 0L)

  p2   <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
  dat2 <- read_survival_data(p2, covar_cols = "all", quiet = TRUE)
  expect_true(all(c("age", "thickness", "ulcer") %in% names(dat2)))
})

test_that("named covariates are retained and the rest reported as dropped", {
  p   <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, covar_cols = "age", quiet = TRUE)
  expect_named(dat, c("time", "status", "age"))
  expect_setequal(attr(dat, "bd_data_report")$dropped_columns,
                  c("thickness", "ulcer"))
})

test_that("dropped columns are announced unless quiet", {
  p <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")

  ## Several messages are emitted, so collect them all rather than relying on
  ## expect_message matching whichever comes first.
  msgs <- character(0)
  withCallingHandlers(
    read_survival_data(p),
    message = function(cond) {
      msgs <<- c(msgs, conditionMessage(cond))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("Not retained", msgs)))
  expect_true(any(grepl("group", msgs)))

  expect_silent(read_survival_data(p, quiet = TRUE))
})

test_that("an unguessable time column is a clear error", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(zzz = 1:5, qqq = 1L), tmp, row.names = FALSE)
  expect_error(read_survival_data(tmp, quiet = TRUE), "time column")
})

test_that("a competing-risks cause column is read", {
  p   <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
  dat <- read_survival_data(p, time_col = "time", cause_col = "cause",
                            quiet = TRUE)
  expect_true("cause" %in% names(dat))
  expect_true(all(dat$cause %in% 0:2))
})

test_that("semicolon-delimited European-format input works", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  writeLines(c("time;status", "1,5;1", "2,25;0", "3,75;1", "4,5;1", "5,25;1"),
             tmp)
  dat <- read_survival_data(tmp, sep = ";", dec = ",", quiet = TRUE)
  expect_equal(dat$time, c(1.5, 2.25, 3.75, 4.5, 5.25))
  expect_equal(dat$status, c(1, 0, 1, 1, 1))
})

test_that("a genuinely single-column file does not trigger a separator warning", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3, 4, 5)), tmp, row.names = FALSE)
  expect_warning(read_survival_data(tmp, quiet = TRUE), regexp = NA)
})

test_that("the separator heuristic flags only a mis-delimited header", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)

  ## Single column, but the header betrays a semicolon: flag it.
  writeLines(c("time;status", "1;1"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 1L, ","),
                 "another delimiter")

  ## Single column with a clean header: perfectly normal, stay quiet.
  writeLines(c("time", "1", "2"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 1L, ","), regexp = NA)

  ## Several columns were parsed, so the separator was clearly right.
  writeLines(c("time;status", "1;1"), tmp)
  expect_warning(BetaDanish:::.bd_check_separator(tmp, 3L, ","), regexp = NA)
})

test_that("a censoring indicator is rejected rather than guessed", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3, 4, 5), censor = c(0, 1, 0, 0, 1)),
                   tmp, row.names = FALSE)
  ## "censor" is not in the candidate list, so status is not guessed from it.
  dat <- read_survival_data(tmp, quiet = TRUE)
  expect_true(all(dat$status == 1))
})

test_that("grid-recorded times are detected and reported", {
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(3, 6, 9, 12, 18, 24, 36),
                              status = c(1, 1, 0, 1, 1, 0, 1)),
                   tmp, row.names = FALSE)
  dat <- read_survival_data(tmp, quiet = TRUE)
  expect_equal(attr(dat, "bd_data_report")$grid_step, 3)
})
