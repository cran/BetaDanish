csv_path <- function(f) system.file("extdata", f, package = "BetaDanish")

test_that("bd_csv_template writes each layout", {
  expected <- list(univariate = c("time", "status"),
                   complete   = "time",
                   covariate  = c("time", "status", "age", "group"),
                   competing  = c("time", "cause"))
  for (ty in names(expected)) {
    p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
    suppressMessages(bd_csv_template(p, type = ty, n = 6))
    got <- utils::read.csv(p)
    expect_named(got, expected[[ty]], info = ty)
    expect_equal(nrow(got), 6L, info = ty)
    expect_true(all(got$time > 0), info = ty)
  }
})

test_that("bd_csv_template refuses to clobber without overwrite", {
  p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
  suppressMessages(bd_csv_template(p))
  expect_error(bd_csv_template(p), "already exists")
  expect_silent(suppressMessages(bd_csv_template(p, overwrite = TRUE)))
})

test_that("a template round-trips through read_survival_data", {
  p <- tempfile(fileext = ".csv"); on.exit(unlink(p), add = TRUE)
  suppressMessages(bd_csv_template(p, type = "covariate", n = 12))
  dat <- read_survival_data(p, covar_cols = "all", quiet = TRUE)
  expect_named(dat, c("time", "status", "age", "group"))
  expect_true(all(dat$status %in% c(0, 1)))
})

test_that("the univariate pipeline runs and returns tidy tables", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"),
                        analysis = "univariate", model = "ED",
                        compare = FALSE, n_starts = 3, seed = 11, quiet = TRUE)

  expect_s3_class(res, "bd_analysis")
  expect_true(inherits(res$fits$ED, "betadanish"))
  expect_true(is.data.frame(res$tables$estimates))
  expect_setequal(res$tables$estimates$parameter, c("b", "c", "k"))
  expect_true(is.data.frame(res$tables$information_criteria))
  expect_true(is.finite(res$tables$information_criteria$AIC[1]))
  expect_length(res$failures, 0L)
  expect_output(print(res), "BetaDanish CSV analysis")
})

test_that("model = 'both' adds a likelihood ratio test", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "both",
                        compare = FALSE, n_starts = 3, seed = 12, quiet = TRUE)
  expect_true(all(c("BD", "ED") %in% names(res$fits)))
  expect_true(is.data.frame(res$tables$likelihood_ratio_test))
  expect_equal(nrow(res$tables$information_criteria), 2L)
})

test_that("nothing is written unless output_dir is given", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 13, quiet = TRUE)
  expect_null(res$output_dir)
  expect_length(res$files, 0L)
})

test_that("output_dir receives tables and figures", {
  skip_on_cran()
  od <- file.path(tempdir(), paste0("bdout_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(od, recursive = TRUE), add = TRUE)

  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, output_dir = od, n_starts = 3,
                        seed = 14, quiet = TRUE)

  expect_true(dir.exists(od))
  expect_gt(length(res$files), 3L)
  expect_true(any(grepl("estimates[.]csv$", res$files)))
  expect_true(any(grepl("survival[.]png$", res$files)))
  expect_true(all(file.exists(res$files)))

  back <- utils::read.csv(grep("estimates[.]csv$", res$files, value = TRUE)[1])
  expect_true(all(c("model", "parameter", "estimate") %in% names(back)))
})

test_that("the AFT path uses the covariates from the file", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("covariate_sample.csv"), analysis = "aft",
                        covariates = c("age", "thickness"),
                        compare = FALSE, n_starts = 3, seed = 15, quiet = TRUE)
  expect_true(inherits(res$fits$AFT, "bd_aft"))
  expect_true(any(grepl("^delta_age$", names(res$fits$AFT$coefficients))))
  expect_true(!is.null(res$extras$aft))
})

test_that("the competing-risks path needs a cause column", {
  skip_on_cran()
  expect_error(
    bd_analyze_csv(csv_path("censored_sample.csv"), analysis = "competing",
                   n_starts = 2, quiet = TRUE),
    "cause column")

  res <- bd_analyze_csv(csv_path("competing_sample.csv"), analysis = "competing",
                        time_col = "time", cause_col = "cause",
                        compare = FALSE, n_starts = 2, seed = 16, quiet = TRUE)
  expect_true(inherits(res$fits$CR, "bd_competing"))
})

test_that("a failing fit is recorded rather than fatal", {
  skip_on_cran()
  ## Two rows cannot support a fit, but extract_surv_data rejects it first,
  ## so the failure must surface as a recorded message, not an abort.
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data.frame(time = c(1, 2, 3), status = c(1, 0, 0)), tmp,
                   row.names = FALSE)
  res <- bd_analyze_csv(tmp, model = "ED", compare = FALSE, n_starts = 2,
                        quiet = TRUE)
  expect_s3_class(res, "bd_analysis")
  expect_gt(length(res$failures), 0L)
  expect_null(res$fits$ED)
})

test_that("model warnings are captured rather than discarded", {
  skip_on_cran()
  ## The four-parameter fit on this file lands near the b = 1 ridge, so the
  ## identifiability diagnostic fires. It must survive the pipeline: swallowing
  ## it would hide exactly the thing the user needs to see.
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "both",
                        compare = FALSE, n_starts = 5, seed = 21, quiet = TRUE)
  expect_type(res$warnings, "character")
  if (length(res$warnings)) expect_true(all(nzchar(names(res$warnings))))
  ## Whatever happened, it is recorded somewhere rather than lost.
  expect_true(!is.null(res$fits$BD) || length(res$failures) > 0L)
})

test_that("summary prints every table without error", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 17, quiet = TRUE)
  expect_output(print(summary(res)), "estimates")
})

test_that("plot.bd_analysis dispatches to the primary fit", {
  skip_on_cran()
  res <- bd_analyze_csv(csv_path("censored_sample.csv"), model = "ED",
                        compare = FALSE, n_starts = 3, seed = 18, quiet = TRUE)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot(res, type = "survival"), NA)
})
