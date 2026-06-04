test_that("read_survival_data handles CSVs and missing values", {
  tmp_file <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(surv_time = c(1, 2, NA, 4),
               event     = c(1, 0, 1, 1)),
    tmp_file, row.names = FALSE
  )
  # The function correctly warns about the dropped NA row;
  # the test acknowledges and silences that warning.
  expect_warning(
    d <- read_survival_data(tmp_file,
                            time_col   = "surv_time",
                            status_col = "event"),
    regexp = "missing values"
  )
  expect_equal(nrow(d), 3)
  expect_setequal(colnames(d), c("time", "status"))
  unlink(tmp_file)
})

test_that("read_survival_data assumes uncensored when status_col is NULL", {
  tmp_file <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(surv_time = c(1, 2, 3, 4)),
    tmp_file, row.names = FALSE
  )
  expect_message(
    d <- read_survival_data(tmp_file, time_col = "surv_time"),
    regexp = "uncensored"
  )
  expect_true(all(d$status == 1))
  unlink(tmp_file)
})
