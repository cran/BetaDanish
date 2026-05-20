test_that("read_survival_data handles CSVs and missing values", {
  # Create a temporary CSV file
  tmp_file <- tempfile(fileext = ".csv")
  test_data <- data.frame(
    surv_time = c(10, 20, NA, 40),
    event_stat = c(1, 0, 1, 1),
    age = c(50, 60, 70, 80)
  )
  write.csv(test_data, tmp_file, row.names = FALSE)

  # Test reading with covariates
  expect_warning(
    clean_dat <- read_survival_data(tmp_file, time_col = "surv_time", status_col = "event_stat", covar_cols = "age"),
    "Dropped 1 rows due to missing values"
  )

  expect_equal(nrow(clean_dat), 3)
  expect_equal(colnames(clean_dat), c("time", "status", "age"))

  # Test reading without status (assumes all 1s)
  clean_dat2 <- read_survival_data(tmp_file, time_col = "surv_time")
  expect_equal(clean_dat2$status, c(1, 1, 1)) # NA row was dropped

  unlink(tmp_file)
})
