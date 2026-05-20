#' Read and Prepare Survival Data
#'
#' A helper function to read survival data from CSV or Excel files and prepare
#' it for analysis with the Beta-Danish package. It automatically handles missing
#' status columns by assuming all observations are complete (uncensored).
#'
#' @param file Character string specifying the path to the file.
#' @param time_col Character string specifying the name of the time/survival column.
#' @param status_col Character string specifying the name of the event/censoring
#'   indicator column. If `NULL` (default), the function assumes all observations
#'   are uncensored and creates a status column filled with 1s.
#' @param covar_cols Character vector specifying the names of covariate columns
#'   to keep. If `NULL` (default), no covariates are kept.
#'
#' @return A clean `data.frame` containing the `time`, `status`, and any
#'   specified covariates, ready to be passed to `fit_betadanish()`.
#'
#' @details
#' The function checks the file extension to determine how to read the data.
#' For `.xlsx` or `.xls` files, the `readxl` package must be installed.
#' Missing values (`NA`) in the specified columns will cause those rows to be
#' dropped with a warning.
#'
#' @export
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' example_data <- data.frame(
#'   survival_time = c(5, 8, 12, 16),
#'   status = c(1, 1, 0, 1)
#' )
#' write.csv(example_data, tmp, row.names = FALSE)
#' dat <- read_survival_data(tmp, time_col = "survival_time", status_col = "status")
#' unlink(tmp)
read_survival_data <- function(file, time_col, status_col = NULL, covar_cols = NULL) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Determine file type and read
  ext <- tolower(tools::file_ext(file))

  if (ext == "csv") {
    dat <- utils::read.csv(file, stringsAsFactors = FALSE)
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("The 'readxl' package is required to read Excel files. Please install it.")
    }
    dat <- as.data.frame(readxl::read_excel(file))
  } else {
    stop("Unsupported file extension. Please provide a .csv or .xlsx file.")
  }

  # Validate time column
  if (!time_col %in% colnames(dat)) {
    stop("Time column '", time_col, "' not found in the data.")
  }

  # Validate and handle status column
  if (is.null(status_col)) {
    message("No status column provided. Assuming all observations are uncensored (status = 1).")
    dat$status <- 1
    status_col <- "status"
  } else {
    if (!status_col %in% colnames(dat)) {
      stop("Status column '", status_col, "' not found in the data.")
    }
  }

  # Validate covariates
  if (!is.null(covar_cols)) {
    missing_covars <- setdiff(covar_cols, colnames(dat))
    if (length(missing_covars) > 0) {
      stop("Covariate column(s) not found: ", paste(missing_covars, collapse = ", "))
    }
  }

  # Subset data
  cols_to_keep <- c(time_col, status_col, covar_cols)
  clean_dat <- dat[, cols_to_keep, drop = FALSE]

  # Standardize core column names for easier internal use
  colnames(clean_dat)[1:2] <- c("time", "status")

  # Ensure numeric types
  clean_dat$time <- as.numeric(clean_dat$time)
  clean_dat$status <- as.numeric(clean_dat$status)

  # Handle NAs
  initial_rows <- nrow(clean_dat)
  clean_dat <- stats::na.omit(clean_dat)
  final_rows <- nrow(clean_dat)

  if (final_rows < initial_rows) {
    warning("Dropped ", initial_rows - final_rows, " rows due to missing values (NA).")
  }

  # Validate time and status values
  if (any(clean_dat$time <= 0)) {
    warning("Some time values are <= 0. Survival models typically require strictly positive times.")
  }

  if (!all(clean_dat$status %in% c(0, 1))) {
    stop("Status column contains values other than 0 and 1. Ensure 1 = event, 0 = censored.")
  }

  return(clean_dat)
}
