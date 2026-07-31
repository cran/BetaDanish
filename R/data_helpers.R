## Candidate column names used when time_col or status_col is left NULL.
## Deliberately excludes "censor" and "cens": those names are used for both
## codings in practice, and guessing wrong silently inverts every event.
.BD_TIME_NAMES   <- c("time", "times", "t", "survtime", "surv_time",
                      "survival_time", "survival_days", "duration", "lifetime",
                      "futime", "days", "months", "years")
.BD_STATUS_NAMES <- c("status", "event", "died", "death", "delta", "survstatus",
                      "failure", "observed")
.BD_CAUSE_NAMES  <- c("cause", "causes", "failtype", "fail_type", "event_type",
                      "risk", "failcause")

#' Match a Column Name Case-Insensitively Against a Candidate List
#' @noRd
.bd_guess_col <- function(nms, candidates) {
  lower <- tolower(nms)
  for (cand in candidates) {
    hit <- which(lower == cand)
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}

#' Detect Whether Times Look Recorded on a Coarse Grid
#'
#' Returns the inferred recording increment, or `NA` if the times do not look
#' grid-recorded. Whole-month or whole-day recording matters because the
#' point-density likelihood is not appropriate for coarsely rounded times.
#'
#' @noRd
.bd_grid_step <- function(t) {
  t <- t[is.finite(t) & t > 0]
  if (length(t) < 5L) return(NA_real_)
  u <- sort(unique(t))
  if (length(u) < 3L) return(NA_real_)
  d <- diff(u); d <- d[d > 0]
  if (!length(d)) return(NA_real_)
  g <- min(d)
  if (!is.finite(g) || g <= 0) return(NA_real_)
  r <- t / g
  if (max(abs(r - round(r))) < 1e-8) g else NA_real_
}

#' Warn Only When the Header Really Suggests a Different Separator
#'
#' A single-column result is perfectly normal for a file of bare times. It is
#' only suspicious when the header line also contains a tab, semicolon or pipe.
#'
#' @noRd
.bd_check_separator <- function(file, n_col, sep) {
  if (n_col > 1L || !identical(sep, ",")) return(invisible(NULL))
  first <- tryCatch(readLines(file, n = 1L, warn = FALSE),
                    error = function(e) character(0))
  if (!length(first)) return(invisible(NULL))
  if (grepl("[\t;|]", first))
    warning("Only one column was read, but the header contains another ",
            "delimiter. If the file is not comma separated, pass sep = ",
            "'\\t', sep = ';' or sep = '|'.", call. = FALSE)
  invisible(NULL)
}

#' Read and Prepare Survival Data
#'
#' Reads survival data from a delimited text file or an Excel workbook and
#' returns a clean data frame ready for [fit_betadanish()] and the regression
#' fitters. Columns are selected by name, covariates keep their original type,
#' and the result carries a report describing what was read.
#'
#' @param file Path to a `.csv`, `.txt`, `.tsv`, `.xls` or `.xlsx` file.
#' @param time_col Name of the time column. If `NULL` (default), the column is
#'   guessed from a list of common names and the choice is reported.
#' @param status_col Name of the event indicator (1 = event, 0 = censored). If
#'   `NULL`, guessed the same way; if no candidate is found, all observations
#'   are treated as uncensored.
#' @param covar_cols Covariate columns to retain. Either a character vector of
#'   column names, the string `"all"` to retain every column that is not part
#'   of the response, or `NULL` (default) to retain none. Whatever is dropped
#'   is named in a message unless `quiet = TRUE`.
#' @param cause_col Name of a competing-risks cause column, or `NULL`. When
#'   supplied, code 0 means censored and positive integers index the causes.
#' @param sep Field separator. Defaults to a comma. Use `"\\t"` for
#'   tab-delimited input, `";"` for semicolon-delimited.
#' @param dec Decimal mark. Use `","` for European-format numerics.
#' @param encoding File encoding passed to [utils::read.table()], for example
#'   `"UTF-8"` or `"latin1"`.
#' @param na.strings Strings to treat as missing.
#' @param drop_na Logical; drop incomplete rows. Default `TRUE`.
#' @param quiet Logical; suppress the informational messages.
#'
#' @return A data frame with columns `time`, `status`, an optional `cause`, and
#'   any retained covariates. The `"bd_data_report"` attribute is a list
#'   recording the file, which columns were used, rows read and kept, number of
#'   events, censoring proportion, retained and dropped columns, and the
#'   inferred recording grid.
#'
#' @details
#' Column guessing is deliberately conservative. Names meaning "censoring
#' indicator" are excluded from the candidate list, because `censor = 1` means
#' *censored* in some conventions and *observed* in others, and guessing wrong
#' would silently invert every event. Name the column explicitly if in doubt.
#'
#' A covariate that happens to be called `time`, `status` or `cause` is renamed
#' with a `_cov` suffix and a warning rather than colliding with the response.
#'
#' Excel input requires the `readxl` package.
#'
#' @seealso [fit_betadanish()], [fit_bd_aft()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' # A complete uncensored sample with a non-standard column name
#' f <- system.file("extdata", "complete_sample.csv", package = "BetaDanish")
#' dat <- read_survival_data(f, time_col = "survival_days", quiet = TRUE)
#' attr(dat, "bd_data_report")$rows_kept
#'
#' # Column names guessed, and every covariate retained
#' f2 <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
#' dat2 <- read_survival_data(f2, covar_cols = "all", quiet = TRUE)
#' names(dat2)
#'
#' # Competing risks
#' f3 <- system.file("extdata", "competing_sample.csv", package = "BetaDanish")
#' dat3 <- read_survival_data(f3, time_col = "time", cause_col = "cause",
#'                            quiet = TRUE)
#' table(dat3$cause)
read_survival_data <- function(file, time_col = NULL, status_col = NULL,
                               covar_cols = NULL, cause_col = NULL,
                               sep = ",", dec = ".", encoding = "unknown",
                               na.strings = c("NA", "", ".", "#N/A"),
                               drop_na = TRUE, quiet = FALSE) {

  if (!is.character(file) || length(file) != 1L)
    stop("'file' must be a single file path.", call. = FALSE)
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  say <- function(...) if (!isTRUE(quiet)) message(...)

  ## ---- read -----------------------------------------------------------------
  ext <- tolower(tools::file_ext(file))
  dat <- if (ext %in% c("csv", "txt", "tsv", "dat")) {
    if (ext == "tsv" && identical(sep, ",")) sep <- "\t"
    utils::read.table(file, header = TRUE, sep = sep, dec = dec,
                      na.strings = na.strings, encoding = encoding,
                      stringsAsFactors = FALSE, check.names = TRUE,
                      comment.char = "")
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Excel files requires the 'readxl' package.", call. = FALSE)
    as.data.frame(readxl::read_excel(file))
  } else {
    stop("Unsupported file extension '", ext, "'. Supply a .csv, .txt, .tsv, ",
         ".xls or .xlsx file.", call. = FALSE)
  }

  if (!nrow(dat)) stop("The file contains no data rows.", call. = FALSE)
  .bd_check_separator(file, ncol(dat), sep)

  nms <- names(dat)

  ## ---- resolve the time column ---------------------------------------------
  if (is.null(time_col)) {
    time_col <- .bd_guess_col(nms, .BD_TIME_NAMES)
    if (is.null(time_col))
      stop("Could not identify a time column. Available columns: ",
           paste(nms, collapse = ", "), ". Pass time_col explicitly.",
           call. = FALSE)
    say("Using '", time_col, "' as the time column.")
  }
  if (!time_col %in% nms)
    stop("Time column '", time_col, "' not found. Available columns: ",
         paste(nms, collapse = ", "), call. = FALSE)

  clean <- data.frame(time = suppressWarnings(as.numeric(dat[[time_col]])))
  if (all(is.na(clean$time)))
    stop("Column '", time_col, "' contains no numeric values. Check 'dec' if ",
         "the file uses a comma as the decimal mark.", call. = FALSE)

  ## ---- resolve the status column -------------------------------------------
  guessed_status <- FALSE
  if (is.null(status_col)) {
    status_col     <- .bd_guess_col(nms, .BD_STATUS_NAMES)
    guessed_status <- !is.null(status_col)
  }
  if (is.null(status_col)) {
    say("No status column found; treating all observations as uncensored.")
    clean$status <- 1
  } else {
    if (!status_col %in% nms)
      stop("Status column '", status_col, "' not found.", call. = FALSE)
    if (guessed_status)
      say("Using '", status_col, "' as the event indicator ",
          "(1 = event, 0 = censored).")
    clean$status <- suppressWarnings(as.numeric(dat[[status_col]]))
  }

  ## ---- optional cause column ------------------------------------------------
  if (is.null(cause_col)) {
    cc <- .bd_guess_col(nms, .BD_CAUSE_NAMES)
    if (!is.null(cc) && !identical(cc, status_col)) {
      cause_col <- cc
      say("Using '", cause_col, "' as the competing-risks cause column.")
    }
  }
  if (!is.null(cause_col)) {
    if (!cause_col %in% nms)
      stop("Cause column '", cause_col, "' not found.", call. = FALSE)
    clean$cause <- suppressWarnings(as.integer(dat[[cause_col]]))
    if (anyNA(clean$cause) && !all(is.na(dat[[cause_col]])))
      warning("Some cause codes could not be read as integers.", call. = FALSE)
  }

  ## ---- covariates -----------------------------------------------------------
  response_cols <- c(time_col, status_col, cause_col)

  if (identical(covar_cols, "all")) {
    covar_cols <- setdiff(nms, response_cols)
    if (length(covar_cols))
      say("Retaining ", length(covar_cols), " covariate(s): ",
          paste(covar_cols, collapse = ", "), ".")
  }

  if (length(covar_cols)) {
    missing_cov <- setdiff(covar_cols, nms)
    if (length(missing_cov))
      stop("Covariate column(s) not found: ", paste(missing_cov, collapse = ", "),
           call. = FALSE)
    extra    <- dat[, covar_cols, drop = FALSE]
    reserved <- c("time", "status", "cause")
    clash    <- intersect(names(extra), reserved)
    if (length(clash)) {
      warning("Covariate(s) named ", paste(clash, collapse = ", "),
              " collide with the response; renamed with a '_cov' suffix.",
              call. = FALSE)
      names(extra)[names(extra) %in% clash] <-
        paste0(names(extra)[names(extra) %in% clash], "_cov")
    }
    clean <- cbind(clean, extra)
  }

  dropped_cols <- setdiff(nms, c(response_cols, covar_cols))
  if (length(dropped_cols))
    say("Not retained: ", paste(dropped_cols, collapse = ", "),
        ". Pass covar_cols to keep specific columns, or covar_cols = \"all\".")

  ## ---- validate -------------------------------------------------------------
  n_read <- nrow(clean)
  if (isTRUE(drop_na)) {
    clean <- clean[stats::complete.cases(clean), , drop = FALSE]
    if (nrow(clean) < n_read)
      warning("Dropped ", n_read - nrow(clean), " row(s) with missing values.",
              call. = FALSE)
  }
  n_kept <- nrow(clean)
  if (!n_kept)
    stop("No complete rows remain after removing missing values.", call. = FALSE)

  if (any(clean$time <= 0, na.rm = TRUE))
    warning("Some times are <= 0; survival models require strictly positive ",
            "times.", call. = FALSE)
  if (!all(clean$status %in% c(0, 1)))
    stop("The status column contains values other than 0 and 1. Recode so ",
         "that 1 = event and 0 = censored. If '", status_col, "' is a ",
         "censoring indicator, invert it first.", call. = FALSE)

  grid <- .bd_grid_step(clean$time)
  if (!is.na(grid))
    say("Times appear recorded on a grid of ", format(grid),
        ". For coarsely rounded data the point-density likelihood ",
        "understates uncertainty.")

  row.names(clean) <- NULL
  attr(clean, "bd_data_report") <- list(
    file              = normalizePath(file, mustWork = FALSE),
    time_col          = time_col,
    status_col        = status_col,
    cause_col         = cause_col,
    available_columns = nms,
    covariates        = setdiff(names(clean), c("time", "status", "cause")),
    dropped_columns   = dropped_cols,
    rows_read         = n_read,
    rows_kept         = n_kept,
    rows_dropped      = n_read - n_kept,
    n_events          = sum(clean$status == 1),
    censoring_prop    = mean(clean$status == 0),
    grid_step         = grid
  )
  clean
}
