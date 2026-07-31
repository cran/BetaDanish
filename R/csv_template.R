#' Write a Skeleton CSV for BetaDanish
#'
#' Writes a small, correctly shaped CSV showing the layout [bd_analyze_csv()]
#' and [read_survival_data()] expect. Fill it in with your own data, or use it
#' as a reference for renaming the columns of an existing file.
#'
#' @param path Destination file path. Use [tempfile()] if you only want to
#'   inspect the result.
#' @param type Which layout to write:
#'   `"univariate"` gives `time` and `status`;
#'   `"complete"` gives `time` alone, for an uncensored sample;
#'   `"covariate"` adds two example covariate columns;
#'   `"competing"` gives `time` and `cause`.
#' @param n Number of illustrative rows. Default 10.
#' @param overwrite Logical; overwrite an existing file. Default `FALSE`.
#'
#' @return The path, invisibly.
#'
#' @details
#' The illustrative values are plausible but arbitrary. `status` is coded
#' 1 = event observed, 0 = right-censored. For `"competing"`, `cause` is coded
#' 0 = censored and 1, 2, ... for the competing causes.
#'
#' @seealso [bd_analyze_csv()], [read_survival_data()]
#'
#' @export
#'
#' @examples
#' p <- bd_csv_template(tempfile(fileext = ".csv"), type = "covariate", n = 5)
#' read.csv(p)
#' unlink(p)
bd_csv_template <- function(path,
                            type = c("univariate", "complete", "covariate",
                                     "competing"),
                            n = 10, overwrite = FALSE) {
  type <- match.arg(type)
  if (!is.character(path) || length(path) != 1L)
    stop("'path' must be a single file path.", call. = FALSE)
  if (file.exists(path) && !isTRUE(overwrite))
    stop("'", path, "' already exists. Pass overwrite = TRUE to replace it.",
         call. = FALSE)
  n <- as.integer(n)
  if (is.na(n) || n < 1L) stop("'n' must be a positive integer.", call. = FALSE)

  tm <- round(seq(4, 4 + 6 * (n - 1), length.out = n) +
                stats::runif(n, -1, 1), 1)
  tm <- pmax(tm, 0.5)

  out <- switch(
    type,
    univariate = data.frame(time = tm,
                            status = rep_len(c(1L, 1L, 0L), n)),
    complete   = data.frame(time = tm),
    covariate  = data.frame(time = tm,
                            status = rep_len(c(1L, 1L, 0L), n),
                            age = rep_len(c(45L, 62L, 51L, 70L), n),
                            group = rep_len(c("treated", "control"), n)),
    competing  = data.frame(time = tm,
                            cause = rep_len(c(1L, 2L, 0L), n))
  )

  utils::write.csv(out, path, row.names = FALSE)
  message("Wrote a '", type, "' template with ", n, " row(s) to: ", path)
  invisible(path)
}
