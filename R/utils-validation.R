#' Extract Survival Data from Formula
#'
#' @param formula A survival formula (e.g. `Surv(time, status) ~ 1`).
#' @param data A data frame.
#'
#' @return A list with `time`, `status`, the design matrix `X`, and the model
#'   frame.
#' @noRd
extract_surv_data <- function(formula, data) {
  if (missing(data)) data <- environment(formula)

  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y  <- stats::model.extract(mf, "response")

  if (!survival::is.Surv(Y))
    stop("The response must be a 'Surv' object. Example: Surv(time, status) ~ 1",
         call. = FALSE)

  time   <- Y[, 1]
  status <- Y[, 2]

  if (any(!is.finite(time)) || any(time <= 0))
    stop("All survival times must be finite and strictly positive.", call. = FALSE)
  if (!all(status %in% c(0, 1)))
    stop("The event indicator must be coded 0 (censored) or 1 (event).",
         call. = FALSE)
  if (sum(status == 1) < 2L)
    stop("At least two events are needed to fit the model.", call. = FALSE)

  X <- stats::model.matrix(stats::terms(formula), mf)

  list(time = time, status = status, X = X, data_frame = mf)
}

#' Validate Beta-Danish Parameters
#'
#' Scalar validity check retained for the structural-property functions, which
#' take scalar parameters. The vectorised distribution functions validate
#' element-wise internally; see `.bd_bad()`.
#'
#' @param a,b,c,k Numeric parameters.
#' @return `TRUE` if every value is a finite strictly positive number.
#' @noRd
check_positive_params <- function(a, b, c, k) {
  v <- c(a, b, c, k)
  if (!is.numeric(v)) return(FALSE)
  if (anyNA(v)) return(FALSE)
  if (any(!is.finite(v))) return(FALSE)
  if (any(v <= 0)) return(FALSE)
  TRUE
}
