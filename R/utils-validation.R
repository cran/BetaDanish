#' Extract Survival Data from Formula
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ 1)
#' @param data A data frame
#'
#' @return A list containing `time`, `status`, and `covariates` (if any)
#' @noRd
extract_surv_data <- function(formula, data) {
  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y <- stats::model.extract(mf, "response")

  if (!survival::is.Surv(Y)) {
    stop("The response must be a 'Surv' object. Example: Surv(time, status) ~ 1")
  }

  time <- Y[, 1]
  status <- Y[, 2] # 1 = event, 0 = censored

  # Extract covariates if present (excluding intercept)
  X <- stats::model.matrix(stats::terms(formula), mf)

  list(
    time = time,
    status = status,
    X = X,
    data_frame = mf
  )
}

#' Validate Beta-Danish Parameters
#'
#' @param a,b,c,k Numeric parameters
#' @return Logical TRUE if valid, FALSE otherwise
#' @noRd
check_positive_params <- function(a, b, c, k) {
  if (any(is.na(c(a, b, c, k)))) return(FALSE)
  if (any(c(a, b, c, k) <= 0)) return(FALSE)
  TRUE
}
