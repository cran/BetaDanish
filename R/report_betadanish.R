#' Create a Compact Report from a Beta-Danish Model Fit
#'
#' Collects the headline quantities from a fitted model into a small object with
#' a `print` method. Information criteria are computed from the fitted
#' log-likelihood via the `logLik` method, so they cannot fall out of step with
#' the fit.
#'
#' @param fit A fitted `"betadanish"` object.
#'
#' @return An object of class `"betadanish_report"`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- simulate_bd_data(120, a = 1, b = 3, c = 2, k = 0.5)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE)
#' report_betadanish(fit)
#' }
#'
#' @export
report_betadanish <- function(fit) {
  if (is.null(fit) || !inherits(fit, "betadanish"))
    stop("'fit' must be a fitted betadanish object.", call. = FALSE)

  out <- list(
    call         = fit$call,
    coefficients = fit$coefficients,
    submodel     = isTRUE(fit$submodel),
    logLik       = as.numeric(fit$logLik),
    npar         = .bd_or(fit$npar, length(fit$coefficients)),
    nobs         = .bd_or(fit$nobs, length(fit$data$time)),
    AIC          = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
    BIC          = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
    convergence  = fit$convergence,
    diagnostics  = fit$diagnostics
  )

  class(out) <- "betadanish_report"
  out
}

## Not `%||%`: base R gained that operator in 4.4.0, and defining it here would
## mask it for anyone attaching the package.
#' @noRd
.bd_or <- function(x, y) if (is.null(x)) y else x

#' @param x A `"betadanish_report"` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @rdname report_betadanish
#' @export
print.betadanish_report <- function(x, ...) {
  cat("Beta-Danish Model Report\n")
  cat("------------------------\n")
  cat("Model:          ",
      if (x$submodel) "3-parameter ED submodel (a = 1)" else "4-parameter Beta-Danish",
      "\n", sep = "")
  cat("Observations:   ", x$nobs, "  Parameters: ", x$npar, "\n", sep = "")
  cat("Log-likelihood: ", format(round(x$logLik, 4), nsmall = 4), "\n", sep = "")
  cat("AIC:            ", format(round(x$AIC, 4), nsmall = 4), "\n", sep = "")
  cat("BIC:            ", format(round(x$BIC, 4), nsmall = 4), "\n", sep = "")
  cat("Convergence:    ", x$convergence, "\n", sep = "")
  cat("\nEstimates:\n")
  print(round(x$coefficients, 4))

  d <- x$diagnostics
  if (!is.null(d)) {
    flags <- character(0)
    if (isTRUE(d$vcov_singular)) flags <- c(flags, "singular information matrix")
    if (isTRUE(d$near_b_ridge))  flags <- c(flags, "near the b = 1 ridge")
    if (is.finite(d$ac_correlation) && abs(d$ac_correlation) > 0.95)
      flags <- c(flags, "(a, c) confounded")
    if (!is.null(d$starts_rejected) && !is.na(d$starts_rejected) &&
        d$starts_rejected > 0)
      flags <- c(flags, sprintf("%d degenerate start(s) discarded",
                                d$starts_rejected))
    if (!is.null(d$loglik_spread) && !is.na(d$loglik_spread) &&
        d$loglik_spread > 2)
      flags <- c(flags, sprintf("local optima span %.2f log-lik units",
                                d$loglik_spread))
    if (length(flags))
      cat("\nDiagnostic flags: ", paste(flags, collapse = "; "),
          "\n  See the Identifiability section of ?fit_betadanish.\n", sep = "")
  }
  invisible(x)
}
