#' Create a compact report from a BetaDanish model fit
#'
#' @param fit A fitted BetaDanish model object.
#' @return A list containing model summary information.
#' @export
report_betadanish <- function(fit) {
  if (is.null(fit)) stop('fit must be a fitted BetaDanish object')

  out <- list(
    call = fit$call,
    coefficients = fit$coefficients,
    logLik = fit$logLik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    convergence = fit$convergence
  )

  class(out) <- 'betadanish_report'
  out
}

#' @export
print.betadanish_report <- function(x, ...) {
  cat('BetaDanish Model Report\n')
  cat('-----------------------\n')
  if (!is.null(x$AIC)) cat('AIC:', x$AIC, '\n')
  if (!is.null(x$BIC)) cat('BIC:', x$BIC, '\n')
  if (!is.null(x$logLik)) cat('Log-likelihood:', x$logLik, '\n')
  if (!is.null(x$convergence)) cat('Convergence:', x$convergence, '\n')
  invisible(x)
}
