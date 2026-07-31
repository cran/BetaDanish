#' Print Method for Beta-Danish Fit
#'
#' @param x An object of class `betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return Invisibly returns the input betadanish object. Called mainly for its side effect of printing the fitted model summary.
#' @export
print.betadanish <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nBeta-Danish Distribution Fit\n")
  cat("Model:", if(x$submodel) "3-Parameter Submodel (a=1)" else "Full 4-Parameter Model", "\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients:\n")
  print(round(x$coefficients, 4))
  cat("\n")
  invisible(x)
}

#' Summary Method for Beta-Danish Fit
#'
#' @param object An object of class `betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class summary.betadanish containing coefficient estimates, standard errors, test statistics, p-values, log-likelihood, and model selection criteria.
#' @export
summary.betadanish <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(pmax(diag(object$vcov), 0))
  z_val <- est / se
  p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)

  lower_95 <- est - 1.96 * se
  upper_95 <- est + 1.96 * se

  coef_table <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `Lower 95%` = lower_95,
    `Upper 95%` = upper_95,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )

  rownames(coef_table) <- names(est)

  k <- length(est)
  n <- length(object$data$time)
  aic <- 2 * k - 2 * object$logLik
  bic <- k * log(n) - 2 * object$logLik

  res <- list(
    call = object$call,
    submodel = object$submodel,
    coefficients = coef_table,
    logLik = object$logLik,
    aic = aic,
    bic = bic,
    convergence = object$convergence
  )

  class(res) <- "summary.betadanish"
  return(res)
}

#' Print Summary Method for Beta-Danish Fit
#'
#' @param x An object of class `summary.betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return Invisibly returns the input summary.betadanish object. Called mainly for its side effect of printing the coefficient table and fit statistics.
#' @export
print.summary.betadanish <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nBeta-Danish Distribution Fit\n")
  cat("Model:", if(x$submodel) "3-Parameter Submodel (a=1)" else "Full 4-Parameter Model", "\n\n")

  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("---\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n")
  cat("AIC:", round(x$aic, 4), " | BIC:", round(x$bic, 4), "\n")

  if (x$convergence != 0) {
    cat("\nWarning: Optimization may not have converged (Code:", x$convergence, ")\n")
  }
  invisible(x)
}

#' Extract Log-Likelihood
#'
#' @param object An object of class `betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return An object of class logLik containing the maximized log-likelihood value, with degrees of freedom and number of observations stored as attributes.
#' @export
logLik.betadanish <- function(object, ...) {
  val <- object$logLik
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- length(object$data$time)
  class(val) <- "logLik"
  return(val)
}

#' Extract Variance-Covariance Matrix
#'
#' @param object An object of class `betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return A numeric variance-covariance matrix for the estimated model parameters.
#' @export
vcov.betadanish <- function(object, ...) {
  return(object$vcov)
}

#' Extract Coefficients
#'
#' @param object An object of class `betadanish`.
#'
#' @param ... Further arguments passed to or from other methods.
#' @return A named numeric vector of maximum likelihood parameter estimates.
#' @export
coef.betadanish <- function(object, ...) {
  return(object$coefficients)
}
