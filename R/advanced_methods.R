#' @export
print.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients:\n"); print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}
#' @export
summary.bd_aft <- function(object, ...) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- object$coefficients / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  res <- list(call = object$call, coefficients = cbind(Estimate=object$coefficients, `Std. Error`=se, `z value`=z, `Pr(>|z|)`=p), logLik = object$logLik)
  class(res) <- "summary.bd_aft"; return(res)
}
#' @export
print.summary.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model\n\n")
  stats::printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}
#' @export
print.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n", sep="")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients:\n"); print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}
#' @export
summary.bd_cure <- function(object, ...) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- object$coefficients / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  res <- list(call = object$call, type = object$type, coefficients = cbind(Estimate=object$coefficients, `Std. Error`=se, `z value`=z, `Pr(>|z|)`=p), logLik = object$logLik)
  class(res) <- "summary.bd_cure"; return(res)
}
#' @export
print.summary.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n\n", sep="")
  stats::printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}
#' @export
print.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients:\n"); print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}
#' @export
summary.bd_competing <- function(object, ...) {
  res <- list(call = object$call, coefficients = object$coefficients, se = object$se, logLik = object$logLik)
  class(res) <- "summary.bd_competing"; return(res)
}
#' @export
print.summary.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n\n")
  cat("Coefficients:\n"); print(x$coefficients)
  cat("\nStandard Errors:\n"); print(x$se)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}
