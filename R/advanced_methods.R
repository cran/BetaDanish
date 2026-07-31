## S3 methods for the AFT, cure and competing-risks fits.
##
## Shape parameters are estimated on the log scale. They are reported here on
## the natural scale with delta-method standard errors and an exponentiated
## log-scale confidence interval, which respects positivity. No Wald test is
## reported for them: the hypothesis "b = 0" is outside the parameter space, so
## a z statistic against zero would be meaningless. Regression coefficients are
## already on their natural scale and are reported with the usual Wald test.

#' Coefficient Table for a Shape/Regression Split
#' @noRd
.bd_split_coef <- function(est, se) {
  shape_idx <- grep("^log_", names(est))
  reg_idx   <- setdiff(seq_along(est), shape_idx)

  shape <- NULL
  if (length(shape_idx)) {
    l <- est[shape_idx]; s <- se[shape_idx]
    shape <- cbind(Estimate     = exp(l),
                   `Std. Error` = exp(l) * s,
                   `Lower 95%`  = exp(l - stats::qnorm(0.975) * s),
                   `Upper 95%`  = exp(l + stats::qnorm(0.975) * s))
    rownames(shape) <- sub("^log_", "", names(l))
  }

  reg <- NULL
  if (length(reg_idx)) {
    e <- est[reg_idx]; s <- se[reg_idx]
    z <- e / s
    reg <- cbind(Estimate     = e,
                 `Std. Error` = s,
                 `z value`    = z,
                 `Pr(>|z|)`   = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
    rownames(reg) <- names(e)
  }
  list(shape = shape, regression = reg)
}

.bd_print_split <- function(x) {
  if (!is.null(x$shape)) {
    cat("Shape parameters (natural scale, delta-method SE):\n")
    print(round(x$shape, 4))
    cat("\n")
  }
  if (!is.null(x$regression)) {
    cat("Regression coefficients (log-scale link):\n")
    stats::printCoefmat(x$regression, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  }
  invisible(NULL)
}

.bd_se <- function(object) {
  se <- sqrt(pmax(diag(object$vcov), 0))
  names(se) <- names(object$coefficients)
  se
}

#' @export
print.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model (Exponentiated Danish kernel, a = 1)\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients (optimisation scale):\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_aft <- function(object, ...) {
  res <- list(call = object$call,
              tables = .bd_split_coef(object$coefficients, .bd_se(object)),
              logLik = object$logLik)
  class(res) <- "summary.bd_aft"
  res
}

#' @export
print.summary.bd_aft <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish AFT Model (Exponentiated Danish kernel, a = 1)\n\n")
  .bd_print_split(x$tables)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}

#' @export
print.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n", sep = "")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Coefficients (optimisation scale):\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_cure <- function(object, ...) {
  res <- list(call = object$call, type = object$type,
              tables = .bd_split_coef(object$coefficients, .bd_se(object)),
              logLik = object$logLik)
  class(res) <- "summary.bd_cure"
  res
}

#' @export
print.summary.bd_cure <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Cure Model (", x$type, ")\n\n", sep = "")
  .bd_print_split(x$tables)
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  invisible(x)
}

#' @export
print.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n")
  cat("Log-Likelihood:", round(x$logLik, 4), "\n\n")
  cat("Cause-specific estimates:\n")
  print(round(x$coefficients, 4)); cat("\n")
  invisible(x)
}

#' @export
summary.bd_competing <- function(object, ...) {
  res <- list(call = object$call, coefficients = object$coefficients,
              se = object$se, logLik = object$logLik,
              causes = object$causes)
  class(res) <- "summary.bd_competing"
  res
}

#' @export
print.summary.bd_competing <- function(x, ...) {
  cat("\nCall:\n"); print(x$call)
  cat("\nBeta-Danish Competing Risks Model\n\n")
  cat("Cause-specific estimates (natural scale):\n"); print(round(x$coefficients, 4))
  cat("\nStandard errors:\n"); print(round(x$se, 4))
  cat("---\nLog-Likelihood:", round(x$logLik, 4), "\n")
  cat("\nNote: the cause-specific marginals rest on an assumption of\n")
  cat("independent latent failure times, which is an identifying\n")
  cat("assumption and not testable from the observed data. See\n")
  cat("?fit_bd_competing.\n")
  invisible(x)
}
