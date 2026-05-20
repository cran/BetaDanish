#' Robust Multi-Start Optimization using maxLik
#'
#' @param ll_fun Log-likelihood function to maximize
#' @param start_list A list of named numeric vectors for starting values
#' @param method Optimization method (default "BFGS")
#'
#' @return The best fitted maxLik object, or NULL if all fail
#' @noRd
optim_multistart <- function(ll_fun, start_list, method = "BFGS") {
  best_fit <- NULL
  best_ll <- -Inf

  for (start_vals in start_list) {
    fit <- tryCatch({
      maxLik::maxLik(logLik = ll_fun, start = start_vals, method = method)
    }, error = function(e) NULL)

    # Check if fit was successful and yielded a valid log-likelihood
    if (!is.null(fit) && !is.na(fit$maximum) && is.finite(fit$maximum)) {
      if (fit$maximum > best_ll) {
        best_ll <- fit$maximum
        best_fit <- fit
      }
    }
  }

  return(best_fit)
}

#' Log-Space Addition
#' Computes log(exp(x) + exp(y)) stably
#' @noRd
logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}
