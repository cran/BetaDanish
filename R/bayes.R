#' Bayesian Estimation for the Beta-Danish Distribution
#'
#' Samples from the posterior of the Beta-Danish or Exponentiated Danish
#' parameters using a random-walk Metropolis sampler with vague
#' \eqn{\Gamma(0.01, 0.01)} priors on the positive parameters.
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators (1 = event, 0 = right-censored).
#' @param submodel Logical; \code{TRUE} for the 3-parameter ED submodel,
#'   \code{FALSE} for the 4-parameter full model.
#' @param burnin Burn-in iterations.
#' @param mcmc Post-burnin iterations.
#' @param tune Random-walk tuning parameter.
#' @param theta_init Optional starting values on the log scale.
#' @param seed Optional integer seed.
#' @param verbose Integer; passed to MCMCmetrop1R (0 = silent).
#'
#' @return An object of class \code{"bd_bayes"} with components
#'   \code{draws} (mcmc object), \code{summary}, \code{HPD},
#'   \code{submodel}, \code{call}.
#'
#' @details Requires \pkg{MCMCpack} and \pkg{coda} (Suggests).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' dat <- rbetadanish(100, a = 1.5, b = 2, c = 3, k = 0.5)
#' fit <- bayes_betadanish(time = dat, submodel = TRUE,
#'                         burnin = 500, mcmc = 1500)
#' fit$summary
#' }
#'
#' @export
bayes_betadanish <- function(time, status = NULL, submodel = TRUE,
                             burnin = 5000, mcmc = 15000, tune = 0.5,
                             theta_init = NULL, seed = NULL, verbose = 0) {
  if (!requireNamespace("MCMCpack", quietly = TRUE))
    stop("Bayesian estimation requires MCMCpack. Install with install.packages('MCMCpack').")
  if (!requireNamespace("coda", quietly = TRUE))
    stop("Bayesian estimation requires coda. Install with install.packages('coda').")
  time <- as.numeric(time)
  if (is.null(status)) status <- rep(1, length(time))
  status <- as.integer(status)
  if (length(status) != length(time))
    stop("time and status must have the same length.")
  dgamma_log <- function(x, shape, rate)
    stats::dgamma(x, shape = shape, rate = rate, log = TRUE)
  if (submodel) {
    logpost <- function(theta) {
      b <- exp(theta[1]); c <- exp(theta[2]); k <- exp(theta[3])
      if (any(!is.finite(c(b, c, k))) || any(c(b, c, k) <= 0)) return(-Inf)
      log_f <- dbetadanish(time, a = 1, b = b, c = c, k = k, log = TRUE)
      log_S <- pbetadanish(time, a = 1, b = b, c = c, k = k,
                           lower.tail = FALSE, log.p = TRUE)
      ll <- sum(status * log_f + (1 - status) * log_S)
      lp <- dgamma_log(b, 0.01, 0.01) + dgamma_log(c, 0.01, 0.01) +
            dgamma_log(k, 0.01, 0.01)
      ll + lp + theta[1] + theta[2] + theta[3]
    }
    if (is.null(theta_init)) theta_init <- log(c(2, 2, 1))
    names(theta_init) <- c("log_b", "log_c", "log_k")
  } else {
    logpost <- function(theta) {
      a <- exp(theta[1]); b <- exp(theta[2])
      c <- exp(theta[3]); k <- exp(theta[4])
      if (any(!is.finite(c(a, b, c, k))) || any(c(a, b, c, k) <= 0)) return(-Inf)
      log_f <- dbetadanish(time, a = a, b = b, c = c, k = k, log = TRUE)
      log_S <- pbetadanish(time, a = a, b = b, c = c, k = k,
                           lower.tail = FALSE, log.p = TRUE)
      ll <- sum(status * log_f + (1 - status) * log_S)
      lp <- dgamma_log(a, 0.01, 0.01) + dgamma_log(b, 0.01, 0.01) +
            dgamma_log(c, 0.01, 0.01) + dgamma_log(k, 0.01, 0.01)
      ll + lp + sum(theta)
    }
    if (is.null(theta_init)) theta_init <- log(c(1.5, 2, 2, 1))
    names(theta_init) <- c("log_a", "log_b", "log_c", "log_k")
  }
  if (!is.null(seed)) set.seed(seed)
  post <- MCMCpack::MCMCmetrop1R(logpost, theta.init = theta_init,
                                 burnin = burnin, mcmc = mcmc,
                                 tune = tune, verbose = verbose)
  draws_nat <- if (submodel)
    cbind(b = exp(post[, 1]), c = exp(post[, 2]), k = exp(post[, 3]))
  else
    cbind(a = exp(post[, 1]), b = exp(post[, 2]),
          c = exp(post[, 3]), k = exp(post[, 4]))
  draws_mcmc <- coda::as.mcmc(draws_nat)
  out <- list(
    draws    = draws_mcmc,
    summary  = summary(draws_mcmc),
    HPD      = coda::HPDinterval(draws_mcmc),
    submodel = submodel,
    call     = match.call()
  )
  class(out) <- "bd_bayes"
  out
}

#' @export
print.bd_bayes <- function(x, ...) {
  cat("\nBeta-Danish Bayesian Fit (",
      if (x$submodel) "3-Parameter ED Submodel"
      else "Full 4-Parameter Beta-Danish", ")\n", sep = "")
  cat("Random-walk Metropolis via MCMCpack::MCMCmetrop1R\n\n")
  print(x$summary)
  cat("\n95% HPD intervals:\n")
  print(x$HPD)
  invisible(x)
}
