#' Fit the Beta-Danish Distribution to Survival Data
#'
#' Fits the Beta-Danish distribution using Maximum Likelihood Estimation (MLE).
#' Supports both complete and right-censored data via `survival::Surv` objects.
#'
#' @param formula A formula object, with the response on the left of a `~` operator,
#'   and the terms on the right. The response must be a survival object as returned
#'   by the `Surv` function. Use `~ 1` for models without covariates.
#' @param data A data frame containing the variables named in the formula.
#' @param submodel Logical; if `TRUE`, fits the 3-parameter submodel by fixing `a = 1`.
#' @param n_starts Integer; the number of random starting points to use for the
#'   optimization to ensure global convergence. Default is 10.
#' @param method Character; the optimization method passed to `maxLik`. Default is "BFGS".
#'
#' @return An object of S3 class `"betadanish"`, containing the parameter estimates,
#'   log-likelihood, variance-covariance matrix, and convergence diagnostics.
#'
#' @details
#' The optimization is performed on the log-transformed parameters to strictly
#' enforce positivity constraints. The returned coefficients and variance-covariance
#' matrix are transformed back to the natural scale using the Delta method.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Simulate some data
#' set.seed(123)
#' sim_time <- rbetadanish(100, a = 1.5, b = 2, c = 3, k = 0.5)
#' sim_status <- sample(c(0, 1), 100, replace = TRUE, prob = c(0.2, 0.8))
#' dat <- data.frame(time = sim_time, status = sim_status)
#'
#' # Fit the 4-parameter model
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat)
#' summary(fit)
#'
#' # Fit the 3-parameter submodel
#' fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat, submodel = TRUE)
#' summary(fit_sub)
#' }
fit_betadanish <- function(formula, data, submodel = FALSE, n_starts = 10, method = "BFGS") {

  # Extract data
  surv_data <- extract_surv_data(formula, data)
  time <- surv_data$time
  status <- surv_data$status

  # Define Log-Likelihood Function (Optimizing over log-parameters)
  ll_fun <- function(pars) {
    # Extract and transform back to natural scale
    if (submodel) {
      a_par <- 1.0
      b_par <- exp(pars["log_b"])
      c_par <- exp(pars["log_c"])
      k_par <- exp(pars["log_k"])
    } else {
      a_par <- exp(pars["log_a"])
      b_par <- exp(pars["log_b"])
      c_par <- exp(pars["log_c"])
      k_par <- exp(pars["log_k"])
    }

    # Log-PDF for events (uncensored)
    lp <- dbetadanish(time, a_par, b_par, c_par, k_par, log = TRUE)

    # Log-Survival for censored
    ls <- pbetadanish(time, a_par, b_par, c_par, k_par, lower.tail = FALSE, log.p = TRUE)

    # Total Log-Likelihood
    loglik <- sum(status * lp + (1 - status) * ls)

    if (!is.finite(loglik)) return(-1e10)
    return(loglik)
  }

  # Generate Multi-Start Initial Values
  avg_t <- mean(time[status == 1], na.rm = TRUE)
  if (is.na(avg_t) || avg_t <= 0) avg_t <- mean(time, na.rm = TRUE)
  k_base <- 1 / avg_t

  start_list <- list()
  for (i in 1:n_starts) {
    if (submodel) {
      start_list[[i]] <- c(
        log_b = log(stats::runif(1, 0.5, 5)),
        log_c = log(stats::runif(1, 0.5, 5)),
        log_k = log(k_base * stats::runif(1, 0.5, 2))
      )
    } else {
      start_list[[i]] <- c(
        log_a = log(stats::runif(1, 0.5, 5)),
        log_b = log(stats::runif(1, 0.5, 5)),
        log_c = log(stats::runif(1, 0.5, 5)),
        log_k = log(k_base * stats::runif(1, 0.5, 2))
      )
    }
  }

  # Optimize
  fit <- optim_multistart(ll_fun, start_list, method = method)

  if (is.null(fit)) {
    stop("Optimization failed to converge for all starting points.")
  }

  # Transform estimates back to natural scale
  est_log <- fit$estimate
  est_nat <- exp(est_log)
  names(est_nat) <- gsub("log_", "", names(est_log))

  # Transform Variance-Covariance Matrix using Delta Method
  # V_nat = J * V_log * J^T, where J is the Jacobian (diagonal matrix of exp(est_log))
  vcov_log <- tryCatch(solve(-fit$hessian), error = function(e) matrix(NA, length(est_log), length(est_log)))
  J <- diag(est_nat, nrow = length(est_nat))
  vcov_nat <- J %*% vcov_log %*% t(J)
  rownames(vcov_nat) <- colnames(vcov_nat) <- names(est_nat)

  # Construct Output Object
  out <- list(
    coefficients = est_nat,
    logLik = fit$maximum,
    vcov = vcov_nat,
    convergence = fit$code,
    message = fit$message,
    submodel = submodel,
    data = list(time = time, status = status),
    formula = formula,
    call = match.call()
  )

  class(out) <- "betadanish"
  return(out)
}
