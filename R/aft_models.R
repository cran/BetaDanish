#' Fit Beta-Danish AFT Regression Model
#'
#' Fits an Accelerated Failure Time (AFT) regression model using the
#' Exponentiated Danish (ED) kernel, that is the Beta-Danish distribution with a = 1.
#'
#' @param formula A survival formula (e.g., `Surv(time, status) ~ age + treatment`).
#' @param data A data frame containing the variables.
#' @param n_starts Integer; number of random starts for optimization.
#' @param method Optimization method passed to `maxLik`.
#'
#' @return An object of class `bd_aft`.
#'
#' @details
#' To ensure identifiability, the shape parameter `a` is fixed to 1. The scale
#' parameter `k` is linked to covariates via `k_i = exp(X_i %*% delta)`.
#' Positive coefficients in `delta` indicate a larger `k`, which corresponds
#' to shorter survival times (accelerated failure).
#'
#' @export
fit_bd_aft <- function(formula, data, n_starts = 10, method = "BFGS") {

  # Extract data
  surv_data <- extract_surv_data(formula, data)
  time <- surv_data$time
  status <- surv_data$status
  X <- surv_data$X

  # Define Log-Likelihood
  ll_fun <- function(pars) {
    log_b <- pars[1]
    log_c <- pars[2]
    delta <- pars[3:length(pars)]

    b <- exp(log_b)
    c <- exp(log_c)

    # k_i = exp(X %*% delta)
    eta_k <- as.numeric(X %*% delta)
    k_i <- exp(eta_k)

    # PDF and Survival
    log_f <- dbetadanish(time, a = 1, b = b, c = c, k = k_i, log = TRUE)
    log_S <- pbetadanish(time, a = 1, b = b, c = c, k = k_i, lower.tail = FALSE, log.p = TRUE)

    ll <- sum(status * log_f + (1 - status) * log_S)
    if (!is.finite(ll)) return(-1e10)
    return(ll)
  }

  # Multi-start initialization
  avg_t <- mean(time[status == 1], na.rm = TRUE)
  if (is.na(avg_t) || avg_t <= 0) avg_t <- mean(time, na.rm = TRUE)

  start_list <- list()
  for (i in 1:n_starts) {
    init_delta <- rep(0, ncol(X))
    init_delta[1] <- log(1 / avg_t) + stats::rnorm(1, 0, 0.5) # Jitter intercept

    start_list[[i]] <- c(
      log_b = log(stats::runif(1, 0.5, 3)),
      log_c = log(stats::runif(1, 0.5, 3)),
      init_delta
    )
    names(start_list[[i]])[3:length(start_list[[i]])] <- paste0("delta_", colnames(X))
  }

  # Optimize
  ## The master script's recorded degenerate case was an ED AFT fit that
  ## returned c = 296209.74. Shape explosion is rejected; large regression
  ## coefficients are not, following .ch6_separated: near-separation is
  ## information, not an error.
  fit <- optim_multistart(
    ll_fun, start_list, method = method,
    accept = .bd_make_accept(n = length(surv_data$time),
                             log_shape = c("log_b", "log_c"),
                             log_scale = character(0)))
  if (is.null(fit)) stop("AFT optimization failed to converge.")

  # Extract and format results
  est <- fit$estimate
  vcov_mat <- tryCatch(solve(-fit$hessian), error = function(e) matrix(NA, length(est), length(est)))
  rownames(vcov_mat) <- colnames(vcov_mat) <- names(est)

  out <- list(
    coefficients = est,
    logLik = fit$maximum,
    vcov = vcov_mat,
    convergence = fit$code,
    data = list(time = time, status = status, X = X),
    formula = formula,
    call = match.call(),
    submodel = TRUE
  )
  class(out) <- "bd_aft"
  return(out)
}
