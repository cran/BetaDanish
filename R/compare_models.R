#' Goodness-of-Fit Statistics for Beta-Danish Models
#'
#' Information criteria (AIC, BIC, HQIC, AICC) and the Kolmogorov-Smirnov,
#' Cramer-von Mises and Anderson-Darling statistics for a fitted Beta-Danish
#' model, on one scale for complete and right-censored samples.
#'
#' @param object A fitted `betadanish` object.
#'
#' @return A list with components
#'   \describe{
#'     \item{`IC`}{Named vector `AIC`, `BIC`, `HQIC`, `AICC`.}
#'     \item{`KS_Statistic`}{Kolmogorov-Smirnov statistic \eqn{D}.}
#'     \item{`CvM_Statistic`}{Cramer-von Mises statistic \eqn{W^2}.}
#'     \item{`AD_Statistic`}{Anderson-Darling statistic \eqn{A^2}.}
#'     \item{`GOF`}{The three statistics as a named vector `c(D, W2, A2)`.}
#'     \item{`censored`}{Logical; whether the sample contains censoring.}
#'     \item{`at_boundary`}{Logical; whether the fit is on the edge of the
#'       parameter box, in which case every entry above is `NA`.}
#'   }
#'
#' @details
#' The statistics are computed on the probability-integral-transform scale
#' \eqn{u = F(t; \hat\theta)}, on which the fitted distribution function is
#' the identity, with the empirical distribution function replaced by the
#' Kaplan-Meier estimator \eqn{\hat F} of the transformed times:
#' \deqn{D = \sup_{0 \le u \le U} |\hat F(u) - u|, \quad
#'       W^2 = n\int_0^U (\hat F(u) - u)^2\,du, \quad
#'       A^2 = n\int_0^U \frac{(\hat F(u) - u)^2}{u(1-u)}\,du.}
#' The upper limit \eqn{U} is 1 when the largest observation is an event and
#' the transformed largest observation when it is censored. \eqn{\hat F} is a
#' step function, so the integrals are evaluated exactly. Without censoring
#' the three statistics reduce exactly to the classical formulae, so complete
#' and censored samples are measured on one scale.
#'
#' The supremum in \eqn{D} is taken at both ends of every step of
#' \eqn{\hat F}. Versions before 1.0.0 compared the two survivor curves only
#' at the right-hand end of each step and so could understate \eqn{D} for
#' censored data.
#'
#' A fit on the edge of the parameter box (see [fit_betadanish()]) is a
#' supremum rather than an attained maximum, and no information criteria or
#' goodness-of-fit statistics are reported for it.
#'
#' @references
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan.
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5,
#'                         censor_rate = 0.2, seed = 10)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 2)
#' gof_betadanish(fit)
#' }
#'
#' @export
gof_betadanish <- function(object) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  k <- length(object$coefficients)
  time   <- object$data$time
  status <- object$data$status
  n  <- length(time)
  ll <- as.numeric(object$logLik)
  boundary <- isTRUE(object$at_boundary)

  ic <- c(AIC  = 2 * k - 2 * ll,
          BIC  = k * log(n) - 2 * ll,
          HQIC = 2 * k * log(log(n)) - 2 * ll,
          AICC = 2 * k - 2 * ll + (2 * k * (k + 1)) / (n - k - 1))

  pars <- as.list(object$coefficients)
  if (isTRUE(object$submodel)) pars$a <- 1

  st <- c(D = NA_real_, W2 = NA_real_, A2 = NA_real_)
  if (boundary) {
    ic[] <- NA_real_
  } else {
    st <- .bd_gof_pit(time, status, function(t)
      pbetadanish(t, pars$a, pars$b, pars$c, pars$k))
  }

  list(
    IC            = ic,
    KS_Statistic  = unname(st[["D"]]),
    CvM_Statistic = unname(st[["W2"]]),
    AD_Statistic  = unname(st[["A2"]]),
    GOF           = st,
    censored      = any(status == 0),
    at_boundary   = boundary
  )
}

#' Compare Nested Beta-Danish Models
#'
#' Performs a Likelihood Ratio Test (LRT) between the 4-parameter full model 
#' and the 3-parameter submodel.
#'
#' @param full_model A fitted 4-parameter `betadanish` object.
#' @param sub_model A fitted 3-parameter `betadanish` object.
#' @return A data frame with the test statistic and p-value.
#' @export
compare_models <- function(full_model, sub_model) {
  if (full_model$submodel || !sub_model$submodel) {
    stop("Please provide the Full model as the first argument and the Submodel as the second.")
  }
  
  l1 <- full_model$logLik
  l0 <- sub_model$logLik
  
  if (l1 < l0) {
    warning("Full model log-likelihood is lower than submodel. Check convergence.")
  }
  if (isTRUE(full_model$at_boundary)) {
    warning("The four-parameter fit is on the edge of the parameter box, so ",
            "its log-likelihood is a supremum and the chi-squared reference ",
            "distribution is only approximate.", call. = FALSE)
  }
  
  test_stat <- max(0, -2 * (l0 - l1))
  p_val <- stats::pchisq(test_stat, df = 1, lower.tail = FALSE)
  
  res <- data.frame(
    Model = c("Submodel (3-param)", "Full Model (4-param)"),
    LogLik = c(l0, l1),
    Chisq = c(NA, test_stat),
    Df = c(NA, 1),
    `Pr(>Chisq)` = c(NA, p_val),
    check.names = FALSE
  )
  
  cat("Likelihood Ratio Test (a = 1 vs a != 1)\n\n")
  print(res, na.print = "")
  return(invisible(res))
}

#' Compare Beta-Danish with Standard Distributions
#'
#' Fits standard lifetime distributions (Weibull, Log-Normal, Log-Logistic, Gamma, 
#' Exponential) using the `flexsurv` package and compares them to the Beta-Danish fit.
#'
#' @param object A fitted `betadanish` object.
#' @return A ranked data frame of model comparisons.
#' @export
compare_distributions <- function(object) {
  if (!requireNamespace("flexsurv", quietly = TRUE)) {
    stop("The 'flexsurv' package is required for this function. Please install it.")
  }
  
  dat <- data.frame(time = object$data$time, status = object$data$status)
  form <- survival::Surv(time, status) ~ 1
  
  # Fit standard models safely
  safe_fit <- function(dist) {
    tryCatch({
      fit <- flexsurv::flexsurvreg(form, data = dat, dist = dist)
      ll <- fit$loglik
      k <- fit$npars
      n <- nrow(dat)
      c(LogLik = ll, AIC = -2*ll + 2*k, BIC = -2*ll + log(n)*k)
    }, error = function(e) c(LogLik = NA, AIC = NA, BIC = NA))
  }
  
  dists <- c("weibull", "lnorm", "llogis", "gamma", "exp")
  names(dists) <- c("Weibull", "Log-Normal", "Log-Logistic", "Gamma", "Exponential")
  
  results <- list()
  
  # Add Beta-Danish
  bd_gof <- gof_betadanish(object)
  results[["Beta-Danish"]] <- c(LogLik = object$logLik, AIC = unname(bd_gof$IC["AIC"]), BIC = unname(bd_gof$IC["BIC"]))
  
  # Add others
  for (nm in names(dists)) {
    results[[nm]] <- safe_fit(dists[nm])
  }
  
  res_df <- do.call(rbind, results)
  res_df <- as.data.frame(res_df)
  res_df$Distribution <- rownames(res_df)
  
  # Rank by AIC. A model without an AIC (a failed fit, or a Beta-Danish fit
  # on the edge of the parameter box) is listed last and not ranked.
  res_df <- res_df[order(res_df$AIC), c("Distribution", "LogLik", "AIC", "BIC")]
  res_df$Rank <- NA_integer_
  ok <- !is.na(res_df$AIC)
  res_df$Rank[ok] <- seq_len(sum(ok))
  rownames(res_df) <- NULL
  
  return(res_df)
}
