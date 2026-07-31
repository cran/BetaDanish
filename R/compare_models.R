#' Goodness-of-Fit Statistics for Beta-Danish Models
#'
#' Computes Information Criteria (AIC, BIC, HQIC, AICC) and the Kolmogorov-Smirnov 
#' (K-S) statistic for a fitted Beta-Danish model.
#'
#' @param object A fitted `betadanish` object.
#' @return A list containing the Information Criteria and the K-S statistic.
#' @export
gof_betadanish <- function(object) {
  k <- length(object$coefficients)
  n <- length(object$data$time)
  ll <- object$logLik
  
  aic <- 2 * k - 2 * ll
  bic <- k * log(n) - 2 * ll
  hqic <- 2 * k * log(log(n)) - 2 * ll
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  
  # K-S Statistic using Kaplan-Meier for censored data
  km <- survival::survfit(survival::Surv(object$data$time, object$data$status) ~ 1)
  km_surv <- km$surv
  km_time <- km$time
  
  pars <- as.list(object$coefficients)
  if (object$submodel) pars$a <- 1
  
  fitted_surv <- pbetadanish(km_time, pars$a, pars$b, pars$c, pars$k, lower.tail = FALSE)
  ks_stat <- max(abs(km_surv - fitted_surv), na.rm = TRUE)
  
  list(
    IC = c(AIC = aic, BIC = bic, HQIC = hqic, AICC = aicc),
    KS_Statistic = ks_stat
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
  
  # Rank by AIC
  res_df <- res_df[order(res_df$AIC), c("Distribution", "LogLik", "AIC", "BIC")]
  res_df$Rank <- 1:nrow(res_df)
  rownames(res_df) <- NULL
  
  return(res_df)
}
