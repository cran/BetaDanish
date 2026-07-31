#' Fit Beta-Danish Cure Models
#'
#' Fits mixture and promotion-time (non-mixture) cure models using the 
#' Beta-Danish AFT baseline.
#'
#' @param formula_aft A formula for the latency component (e.g., `Surv(time, status) ~ age`).
#' @param formula_cure A one-sided formula for the incidence/cure component (e.g., `~ treatment`).
#' @param data A data frame containing the variables.
#' @param type Character; either `"mixture"` or `"promotion"`.
#' @param n_starts Integer; number of random starts for optimization.
#' @param method Optimization method passed to `maxLik`.
#'
#' @return An object of class `bd_cure`.
#'
#' @details
#' In the **mixture** model, the population is split into susceptible and cured 
#' fractions. The susceptible probability is modeled via logistic regression: 
#' `pi = exp(Z %*% gamma) / (1 + exp(Z %*% gamma))`.
#' 
#' In the **promotion-time** (non-mixture) model, the cure fraction is derived 
#' from a latent Poisson process of clonogenic cells: `theta = exp(Z %*% gamma)`. 
#' The cure fraction is `exp(-theta)`.
#'
#' @export
fit_bd_cure <- function(formula_aft, formula_cure, data, type = c("mixture", "promotion"), n_starts = 10, method = "BFGS") {
  type <- match.arg(type)
  
  # Ensure complete cases across both formulas
  vars_needed <- unique(c(all.vars(formula_aft), all.vars(formula_cure)))
  data <- data[stats::complete.cases(data[, vars_needed, drop = FALSE]), , drop = FALSE]
  
  # Extract data
  surv_data <- extract_surv_data(formula_aft, data)
  time <- surv_data$time
  status <- surv_data$status
  X <- surv_data$X # Latency covariates
  
  # Extract cure covariates from the ORIGINAL data
  mf_cure <- stats::model.frame(formula_cure, data, na.action = stats::na.pass)
  Z <- stats::model.matrix(formula_cure, mf_cure)
  
  # Define Log-Likelihood
  ll_fun <- function(pars) {
    log_b <- pars[1]
    log_c <- pars[2]
    delta <- pars[3:(2 + ncol(X))]
    gamma <- pars[(3 + ncol(X)):length(pars)]
    
    b <- exp(log_b)
    c <- exp(log_c)
    
    # Latency scale
    k_i <- exp(as.numeric(X %*% delta))
    
    # Baseline functions
    log_f <- dbetadanish(time, a = 1, b = b, c = c, k = k_i, log = TRUE)
    log_S <- pbetadanish(time, a = 1, b = b, c = c, k = k_i, lower.tail = FALSE, log.p = TRUE)
    F_u <- pbetadanish(time, a = 1, b = b, c = c, k = k_i, lower.tail = TRUE, log.p = FALSE)
    
    if (type == "mixture") {
      # pi = susceptible proportion
      eta_pi <- as.numeric(Z %*% gamma)
      log_pi <- -log1p(exp(-eta_pi))       # log(pi)
      log_1m_pi <- -log1p(exp(eta_pi))     # log(1 - pi)
      
      # Event: log(pi * f)
      ll_event <- log_pi + log_f
      
      # Censored: log((1 - pi) + pi * S)
      ll_cens <- logspace_add(log_1m_pi, log_pi + log_S)
      
      ll <- sum(status * ll_event + (1 - status) * ll_cens)
      
    } else {
      # Promotion-time
      theta_i <- exp(as.numeric(Z %*% gamma))
      
      # Event: log(theta * f * exp(-theta * F))
      ll_event <- log(theta_i) + log_f - theta_i * F_u
      
      # Censored: log(exp(-theta * F))
      ll_cens <- -theta_i * F_u
      
      ll <- sum(status * ll_event + (1 - status) * ll_cens)
    }
    
    if (!is.finite(ll)) return(-1e10)
    return(ll)
  }
  
  # Multi-start initialization
  km <- survival::survfit(survival::Surv(time, status) ~ 1)
  km_tail <- max(min(utils::tail(km$surv, 1), 0.99), 0.01)
  
  start_list <- list()
  for (i in 1:n_starts) {
    init_delta <- rep(0, ncol(X))
    init_delta[1] <- log(1 / mean(time)) + stats::rnorm(1, 0, 0.5)
    
    init_gamma <- rep(0, ncol(Z))
    if (type == "mixture") {
      init_gamma[1] <- stats::qlogis(1 - km_tail) # Susceptible proportion
    } else {
      init_gamma[1] <- log(-log(km_tail))
    }
    
    start_list[[i]] <- c(
      log_b = log(stats::runif(1, 0.5, 3)),
      log_c = log(stats::runif(1, 0.5, 3)),
      init_delta,
      init_gamma
    )
    names(start_list[[i]])[3:(2 + ncol(X))] <- paste0("delta_", colnames(X))
    names(start_list[[i]])[(3 + ncol(X)):length(start_list[[i]])] <- paste0("gamma_", colnames(Z))
  }
  
  # Optimize
  fit <- optim_multistart(
    ll_fun, start_list, method = method,
    accept = .bd_make_accept(n = length(time),
                             log_shape = c("log_b", "log_c"),
                             log_scale = character(0)))
  if (is.null(fit)) stop("Cure model optimization failed to converge.")
  
  est <- fit$estimate
  vcov_mat <- tryCatch(solve(-fit$hessian), error = function(e) matrix(NA, length(est), length(est)))
  rownames(vcov_mat) <- colnames(vcov_mat) <- names(est)
  
  out <- list(
    coefficients = est,
    logLik = fit$maximum,
    vcov = vcov_mat,
    convergence = fit$code,
    type = type,
    data = list(time = time, status = status, X = X, Z = Z),
    formula_aft = formula_aft,
    formula_cure = formula_cure,
    call = match.call(),
    submodel = TRUE
  )
  class(out) <- "bd_cure"
  return(out)
}
