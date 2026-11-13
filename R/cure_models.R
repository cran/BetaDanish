#' Fit Beta-Danish Cure Models
#'
#' Fits mixture and promotion-time (non-mixture) cure models with a
#' Beta-Danish AFT latency distribution: the three-parameter Exponentiated
#' Danish kernel by default, or the four-parameter Beta-Danish kernel with
#' `submodel = FALSE`.
#'
#' @param formula_aft A formula for the latency component (e.g., `Surv(time, status) ~ age`).
#' @param formula_cure A one-sided formula for the incidence/cure component (e.g., `~ treatment`).
#' @param data A data frame containing the variables.
#' @param type Character; either `"mixture"` or `"promotion"`.
#' @param n_starts Integer; number of random starts for optimization.
#' @param method Optimization method passed to `maxLik`.
#' @param submodel Logical; `TRUE` (default) uses the Exponentiated Danish
#'   kernel (\eqn{a = 1}), as in earlier versions. `FALSE` fits the
#'   four-parameter Beta-Danish kernel, which adds `log_a` to the parameter
#'   vector.
#'
#' @return An object of class `bd_cure` with components including
#'   `coefficients`, `logLik`, `vcov`, `npar`, `nobs`, `AIC`, `BIC`,
#'   `separated` (logical; see Details), `submodel` and, for the
#'   four-parameter kernel, `logLik_submodel`, the log-likelihood of the
#'   Exponentiated Danish cure fit it is nested over.
#'
#' @details
#' In the **mixture** model, the population is split into susceptible and cured
#' fractions. The susceptible probability is modeled via logistic regression:
#' `pi = exp(Z \%*\% gamma) / (1 + exp(Z \%*\% gamma))`.
#'
#' In the **promotion-time** (non-mixture) model, the cure fraction is derived
#' from a latent Poisson process of clonogenic cells: `theta = exp(Z \%*\% gamma)`.
#' The cure fraction is `exp(-theta)`.
#'
#' @section Four-parameter kernel:
#' With `submodel = FALSE` the Exponentiated Danish cure model is fitted
#' first. Its optimum, with \eqn{\log a = 0} inserted, is a feasible point of
#' the Beta-Danish model and is used as a start, together with starts along
#' the \eqn{(a, c)} ridge through it (\eqn{a = f}, \eqn{c = c_{ED}/f} for
#' \eqn{f \in \{1/30, 1/8, 1/3, 1, 3, 8, 30\}}) and random starts. A
#' Beta-Danish fit with a lower log-likelihood than its own Exponentiated
#' Danish submodel can only mean that the optimiser failed, so it is refused
#' with an error rather than reported. Information criteria count every
#' element of the parameter vector: the baseline rate is carried by the
#' latency intercept, so nothing is redundant. This follows the Chapter 6
#' covariate cure ladders of the thesis.
#'
#' The latency rate is \eqn{k_i = \exp(x_i^\top\delta)} in both kernels. A fit
#' whose shape parameters exceed 500 is rejected as the degenerate limit of the
#' family. A fit with a regression coefficient (other than an intercept)
#' above 15 in absolute value is flagged `separated = TRUE`: it is a valid
#' maximum of the likelihood but weakly identified in that direction.
#'
#' @references
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan.
#'
#' @examples
#' \donttest{
#' set.seed(5)
#' n <- 200
#' grp <- rep(0:1, length.out = n)
#' cured <- stats::rbinom(n, 1, ifelse(grp == 1, 0.4, 0.2)) == 1
#' t_u <- rbetadanish(n, a = 1, b = 3, c = 1.5, k = 0.3)
#' cens <- stats::runif(n, 0, 30)
#' t_obs <- ifelse(cured, cens, pmin(t_u, cens))
#' dat <- data.frame(time = t_obs, status = as.integer(!cured & t_u <= cens),
#'                   group = grp)
#' fit <- fit_bd_cure(survival::Surv(time, status) ~ 1, ~ group,
#'                    data = dat, type = "mixture", n_starts = 3)
#' summary(fit)
#' }
#'
#' @export
fit_bd_cure <- function(formula_aft, formula_cure, data, type = c("mixture", "promotion"),
                        n_starts = 10, method = "BFGS", submodel = TRUE) {
  type <- match.arg(type)
  if (!is.logical(submodel) || length(submodel) != 1L || is.na(submodel))
    stop("'submodel' must be TRUE or FALSE.", call. = FALSE)

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

  # Log-likelihood. `bd = FALSE`: (log_b, log_c, delta, gamma), a = 1.
  # `bd = TRUE`: (log_a, log_b, log_c, delta, gamma); log_a = 0 reproduces the
  # ED likelihood exactly, which is what makes the nesting check meaningful.
  make_ll <- function(bd) {
    off <- if (bd) 1L else 0L
    function(pars) {
      a <- if (bd) exp(pars[1]) else 1
      b <- exp(pars[1 + off])
      c <- exp(pars[2 + off])
      delta <- pars[(3 + off):(2 + off + ncol(X))]
      gamma <- pars[(3 + off + ncol(X)):length(pars)]

      # Latency scale
      k_i <- exp(as.numeric(X %*% delta))

      # Baseline functions
      log_f <- dbetadanish(time, a = a, b = b, c = c, k = k_i, log = TRUE)
      log_S <- pbetadanish(time, a = a, b = b, c = c, k = k_i, lower.tail = FALSE, log.p = TRUE)
      F_u <- pbetadanish(time, a = a, b = b, c = c, k = k_i, lower.tail = TRUE, log.p = FALSE)

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
  }
  ll_fun <- make_ll(FALSE)

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

  # Optimize the Exponentiated Danish kernel (always: it is the model when
  # submodel = TRUE, and the nesting reference when submodel = FALSE)
  fit <- optim_multistart(
    ll_fun, start_list, method = method,
    accept = .bd_make_accept(n = length(time),
                             log_shape = c("log_b", "log_c"),
                             log_scale = character(0)))
  if (is.null(fit)) stop("Cure model optimization failed to converge.")
  ll_ed <- fit$maximum

  if (!submodel) {
    ed <- fit$estimate
    ridge <- lapply(c(1/30, 1/8, 1/3, 1, 3, 8, 30), function(f) {
      s <- c(log_a = log(f), ed)
      s[["log_c"]] <- ed[["log_c"]] - log(f)
      s
    })
    rand <- lapply(start_list, function(s)
      c(log_a = log(stats::runif(1, 0.5, 3)), s))
    fit_bd <- optim_multistart(
      make_ll(TRUE), c(ridge, rand), method = method,
      accept = .bd_make_accept(n = length(time),
                               log_shape = c("log_a", "log_b", "log_c"),
                               log_scale = character(0)))
    if (is.null(fit_bd))
      stop("The four-parameter cure model found no admissible optimum. ",
           "Use submodel = TRUE.", call. = FALSE)
    if (fit_bd$maximum < ll_ed - 1e-6)
      stop(sprintf(paste0("The four-parameter cure fit (logLik %.4f) is below ",
                          "its Exponentiated Danish submodel (logLik %.4f). ",
                          "The ED optimum is a feasible point of the larger ",
                          "model, so the optimiser has failed; the fit is not ",
                          "reported. Use submodel = TRUE."),
                   fit_bd$maximum, ll_ed), call. = FALSE)
    fit <- fit_bd
  }

  est <- fit$estimate
  vcov_mat <- tryCatch(solve(-fit$hessian), error = function(e) matrix(NA, length(est), length(est)))
  rownames(vcov_mat) <- colnames(vcov_mat) <- names(est)

  npar <- length(est)
  nobs <- length(time)
  reg  <- est[!grepl("^log_", names(est)) & !grepl("(Intercept)", names(est), fixed = TRUE)]

  out <- list(
    coefficients = est,
    logLik = fit$maximum,
    vcov = vcov_mat,
    convergence = fit$code,
    npar = npar,
    nobs = nobs,
    AIC = 2 * npar - 2 * fit$maximum,
    BIC = npar * log(nobs) - 2 * fit$maximum,
    logLik_submodel = if (submodel) NA_real_ else ll_ed,
    separated = length(reg) > 0 && any(abs(reg) > 15),
    type = type,
    data = list(time = time, status = status, X = X, Z = Z),
    formula_aft = formula_aft,
    formula_cure = formula_cure,
    call = match.call(),
    submodel = submodel
  )
  class(out) <- "bd_cure"
  return(out)
}
