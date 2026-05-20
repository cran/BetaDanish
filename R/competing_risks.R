#' Fit Beta-Danish Competing Risks Model
#'
#' Fits a parametric competing risks model assuming independent latent failure
#' times, where each cause-specific baseline follows a 4-parameter Beta-Danish
#' distribution.
#'
#' @param time Numeric vector of observed times.
#' @param cause Integer vector of event causes. `0` indicates right-censored,
#'   and `1, 2, ..., m` indicate specific event causes.
#' @param n_starts Integer; number of random starts for the joint optimization.
#' @param method Optimization method passed to `maxLik`.
#'
#' @return An object of class `bd_competing`.
#'
#' @details
#' Under the assumption of independent latent failure times, the joint likelihood
#' factorizes. The function first fits independent Beta-Danish models for each
#' cause (treating other causes as censored) to find robust starting values,
#' then optimizes the joint likelihood.
#'
#' @export
fit_bd_competing <- function(time, cause, n_starts = 5, method = "BFGS") {

  if (length(time) != length(cause)) stop("time and cause must have the same length.")

  causes <- sort(unique(cause[cause > 0]))
  m <- length(causes)
  if (m < 2) stop("Competing risks analysis requires at least 2 distinct event causes.")

  # 1. Cause-wise initialization
  init_list <- list()
  for (j in causes) {
    # Treat cause j as event (1), everything else as censored (0)
    status_j <- ifelse(cause == j, 1, 0)

    # Fit single-cause model to get starting values
    # We use a small number of starts here to save time
    dat_j <- data.frame(time = time, status = status_j)
    fit_j <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat_j, n_starts = 3)

    init_list[[paste0("cause", j)]] <- log(fit_j$coefficients)
  }

  joint_start <- unlist(init_list)

  # 2. Joint Log-Likelihood
  ll_joint <- function(pars) {
    ll_total <- 0

    for (idx in seq_along(causes)) {
      j <- causes[idx]

      # Extract parameters for cause j
      p_idx <- ((idx - 1) * 4 + 1):(idx * 4)
      p_j <- exp(pars[p_idx])
      a <- p_j[1]; b <- p_j[2]; c <- p_j[3]; k <- p_j[4]

      # Times for events of cause j
      tj <- time[cause == j]
      # Times for everything else (censored for cause j)
      tc <- time[cause != j]

      if (length(tj) > 0) {
        log_f <- dbetadanish(tj, a, b, c, k, log = TRUE)
        ll_total <- ll_total + sum(log_f)
      }
      if (length(tc) > 0) {
        log_S <- pbetadanish(tc, a, b, c, k, lower.tail = FALSE, log.p = TRUE)
        ll_total <- ll_total + sum(log_S)
      }
    }

    if (!is.finite(ll_total)) return(-1e10)
    return(ll_total)
  }

  # 3. Joint Optimization with jittered starts
  start_list_joint <- list(joint_start)
  if (n_starts > 1) {
    for (i in 2:n_starts) {
      start_list_joint[[i]] <- joint_start + stats::rnorm(length(joint_start), 0, 0.2)
    }
  }

  fit <- optim_multistart(ll_joint, start_list_joint, method = method)
  if (is.null(fit)) stop("Competing risks optimization failed to converge.")

  # 4. Format Output
  est_log <- fit$estimate
  est_nat <- exp(est_log)

  # Reshape into a matrix
  par_matrix <- matrix(est_nat, nrow = m, ncol = 4, byrow = TRUE)
  colnames(par_matrix) <- c("a", "b", "c", "k")
  rownames(par_matrix) <- paste0("Cause_", causes)

  # Vcov via Delta method
  vcov_log <- tryCatch(solve(-fit$hessian), error = function(e) matrix(NA, length(est_log), length(est_log)))
  J <- diag(est_nat, nrow = length(est_nat))
  vcov_nat <- J %*% vcov_log %*% t(J)

  se_matrix <- matrix(sqrt(pmax(diag(vcov_nat), 0)), nrow = m, ncol = 4, byrow = TRUE)
  colnames(se_matrix) <- c("SE_a", "SE_b", "SE_c", "SE_k")
  rownames(se_matrix) <- paste0("Cause_", causes)

  out <- list(
    coefficients = par_matrix,
    se = se_matrix,
    logLik = fit$maximum,
    vcov = vcov_nat,
    convergence = fit$code,
    causes = causes,
    data = list(time = time, cause = cause),
    call = match.call(),
    submodel = FALSE
  )
  class(out) <- "bd_competing"
  return(out)
}

#' Compute Cumulative Incidence Function (CIF)
#'
#' Computes the CIF for a specific cause from a fitted Beta-Danish competing
#' risks model using numerical integration.
#'
#' @param fit An object of class `bd_competing`.
#' @param tvec Numeric vector of times at which to evaluate the CIF.
#' @param cause_idx Integer; the specific cause to evaluate (must match one of
#'   the causes in the fitted model).
#'
#' @return A numeric vector of CIF probabilities corresponding to `tvec`.
#' @export
cif_betadanish <- function(fit, tvec, cause_idx) {
  if (!inherits(fit, "bd_competing")) stop("fit must be a bd_competing object.")
  if (!(cause_idx %in% fit$causes)) stop("cause_idx not found in fitted model.")

  m <- length(fit$causes)
  pm <- fit$coefficients

  # Map cause_idx to row index
  row_idx <- which(fit$causes == cause_idx)

  res <- numeric(length(tvec))

  for (i in seq_along(tvec)) {
    t_val <- tvec[i]
    if (t_val <= 0) {
      res[i] <- 0
      next
    }

    # Integrand: h_j(u) * S_overall(u) = f_j(u) * prod_{l != j} S_l(u)
    integrand <- function(u) {
      # Density of cause of interest
      p_j <- pm[row_idx, ]
      f_j <- dbetadanish(u, p_j["a"], p_j["b"], p_j["c"], p_j["k"], log = FALSE)

      # Survival of other causes
      S_other <- rep(1, length(u))
      for (l in 1:m) {
        if (l != row_idx) {
          p_l <- pm[l, ]
          S_other <- S_other * pbetadanish(u, p_l["a"], p_l["b"], p_l["c"], p_l["k"], lower.tail = FALSE)
        }
      }
      return(f_j * S_other)
    }

    # Integrate from 0 to t_val
    int_res <- tryCatch({
      stats::integrate(integrand, lower = 0, upper = t_val, subdivisions = 2000, rel.tol = 1e-6)$value
    }, error = function(e) NA_real_)

    res[i] <- int_res
  }

  return(res)
}
