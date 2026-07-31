#' Fit a Beta-Danish Competing Risks Model
#'
#' Fits a parametric competing risks model assuming independent latent failure
#' times, with each cause-specific baseline following the Beta-Danish
#' distribution. Covariates may be included, entering each cause through an
#' accelerated failure time link.
#'
#' @param time Numeric vector of observed times.
#' @param cause Integer vector of event causes: `0` for right-censored and
#'   `1, 2, ..., m` for the competing causes.
#' @param covariates Optional covariates. Either a one-sided formula such as
#'   `~ age + group`, evaluated in `data`, or a numeric matrix or data frame
#'   with one row per observation. `NULL` (default) fits cause-specific
#'   baselines with no regression structure.
#' @param data Data frame in which to evaluate `covariates` when it is a
#'   formula.
#' @param submodel Logical; if `TRUE`, fix `a = 1` in every cause, giving
#'   Exponentiated Danish cause-specific kernels. Default `FALSE`.
#' @param n_starts Integer; number of starting points for the joint
#'   optimisation.
#' @param method Optimisation method passed to [maxLik::maxLik()].
#'
#' @return An object of class `"bd_competing"`.
#'
#' @details
#' Writing \eqn{\delta_i} for the observed cause and \eqn{d_i} for the event
#' indicator, the log-likelihood is
#' \deqn{\ell = \sum_{i: d_i = 1}\Bigl[\log f_{\delta_i}(y_i)
#'        + \sum_{j \neq \delta_i}\log S_j(y_i)\Bigr]
#'        + \sum_{i: d_i = 0}\sum_{j=1}^{m}\log S_j(y_i).}
#' An event contributes the log-density of its own cause and the log-survival
#' of every other cause at the same time; a censored observation contributes
#' the log-survival of every cause.
#'
#' @section Covariates:
#' Each cause carries its own coefficient vector \eqn{\gamma_j}, acting on the
#' time scale:
#' \deqn{T_j = T_{0j}\exp(x'\gamma_j),}
#' so that a positive coefficient lengthens time to failure from that cause.
#' The likelihood contribution of an event at \eqn{t} is evaluated at the
#' accelerated time \eqn{t\exp(-x'\gamma_j)} with the Jacobian
#' \eqn{-x'\gamma_j} added on the log scale.
#'
#' The same design matrix is used for every cause, but the coefficients are
#' free to differ, so a covariate may accelerate one cause and retard another.
#' With `m` causes, `p` covariates and the full Beta-Danish kernel the model
#' carries `m * (4 + p)` parameters, which grows quickly: check
#' `summary()` for standard errors before interpreting any of them.
#'
#' @section Identifiability of the cause-specific marginals:
#' Mutual independence of the latent failure times is an *identifying*
#' assumption, not an empirically testable property. Tsiatis (1975) showed that
#' without it infinitely many joint distributions generate the same observable
#' law of \eqn{(T, \delta)}, so the marginals cannot be recovered from the
#' observed data alone.
#'
#' Under positive latent dependence the working independence model overstates
#' each cause-specific survival at moderate times, biasing the fitted
#' cumulative incidence functions downward; negative dependence reverses both.
#' The overall survival \eqn{S(t) = \prod_j S_j(t)} is more robust than the
#' individual marginals. Where conclusions rest on absolute cause-specific
#' CIFs rather than on model selection, supplement this fit with a
#' copula-based sensitivity analysis.
#'
#' @references
#' Tsiatis, A. (1975). A nonidentifiability aspect of the problem of competing
#' risks. *Proceedings of the National Academy of Sciences*, 72(1), 20-22.
#' \doi{10.1073/pnas.72.1.20}
#'
#' @seealso [cif_betadanish()], [cif_compare()], [simulate_bd_competing_data()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Without covariates
#' d <- simulate_bd_competing_data(150, seed = 1)
#' fit <- fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 1)
#' fit
#'
#' # With one binary covariate. Note that simulate_bd_competing_data() only
#' # creates the covariate column when 'gammas' is supplied.
#' dx <- simulate_bd_competing_data(150, gammas = c(0.5, -0.5), seed = 2)
#' fit2 <- fit_bd_competing(dx$time, dx$cause, covariates = ~ x, data = dx,
#'                          submodel = TRUE, n_starts = 1)
#' fit2$coefficients
#' }
fit_bd_competing <- function(time, cause, covariates = NULL, data = NULL,
                             submodel = FALSE, n_starts = 5, method = "BFGS") {

  time  <- as.numeric(time)
  cause <- as.integer(cause)
  if (length(time) != length(cause))
    stop("'time' and 'cause' must have the same length.", call. = FALSE)
  if (any(!is.finite(time)) || any(time <= 0))
    stop("All times must be finite and strictly positive.", call. = FALSE)

  causes <- sort(unique(cause[cause > 0]))
  m <- length(causes)
  if (m < 2L)
    stop("Competing risks needs at least two distinct event causes.",
         call. = FALSE)

  X <- .bd_cr_design(covariates, data, n = length(time))
  p <- if (is.null(X)) 0L else ncol(X)
  per <- (if (submodel) 3L else 4L) + p

  ## ---- names for one cause -------------------------------------------------
  base_nm <- if (submodel) c("log_b", "log_c", "log_k") else
    c("log_a", "log_b", "log_c", "log_k")
  gam_nm  <- if (p) paste0("gamma_", colnames(X)) else character(0)
  one_nm  <- c(base_nm, gam_nm)
  all_nm  <- unlist(lapply(causes, function(j) paste0("c", j, ":", one_nm)))

  ## ---- cause-wise starting values -----------------------------------------
  starts0 <- numeric(0)
  for (j in causes) {
    dat_j <- data.frame(time = time, status = as.integer(cause == j))
    fit_j <- tryCatch(
      suppressWarnings(
        fit_betadanish(survival::Surv(time, status) ~ 1, data = dat_j,
                       submodel = submodel, n_starts = 2,
                       check_identifiability = FALSE)),
      error = function(e) NULL)
    base <- if (is.null(fit_j)) {
      med <- stats::median(time)
      if (submodel) log(c(2, 1.5, 1 / med)) else log(c(1, 2, 1.5, 1 / med))
    } else {
      log(fit_j$coefficients)
    }
    starts0 <- c(starts0, base, rep(0, p))
  }
  names(starts0) <- all_nm

  ## ---- joint log-likelihood ------------------------------------------------
  ll_joint <- function(pars) {
    total <- 0
    for (idx in seq_along(causes)) {
      j   <- causes[idx]
      off <- (idx - 1L) * per
      th  <- pars[(off + 1L):(off + per)]

      a <- if (submodel) 1 else exp(th[1])
      o <- if (submodel) 0L else 1L
      b <- exp(th[o + 1]); cc <- exp(th[o + 2]); k <- exp(th[o + 3])

      if (p) {
        lacc <- -as.numeric(X %*% th[(o + 4):(o + 3 + p)])
      } else {
        lacc <- rep(0, length(time))
      }
      tt <- time * exp(lacc)

      is_j <- cause == j
      if (any(is_j))
        total <- total + sum(
          suppressWarnings(dbetadanish(tt[is_j], a, b, cc, k, log = TRUE)) +
            lacc[is_j])
      if (any(!is_j))
        total <- total + sum(suppressWarnings(
          pbetadanish(tt[!is_j], a, b, cc, k, lower.tail = FALSE, log.p = TRUE)))
    }
    if (!is.finite(total)) return(-1e10)
    total
  }

  ## ---- guarded multi-start -------------------------------------------------
  start_list <- list(starts0)
  for (i in seq_len(max(0L, n_starts - 1L))) {
    s <- starts0 + stats::rnorm(length(starts0), 0, 0.25)
    names(s) <- all_nm
    start_list[[length(start_list) + 1L]] <- s
  }

  accept <- .bd_make_accept(
    n = length(time),
    log_shape = grep("log_(a|b|c)$", all_nm, value = TRUE),
    log_scale = grep("log_k$", all_nm, value = TRUE))

  fit <- optim_multistart(ll_joint, start_list, method = method, accept = accept)
  if (is.null(fit))
    stop("No admissible optimum was found for the competing risks model. ",
         "Every start either failed or reached a degenerate ridge. Try ",
         "submodel = TRUE, or fewer covariates.", call. = FALSE)

  ## ---- format --------------------------------------------------------------
  est_log <- fit$estimate
  names(est_log) <- all_nm

  vcov_log <- tryCatch(solve(-fit$hessian),
                       error = function(e) matrix(NA_real_, length(est_log),
                                                  length(est_log)))
  dimnames(vcov_log) <- list(all_nm, all_nm)
  se_log <- sqrt(pmax(diag(vcov_log), 0))

  ## Shapes and scale back to the natural scale; regression coefficients are
  ## already on their own scale and are left alone.
  is_log <- grepl("log_", all_nm, fixed = TRUE)
  est_rep <- est_log
  est_rep[is_log] <- exp(est_log[is_log])
  names(est_rep) <- sub("log_", "", all_nm, fixed = TRUE)

  se_rep <- se_log
  se_rep[is_log] <- se_log[is_log] * est_rep[is_log]   # delta method
  names(se_rep) <- names(est_rep)

  par_matrix <- matrix(est_rep, nrow = m, byrow = TRUE,
                       dimnames = list(paste0("Cause_", causes),
                                       sub("^c[0-9]+:", "", names(est_rep))[
                                         seq_len(per)]))
  se_matrix <- matrix(se_rep, nrow = m, byrow = TRUE,
                      dimnames = dimnames(par_matrix))

  out <- list(coefficients = par_matrix,
              se           = se_matrix,
              estimate_log = est_log,
              vcov         = vcov_log,
              logLik       = fit$maximum,
              causes       = causes,
              m            = m,
              submodel     = submodel,
              covariates   = if (p) colnames(X) else character(0),
              X            = X,
              npar         = length(est_log),
              nobs         = length(time),
              convergence  = fit$code,
              starts_ok       = .bd_or(attr(fit, "bd_starts_ok"), NA_integer_),
              starts_rejected = .bd_or(attr(fit, "bd_starts_rejected"), NA_integer_),
              data         = list(time = time, cause = cause),
              call         = match.call())
  class(out) <- "bd_competing"

  if (!is.na(out$starts_rejected) && out$starts_rejected > 0)
    warning(sprintf(paste0("%d starting point(s) reached a degenerate ridge ",
                           "and were discarded."), out$starts_rejected),
            call. = FALSE)

  out
}

#' Build the Competing-Risks Design Matrix
#'
#' Accepts a one-sided formula plus data, a matrix, or a data frame, and
#' returns a numeric design matrix without an intercept. The intercept is
#' excluded because each cause already has its own scale parameter `k`, which
#' the intercept would be confounded with.
#'
#' @noRd
.bd_cr_design <- function(covariates, data, n) {
  if (is.null(covariates)) return(NULL)

  if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L)
      stop("'covariates' must be a one-sided formula, such as ~ age + group.",
           call. = FALSE)
    if (is.null(data))
      stop("'data' is required when 'covariates' is a formula.", call. = FALSE)
    mf <- stats::model.frame(covariates, data, na.action = stats::na.pass)
    X  <- stats::model.matrix(covariates, mf)
  } else {
    X <- as.matrix(covariates)
    if (is.null(colnames(X)))
      colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]

  if (nrow(X) != n)
    stop("The covariates have ", nrow(X), " row(s) but there are ", n,
         " observation(s).", call. = FALSE)
  if (anyNA(X))
    stop("The covariates contain missing values; remove those rows first.",
         call. = FALSE)
  if (!ncol(X))
    stop("No covariates remain after dropping the intercept.", call. = FALSE)

  storage.mode(X) <- "double"
  X
}

#' Cause-Specific Parameters at a Covariate Value
#'
#' Returns the four Beta-Danish parameters for one cause, with any accelerated
#' failure time effect folded into the scale. Scaling time by \eqn{\lambda}
#' maps \eqn{k} to \eqn{k/\lambda}, so a covariate acting as
#' \eqn{T_j = T_{0j}\exp(x'\gamma_j)} is equivalent to replacing \eqn{k_j} by
#' \eqn{k_j\exp(-x'\gamma_j)}.
#'
#' @noRd
.bd_cr_pars <- function(fit, cause_idx, x = NULL) {
  row <- which(fit$causes == cause_idx)
  pm  <- fit$coefficients[row, ]

  a <- if (isTRUE(fit$submodel)) 1 else as.numeric(pm[["a"]])
  b <- as.numeric(pm[["b"]])
  c <- as.numeric(pm[["c"]])
  k <- as.numeric(pm[["k"]])

  if (length(fit$covariates)) {
    gam <- as.numeric(pm[paste0("gamma_", fit$covariates)])
    xv  <- if (is.null(x)) rep(0, length(gam)) else as.numeric(x)
    if (length(xv) != length(gam))
      stop("'x' must supply one value per covariate: ",
           paste(fit$covariates, collapse = ", "), call. = FALSE)
    k <- k * exp(-sum(xv * gam))
  }
  list(a = a, b = b, c = c, k = k)
}

#' Cumulative Incidence Function
#'
#' Computes the cumulative incidence of one cause from a fitted Beta-Danish
#' competing risks model,
#' \eqn{F_j(t) = \int_0^t f_j(u)\prod_{l \neq j} S_l(u)\,du}.
#'
#' @param fit An object of class `"bd_competing"`.
#' @param tvec Numeric vector of times at which to evaluate the CIF.
#' @param cause_idx The cause to evaluate, matching one of the fitted causes.
#' @param x Optional covariate values at which to evaluate, one per covariate
#'   in the fit. Defaults to zero for every covariate, which is the reference
#'   subject. Ignored when the model has no covariates.
#'
#' @return A numeric vector of probabilities the same length as `tvec`.
#'
#' @details
#' The integrand is the cause-specific density multiplied by the survival of
#' every other cause, so the sum of all cause CIFs approaches the overall
#' failure probability rather than one, and each individual CIF is bounded by
#' it.
#'
#' Any covariate effect is folded into the scale parameter before integration,
#' since accelerating time by \eqn{\lambda} is equivalent to dividing \eqn{k}
#' by \eqn{\lambda}.
#'
#' @seealso [fit_bd_competing()], [cif_compare()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' d <- simulate_bd_competing_data(120, seed = 2)
#' fit <- fit_bd_competing(d$time, d$cause, submodel = TRUE, n_starts = 1)
#' cif_betadanish(fit, tvec = c(1, 5, 10), cause_idx = 1)
#' }
cif_betadanish <- function(fit, tvec, cause_idx, x = NULL) {
  if (!inherits(fit, "bd_competing"))
    stop("'fit' must be a bd_competing object.", call. = FALSE)
  if (!(cause_idx %in% fit$causes))
    stop("cause_idx ", cause_idx, " is not among the fitted causes: ",
         paste(fit$causes, collapse = ", "), call. = FALSE)

  pj  <- .bd_cr_pars(fit, cause_idx, x)
  oth <- setdiff(fit$causes, cause_idx)
  pl  <- lapply(oth, function(l) .bd_cr_pars(fit, l, x))

  integrand <- function(u) {
    v <- dbetadanish(u, pj$a, pj$b, pj$c, pj$k)
    for (q in pl)
      v <- v * pbetadanish(u, q$a, q$b, q$c, q$k, lower.tail = FALSE)
    v[!is.finite(v)] <- 0
    v
  }

  vapply(as.numeric(tvec), function(t_val) {
    if (!is.finite(t_val) || t_val <= 0) return(0)
    tryCatch(stats::integrate(integrand, 0, t_val, subdivisions = 2000L,
                              rel.tol = 1e-8, stop.on.error = FALSE)$value,
             error = function(e) NA_real_)
  }, numeric(1))
}
