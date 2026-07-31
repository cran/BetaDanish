## Interval estimation for Beta-Danish fits.
##
## The Wald interval on the natural scale is the wrong default for this family.
## The thesis pipeline records the problem directly: on the breaking-stress data
## the symmetric interval for the tail index b crosses zero, the profile
## interval has no finite upper bound, and 29 percent of bootstrap resamples
## fail to bound b above. The functions here follow that finding rather than
## papering over it.

#' Reconstruct the Log-Likelihood of a Fitted Model
#'
#' Rebuilds the objective from the fit object so that profiling does not need
#' the original call. Honours the submodel and grouped settings.
#'
#' @noRd
.bd_refit_loglik <- function(object) {
  time     <- object$data$time
  status   <- object$data$status
  submodel <- isTRUE(object$submodel)
  grouped  <- isTRUE(object$grouped)
  delta    <- object$delta

  function(nat) {
    a <- if (submodel) 1 else nat[["a"]]
    b <- nat[["b"]]; c <- nat[["c"]]; k <- nat[["k"]]
    if (!is.finite(a) || !is.finite(b) || !is.finite(c) || !is.finite(k) ||
        a <= 0 || b <= 0 || c <= 0 || k <= 0) return(-Inf)
    lp <- if (grouped) {
      suppressWarnings(.bd_log_cell(time, delta, a, b, c, k))
    } else {
      suppressWarnings(dbetadanish(time, a, b, c, k, log = TRUE))
    }
    ls <- suppressWarnings(
      pbetadanish(time, a, b, c, k, lower.tail = FALSE, log.p = TRUE))
    v <- sum(status * lp + (1 - status) * ls)
    if (is.finite(v)) v else -Inf
  }
}

#' Log-Scale Wald Confidence Intervals
#'
#' Symmetric Wald intervals on the log-parameter scale, exponentiated back.
#'
#' @param object A fitted `"betadanish"` object.
#' @param level Confidence level. Default 0.95.
#'
#' @return A matrix with one row per parameter and columns `estimate`,
#'   `se`, `lower` and `upper`.
#'
#' @details
#' All four parameters are strictly positive, so a symmetric interval on the
#' natural scale can and does cross zero when a parameter is weakly identified.
#' Building the interval on the log scale and exponentiating keeps it inside
#' the parameter space:
#' \eqn{\hat\theta \exp(\pm z\, \mathrm{SE}(\log\hat\theta))}.
#'
#' This is still a Wald interval, and it inherits the usual weakness: it
#' assumes the log-likelihood is approximately quadratic. Where it is not --
#' which is exactly where these intervals matter -- prefer [bd_profile_ci()].
#'
#' @seealso [bd_profile_ci()], [bd_identified_coef()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(150, a = 1, b = 3, c = 2, k = 0.5, seed = 1)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 3)
#' bd_wald_ci(fit)
#' }
bd_wald_ci <- function(object, level = 0.95) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  est <- object$coefficients
  se  <- sqrt(pmax(diag(object$vcov), 0))
  se_log <- se / est                       # delta method, natural -> log
  z <- stats::qnorm(1 - (1 - level) / 2)

  out <- cbind(estimate = as.numeric(est),
               se       = as.numeric(se),
               lower    = as.numeric(est * exp(-z * se_log)),
               upper    = as.numeric(est * exp(z * se_log)))
  rownames(out) <- names(est)
  attr(out, "level") <- level
  attr(out, "scale") <- "log-scale Wald, exponentiated"
  out
}

#' Profile Likelihood Confidence Interval
#'
#' Profiles one parameter, maximising over the others, and inverts the
#' likelihood ratio test to obtain an interval.
#'
#' @param object A fitted `"betadanish"` object.
#' @param parameter Name of the parameter to profile, one of `"a"`, `"b"`,
#'   `"c"` or `"k"`. Default `"b"`.
#' @param level Confidence level. Default 0.95.
#' @param grid Optional numeric vector of values to scan. If `NULL` (default) a
#'   log-spaced grid is built around the estimate.
#' @param n_grid Number of grid points when `grid` is `NULL`.
#' @param method Optimisation method for the nuisance parameters, passed to
#'   [stats::optim()].
#'
#' @return An object of class `"bd_profile"` with the interval, the grid, the
#'   profile log-likelihood, and the critical value.
#'
#' @details
#' The interval is the set of values \eqn{\theta_0} for which
#' \eqn{2\{\ell_{\max} - \ell_p(\theta_0)\} \le \chi^2_{1,\,\mathrm{level}}},
#' with \eqn{\ell_p} the profile log-likelihood.
#'
#' **An open upper bound is a result, not a failure.** For a weakly identified
#' tail index the profile can stay inside the critical region all the way to
#' the top of the grid, in which case `upper` is returned as `Inf` and the
#' honest report is a lower bound rather than a point estimate with an
#' interval. The thesis pipeline records precisely this on the breaking-stress
#' data, where the profile gives \eqn{b \ge 26.1} with no finite upper bound.
#'
#' Widening `grid` will not manufacture a bound; it only confirms the flatness.
#'
#' @seealso [bd_wald_ci()], [bd_identified_coef()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(200, a = 1, b = 3, c = 2, k = 0.5, seed = 2)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 3)
#' p <- bd_profile_ci(fit, "b")
#' p
#' }
bd_profile_ci <- function(object, parameter = "b", level = 0.95,
                          grid = NULL, n_grid = 60L, method = "Nelder-Mead") {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  est <- object$coefficients
  if (!parameter %in% names(est))
    stop("'", parameter, "' is not a parameter of this fit. Available: ",
         paste(names(est), collapse = ", "), call. = FALSE)

  free   <- setdiff(names(est), parameter)
  ll_at  <- .bd_refit_loglik(object)
  ll_max <- as.numeric(object$logLik)

  profile_at <- function(v) {
    if (!length(free)) {
      nat <- est; nat[parameter] <- v
      return(ll_at(nat))
    }
    obj <- function(lp) {
      nat <- est
      nat[free] <- exp(lp)
      nat[parameter] <- v
      val <- ll_at(nat)
      if (is.finite(val)) -val else 1e10
    }
    r <- tryCatch(stats::optim(log(est[free]), obj, method = method,
                               control = list(maxit = 3000, reltol = 1e-9)),
                  error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value)) -Inf else -r$value
  }

  if (is.null(grid)) {
    e <- as.numeric(est[[parameter]])
    grid <- sort(unique(exp(seq(log(e / 20), log(e * 20),
                                length.out = as.integer(n_grid)))))
  }
  grid <- sort(unique(grid[is.finite(grid) & grid > 0]))
  if (length(grid) < 3L)
    stop("'grid' needs at least three positive values.", call. = FALSE)

  prof <- vapply(grid, profile_at, numeric(1))
  crit <- stats::qchisq(level, df = 1) / 2
  dev  <- ll_max - prof

  inside <- grid[is.finite(dev) & dev <= crit]
  if (!length(inside)) {
    lower <- NA_real_; upper <- NA_real_
  } else {
    lower <- min(inside)
    hi    <- max(inside)
    tol   <- max(grid) * (1 - 1e-8)
    upper <- if (hi >= tol) Inf else hi
    if (lower <= min(grid) * (1 + 1e-8)) lower <- min(grid)
  }

  out <- list(parameter = parameter,
              estimate  = as.numeric(est[[parameter]]),
              lower = lower, upper = upper, level = level,
              grid = grid, profile = prof, logLik_max = ll_max,
              critical = crit,
              open_above = is.infinite(upper),
              open_below = length(inside) > 0 && lower <= min(grid) * (1 + 1e-8))
  class(out) <- "bd_profile"
  out
}

#' @param x A `"bd_profile"` object.
#' @param ... Ignored.
#' @rdname bd_profile_ci
#' @export
print.bd_profile <- function(x, ...) {
  cat("Profile likelihood interval for ", x$parameter, "\n", sep = "")
  cat("  estimate: ", signif(x$estimate, 6), "\n", sep = "")
  cat("  ", format(100 * x$level), "% interval: [",
      signif(x$lower, 6), ", ",
      if (is.infinite(x$upper)) "Inf" else signif(x$upper, 6), "]\n", sep = "")
  if (isTRUE(x$open_above))
    cat("\n  The profile stays within the critical region to the top of the\n",
        "  grid, so there is no finite upper bound. Report ", x$parameter,
        " as a\n  lower bound rather than as a point estimate.\n", sep = "")
  if (isTRUE(x$open_below))
    cat("\n  The interval also reaches the bottom of the grid; widen it if a\n",
        "  finite lower bound matters.\n", sep = "")
  invisible(x)
}

#' Coefficients in the Identified Parametrisation
#'
#' Reports the four-parameter fit in terms of the composite \eqn{ac}, which is
#' identified, rather than \eqn{a} and \eqn{c} separately, which are not.
#'
#' @param object A fitted `"betadanish"` object.
#' @param level Confidence level for the reported intervals. Default 0.95.
#'
#' @return An object of class `"bd_identified"` holding the reparametrised
#'   estimates with delta-method standard errors, and the condition numbers of
#'   the covariance matrix before and after.
#'
#' @details
#' Near the lower tail \eqn{a} and \eqn{c} enter the density almost entirely
#' through their product, so the expected information is close to singular in
#' the direction that separates them. Writing \eqn{p = ac} and \eqn{r = a/c},
#' the ratio \eqn{r} carries almost no information while \eqn{p} is estimated
#' precisely.
#'
#' Reporting \eqn{(ac, b, k)} is therefore the honest summary of a
#' four-parameter fit that sits on the flat direction: it has finite standard
#' errors and a far better conditioned covariance matrix, and it says exactly
#' what the data determine. The thesis pipeline records a condition number
#' falling from 2750 to 152 on the remission data, with the log-likelihood
#' unchanged, since this is a reparametrisation and not a different model.
#'
#' For the ED submodel \eqn{a = 1}, so \eqn{ac = c} and nothing is gained;
#' the function returns the coefficients unchanged with a note.
#'
#' @seealso [fit_betadanish()] and its Identifiability section,
#'   [bd_profile_ci()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(200, a = 1.5, b = 3, c = 2, k = 0.5, seed = 4)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       n_starts = 5)
#' bd_identified_coef(fit)
#' }
bd_identified_coef <- function(object, level = 0.95) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)

  est <- object$coefficients
  V   <- object$vcov
  z   <- stats::qnorm(1 - (1 - level) / 2)

  if (isTRUE(object$submodel) || !all(c("a", "c") %in% names(est))) {
    se <- sqrt(pmax(diag(V), 0))
    tab <- cbind(estimate = as.numeric(est), se = as.numeric(se),
                 lower = as.numeric(est * exp(-z * se / est)),
                 upper = as.numeric(est * exp(z * se / est)))
    rownames(tab) <- names(est)
    out <- list(table = tab, submodel = TRUE,
                condition_before = .bd_condition(V),
                condition_after  = .bd_condition(V), level = level)
    class(out) <- "bd_identified"
    return(out)
  }

  keep <- c("a", "b", "c", "k")
  V4 <- V[keep, keep, drop = FALSE]
  a  <- est[["a"]]; cc <- est[["c"]]

  ## (a, b, c, k) -> (ac, b, k)
  J <- rbind(c(cc, 0, a, 0),
             c(0,  1, 0, 0),
             c(0,  0, 0, 1))
  Vn <- J %*% V4 %*% t(J)
  nm <- c("ac", "b", "k")
  dimnames(Vn) <- list(nm, nm)

  e_new <- c(ac = a * cc, b = est[["b"]], k = est[["k"]])
  se    <- sqrt(pmax(diag(Vn), 0))

  tab <- cbind(estimate = as.numeric(e_new), se = as.numeric(se),
               lower = as.numeric(e_new * exp(-z * se / e_new)),
               upper = as.numeric(e_new * exp(z * se / e_new)))
  rownames(tab) <- nm

  out <- list(table = tab, submodel = FALSE,
              a = a, c = cc, ratio = a / cc,
              condition_before = .bd_condition(V4),
              condition_after  = .bd_condition(Vn),
              level = level)
  class(out) <- "bd_identified"
  out
}

#' Condition Number of a Covariance Matrix
#' @noRd
.bd_condition <- function(V) {
  ev <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (anyNA(ev) || min(ev) <= 0) return(NA_real_)
  max(ev) / min(ev)
}

#' @param x A `"bd_identified"` object.
#' @param ... Ignored.
#' @rdname bd_identified_coef
#' @export
print.bd_identified <- function(x, ...) {
  cat("Identified parametrisation\n")
  if (isTRUE(x$submodel)) {
    cat("  ED submodel: a is fixed at 1, so ac = c and nothing is gained.\n\n")
  } else {
    cat("  a and c enter mainly through their product; ac is reported.\n\n")
  }
  print(round(x$table, 5))
  if (!isTRUE(x$submodel)) {
    cat("\n  a = ", signif(x$a, 5), ", c = ", signif(x$c, 5),
        ", ratio a/c = ", signif(x$ratio, 5), "\n", sep = "")
  }
  cb <- x$condition_before; ca <- x$condition_after
  cat("  covariance condition number: ",
      if (is.finite(cb)) signif(cb, 5) else "NA", " -> ",
      if (is.finite(ca)) signif(ca, 5) else "NA", "\n", sep = "")
  invisible(x)
}
