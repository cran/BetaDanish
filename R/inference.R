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
  se  <- .bd_se_from_vcov(object$vcov, names(est))   # NA, never 0, if degenerate
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
#'   log-spaced grid is built around the estimate and extended as needed; a
#'   grid you supply is used exactly as given.
#' @param n_grid Number of points in the initial grid when `grid` is `NULL`.
#' @param method Optimisation method for the nuisance parameters, passed to
#'   [stats::optim()].
#' @param extend Logical; when `grid` is `NULL`, extend the grid outward until
#'   the profile crosses the cut-off or the parameter box is reached. Default
#'   `TRUE`.
#' @param max_extend Maximum number of extra grid points on each side.
#'   Default 150.
#'
#' @return An object of class `"bd_profile"` with the interval (`lower`,
#'   `upper`), the flags `lower_closed` and `upper_closed`, `open_below` and
#'   `open_above`, the grid, the profile log-likelihood, the critical value,
#'   `logLik_fit` (the log-likelihood of the fit) and `logLik_max` (the
#'   highest log-likelihood reached, from which the interval is measured).
#'
#' @details
#' The interval is the set of values \eqn{\theta_0} for which
#' \eqn{2\{\ell_{\max} - \ell_p(\theta_0)\} \le \chi^2_{1,\,\mathrm{level}}},
#' with \eqn{\ell_p} the profile log-likelihood.
#'
#' With the default grid the profile is first evaluated on a log-spaced grid
#' from \eqn{\hat\theta/20} to \eqn{20\hat\theta}, every point started at the
#' maximum-likelihood estimate. Where the profile has not yet fallen below the
#' cut-off at an end of that grid, the grid is extended outward in the same
#' logarithmic steps, each new point started from its neighbour's nuisance
#' estimates, until the profile drops below the cut-off or the parameter box
#' of [fit_betadanish()] is reached. A crossing that is bracketed by two grid
#' points is located by linear interpolation on the log scale, and the
#' endpoint is flagged closed. This follows the profile routine of the thesis,
#' which replaced a fixed grid whose edges had been reported as interval
#' endpoints.
#'
#' **An open upper bound is a result, not a failure.** For a weakly identified
#' tail index the profile can stay inside the critical region all the way to
#' the edge of the parameter box, in which case `upper` is returned as `Inf`
#' and the honest report is a lower bound rather than a point estimate with an
#' interval. An open lower end is returned as the lowest value scanned, with
#' `open_below = TRUE`.
#'
#' The interval is measured from the highest log-likelihood reached, by the fit
#' or by the profile. If the profile rises above the fitted maximum, the fit
#' stopped short of the optimum; a warning says so, except for a penalised fit,
#' whose estimate is by design not the unpenalised maximum.
#'
#' @references
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan.
#'
#' @seealso [bd_wald_ci()], [bd_identified_coef()], [bd_profile_plot()]
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
                          grid = NULL, n_grid = 60L, method = "Nelder-Mead",
                          extend = TRUE, max_extend = 150L) {
  if (!inherits(object, "betadanish"))
    stop("'object' must be a fitted betadanish model.", call. = FALSE)
  est <- object$coefficients
  if (!parameter %in% names(est))
    stop("'", parameter, "' is not a parameter of this fit. Available: ",
         paste(names(est), collapse = ", "), call. = FALSE)

  free   <- setdiff(names(est), parameter)
  ll_at  <- .bd_refit_loglik(object)
  ll_max <- as.numeric(object$logLik)

  ## Profile at one value, starting the nuisance parameters at `start_free`
  ## (natural scale). Returns the profile value and the nuisance optimum, so
  ## that an extension can continue from its neighbour.
  profile_at <- function(v, start_free) {
    if (!length(free)) {
      nat <- est; nat[parameter] <- v
      return(list(pl = ll_at(nat), par = start_free))
    }
    obj <- function(lp) {
      nat <- est
      nat[free] <- exp(lp)
      nat[parameter] <- v
      val <- ll_at(nat)
      if (is.finite(val)) -val else 1e10
    }
    r <- tryCatch(stats::optim(log(start_free), obj, method = method,
                               control = list(maxit = 3000, reltol = 1e-9)),
                  error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value) || r$value >= 1e10)
      return(list(pl = -Inf, par = start_free))
    list(pl = -r$value, par = stats::setNames(exp(r$par), free))
  }

  auto <- is.null(grid)
  if (auto) {
    e <- as.numeric(est[[parameter]])
    grid <- exp(seq(log(e / 20), log(e * 20), length.out = as.integer(n_grid)))
  }
  grid <- sort(unique(grid[is.finite(grid) & grid > 0]))
  if (length(grid) < 3L)
    stop("'grid' needs at least three positive values.", call. = FALSE)

  start0 <- est[free]
  res  <- lapply(grid, function(v) profile_at(v, start0))
  prof <- vapply(res, function(z) z$pl, numeric(1))
  crit <- stats::qchisq(level, df = 1) / 2

  ## The interval is measured from the highest log-likelihood reached. A
  ## profile above the fitted maximum means the fit stopped short of it (or,
  ## for a penalised fit, that its estimate is not the unpenalised maximum).
  ll_fit <- ll_max
  top    <- suppressWarnings(max(prof[is.finite(prof)]))
  if (is.finite(top) && top > ll_max) ll_max <- top
  cut  <- ll_max - crit

  ## ---- extension by continuation (thesis .r2_profile) ------------------------
  n_ext <- 0L
  if (auto && isTRUE(extend)) {
    bx   <- .bd_param_box(object$data$time, .bd_or(object$box, c(1e-4, 1e4)))
    lo_b <- log(bx$lower[[parameter]])
    hi_b <- log(bx$upper[[parameter]])
    lg   <- log(grid)
    step <- lg[2L] - lg[1L]

    par_dn <- res[[1L]]$par; kk <- 0L
    while (is.finite(prof[1L]) && prof[1L] >= cut &&
           lg[1L] - step >= lo_b && kk < max_extend) {
      r <- profile_at(exp(lg[1L] - step), par_dn)
      lg <- c(lg[1L] - step, lg); prof <- c(r$pl, prof)
      par_dn <- r$par; kk <- kk + 1L
    }
    n_ext <- kk

    par_up <- res[[length(res)]]$par; kk <- 0L
    while (is.finite(prof[length(prof)]) && prof[length(prof)] >= cut &&
           lg[length(lg)] + step <= hi_b && kk < max_extend) {
      r <- profile_at(exp(lg[length(lg)] + step), par_up)
      lg <- c(lg, lg[length(lg)] + step); prof <- c(prof, r$pl)
      par_up <- r$par; kk <- kk + 1L
    }
    n_ext <- n_ext + kk
    grid <- exp(lg)
  }

  top <- suppressWarnings(max(prof[is.finite(prof)]))
  if (is.finite(top) && top > ll_max) { ll_max <- top; cut <- ll_max - crit }
  if (ll_max > ll_fit + 1e-6 && !isTRUE(object$penalty > 0))
    warning(sprintf(paste0("The profile reached a log-likelihood %.6g above ",
                           "the fitted maximum, so the fit stopped short of ",
                           "the optimum. The interval is measured from the ",
                           "higher value; refit with more starts."),
                    ll_max - ll_fit), call. = FALSE)

  m   <- length(grid)
  dev <- ll_max - prof
  ins <- which(is.finite(dev) & dev <= crit)

  lower_closed <- FALSE; upper_closed <- FALSE
  if (!length(ins)) {
    lower <- NA_real_; upper <- NA_real_
  } else {
    i1 <- min(ins); i2 <- max(ins)
    lower_closed <- i1 > 1L && is.finite(dev[i1 - 1L]) && dev[i1 - 1L] > crit
    upper_closed <- i2 < m  && is.finite(dev[i2 + 1L]) && dev[i2 + 1L] > crit
    if (auto) {
      lgr <- log(grid)
      lower <- if (lower_closed) {
        exp(lgr[i1 - 1L] + (cut - prof[i1 - 1L]) / (prof[i1] - prof[i1 - 1L]) *
              (lgr[i1] - lgr[i1 - 1L]))
      } else grid[i1]
      upper <- if (upper_closed) {
        exp(lgr[i2 + 1L] + (cut - prof[i2 + 1L]) / (prof[i2] - prof[i2 + 1L]) *
              (lgr[i2] - lgr[i2 + 1L]))
      } else Inf
    } else {
      ## A supplied grid is read exactly as in earlier versions.
      lower <- grid[i1]
      hi    <- grid[i2]
      tol   <- max(grid) * (1 - 1e-8)
      upper <- if (hi >= tol) Inf else hi
      if (lower <= min(grid) * (1 + 1e-8)) lower <- min(grid)
    }
  }

  out <- list(parameter = parameter,
              estimate  = as.numeric(est[[parameter]]),
              lower = lower, upper = upper, level = level,
              grid = grid, profile = prof, logLik_max = ll_max,
              logLik_fit = ll_fit,
              critical = crit,
              lower_closed = lower_closed,
              upper_closed = upper_closed,
              open_above = is.infinite(upper),
              open_below = length(ins) > 0 && !lower_closed &&
                           lower <= min(grid) * (1 + 1e-8),
              extended = auto && isTRUE(extend),
              n_extended = n_ext)
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
  if (!is.null(x$lower_closed) && is.finite(x$lower))
    cat("  endpoints: lower ", if (isTRUE(x$lower_closed)) "closed" else "not closed",
        ", upper ", if (isTRUE(x$upper_closed)) "closed" else "not closed",
        "\n", sep = "")
  if (isTRUE(x$open_above))
    cat("\n  The profile stays within the critical region to the top of the\n",
        "  grid, so there is no finite upper bound. Report ", x$parameter,
        " as a\n  lower bound rather than as a point estimate.\n", sep = "")
  if (isTRUE(x$open_below)) {
    if (isTRUE(x$extended))
      cat("\n  The profile stays within the critical region down to the edge\n",
          "  of the parameter box; the lower limit shown is that edge, not a\n",
          "  crossing.\n", sep = "")
    else
      cat("\n  The interval also reaches the bottom of the grid; widen it if a\n",
          "  finite lower bound matters.\n", sep = "")
  }
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
    se <- .bd_se_from_vcov(V, names(est))
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
  se    <- .bd_se_from_vcov(Vn, nm)

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
