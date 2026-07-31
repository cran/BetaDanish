#' Robust Multi-Start Optimization using maxLik
#'
#' @param ll_fun Log-likelihood function to maximise.
#' @param start_list A list of named numeric vectors of starting values.
#' @param method Optimisation method (default `"BFGS"`).
#' @param accept Optional predicate taking a fitted `maxLik` object and
#'   returning `TRUE` if the optimum is admissible. Optima for which it returns
#'   `FALSE` are discarded rather than competing for the maximum. `NULL`
#'   (default) accepts anything finite, which is the historical behaviour.
#'
#' @details
#' Selecting purely on the highest log-likelihood is unsafe for this family.
#' The Beta-Danish likelihood has degenerate ridges along which the objective
#' can be driven arbitrarily high without a finite maximiser existing, so an
#' unguarded search will prefer a runaway to every honest optimum. When
#' `accept` is supplied, the best *admissible* optimum is returned instead.
#'
#' The returned object carries three attributes for diagnostics:
#' `bd_starts_ok`, `bd_starts_rejected` and `bd_loglik_spread`, the last being
#' the range of the accepted log-likelihoods.
#'
#' @return The best admissible `maxLik` object, or `NULL` if none qualified.
#' @noRd
optim_multistart <- function(ll_fun, start_list, method = "BFGS", accept = NULL) {
  best_fit <- NULL
  best_ll  <- -Inf
  n_ok     <- 0L
  n_rej    <- 0L
  lls      <- numeric(0)

  for (start_vals in start_list) {
    fit <- tryCatch(
      maxLik::maxLik(logLik = ll_fun, start = start_vals, method = method),
      error = function(e) NULL)

    if (is.null(fit) || is.null(fit$maximum) ||
        is.na(fit$maximum) || !is.finite(fit$maximum)) next

    if (!is.null(accept) && !isTRUE(accept(fit))) {
      n_rej <- n_rej + 1L
      next
    }

    n_ok <- n_ok + 1L
    lls  <- c(lls, fit$maximum)
    if (fit$maximum > best_ll) {
      best_ll  <- fit$maximum
      best_fit <- fit
    }
  }

  if (!is.null(best_fit)) {
    attr(best_fit, "bd_starts_ok")       <- n_ok
    attr(best_fit, "bd_starts_rejected") <- n_rej
    attr(best_fit, "bd_loglik_spread")   <-
      if (length(lls) > 1L) diff(range(lls)) else 0
  }
  best_fit
}

#' Admissibility Predicate for a Beta-Danish Optimum
#'
#' Ports the policy of `.ch6_degenerate()` in the thesis master script. That
#' function documents the mechanism: the density contains
#' \eqn{w = kx/(1 + kx)} raised to the power \eqn{c}, and as \eqn{k \to \infty}
#' and \eqn{c \to \infty} with \eqn{c/k} fixed, \eqn{w \to 1} and
#' \eqn{w^{c} \to \exp\{-(c/k)/x\}}: the family collapses to a different,
#' Frechet-type limiting distribution. That limit is a genuine ridge in the
#' likelihood, and an optimiser started in the wrong place slides down it.
#'
#' The master script records a real instance: an ED AFT fit returning
#' `b = 0.82`, `c = 296209.74` with log-likelihood -934.68, against a true
#' optimum of -907.19. Its rule, adopted verbatim here, is that any fit whose
#' shape parameters explode past `lim_shape` is rejected however good its
#' log-likelihood looks.
#'
#' The per-observation log-likelihood bound is additional. The master script
#' does not have it, but a four-parameter fit on the transplant data returned
#' a log-likelihood of `+5.2e77` on 91 observations, which no shape bound alone
#' would have caught.
#'
#' @param n Number of observations, used for the per-observation bound. Pass 0
#'   to disable that check.
#' @param log_shape,log_scale Names of the shape and scale parameters on the
#'   optimisation (log) scale.
#' @param lim_shape Shape parameters above this are degenerate. Default 500,
#'   the master script's value.
#' @param lim_scale Scale parameters outside `[1/lim_scale, lim_scale]` are
#'   degenerate. Default 1e8, the master script's `big`.
#' @param lim_ll_per_obs Mean log-density above this is implausible for real
#'   data and indicates a runaway rather than a fit.
#'
#' @return A function of one argument, suitable as `accept` in
#'   [optim_multistart()].
#' @noRd
.bd_make_accept <- function(n,
                            log_shape = c("log_a", "log_b", "log_c"),
                            log_scale = "log_k",
                            lim_shape = 500,
                            lim_scale = 1e8,
                            lim_ll_per_obs = 50) {
  force(n)
  function(fit) {
    if (is.null(fit) || is.null(fit$maximum) || !is.finite(fit$maximum))
      return(FALSE)
    if (n > 0 && fit$maximum / n > lim_ll_per_obs) return(FALSE)

    est <- fit$estimate
    if (is.null(est) || any(!is.finite(est))) return(FALSE)

    sh <- est[names(est) %in% log_shape]
    if (length(sh) && any(exp(sh) > lim_shape)) return(FALSE)

    sc <- est[names(est) %in% log_scale]
    if (length(sc) && (any(exp(sc) > lim_scale) ||
                       any(exp(sc) < 1 / lim_scale))) return(FALSE)

    TRUE
  }
}

#' Deterministic Starting Grid
#'
#' The shape combinations are those of `default_starts()` in the thesis master
#' script, which were chosen against the datasets in the dissertation. The
#' scale is set from the data rather than fixed, so the grid transfers to
#' arbitrary time units.
#'
#' A deterministic grid matters more than the number of random starts. The
#' master script uses fixed, sanely placed starts and rarely runs away; the
#' package used `runif(0.5, 5)` on every shape and did.
#'
#' @param submodel Logical; `TRUE` for the three-parameter ED submodel.
#' @param k_base Data-driven scale, typically the reciprocal of the mean
#'   observed event time.
#' @return A list of named numeric vectors on the log scale.
#' @noRd
.bd_default_starts <- function(submodel, k_base) {
  shapes <- if (submodel) {
    list(c(b = 5.4, c = 2.3), c(b = 5.0, c = 2.3), c(b = 4.5, c = 2.3),
         c(b = 2.0, c = 2.0), c(b = 4.0, c = 1.5), c(b = 1.5, c = 2.5),
         c(b = 8.0, c = 2.0), c(b = 3.0, c = 1.0))
  } else {
    list(c(a = 3.0, b = 4.0,  c = 2.3), c(a = 3.0, b = 24.0, c = 2.3),
         c(a = 5.7, b = 2.5,  c = 2.3), c(a = 1.0, b = 2.0,  c = 2.0),
         c(a = 2.0, b = 2.0,  c = 1.0), c(a = 0.7, b = 4.0,  c = 1.5),
         c(a = 3.0, b = 5.0,  c = 2.0), c(a = 1.2, b = 10.0, c = 3.0),
         c(a = 0.5, b = 20.0, c = 2.0))
  }

  out <- list()
  for (kf in c(1, 0.25)) {
    for (s in shapes) {
      v <- c(s, k = k_base * kf)
      out[[length(out) + 1L]] <- stats::setNames(log(v), paste0("log_", names(v)))
    }
  }
  out
}

#' Log-Space Addition
#' Computes log(exp(x) + exp(y)) stably
#' @noRd
logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}
