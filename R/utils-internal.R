#' Robust Multi-Start Optimization using maxLik
#'
#' @param ll_fun Log-likelihood function to maximise.
#' @param start_list A list of named numeric vectors of starting values.
#' @param method Optimisation method (default `"BFGS"`).
#' @param accept Optional predicate taking a fitted `maxLik` object and
#'   returning `TRUE` if the optimum is admissible. Optima for which it returns
#'   `FALSE` are discarded rather than competing for the maximum. `NULL`
#'   (default) accepts anything finite, which is the historical behaviour.
#' @param grad Optional analytic gradient of `ll_fun`, passed to
#'   `maxLik::maxLik()`. `NULL` (default) uses numerical derivatives.
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
optim_multistart <- function(ll_fun, start_list, method = "BFGS", accept = NULL,
                             grad = NULL) {
  best_fit <- NULL
  best_ll  <- -Inf
  n_ok     <- 0L
  n_rej    <- 0L
  lls      <- numeric(0)

  for (start_vals in start_list) {
    fit <- tryCatch(
      maxLik::maxLik(logLik = ll_fun, grad = grad, start = start_vals,
                     method = method),
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

#' Analytic Score of the Complete-Sample Log-Likelihood
#'
#' The score vector of thesis Appendix B.1, returned on the log-parameter
#' scale (each component multiplied by its parameter, by the chain rule). It
#' holds for complete, ungrouped samples only; censored and grouped
#' likelihoods are differentiated numerically by the caller.
#'
#' With \eqn{w = kx/(1+kx)}, \eqn{u = w^c} and \eqn{r = u/(1-u)}:
#' \deqn{U_a = n\{\psi(a+b) - \psi(a)\} + c\sum\log w,}
#' \deqn{U_b = n\{\psi(a+b) - \psi(b)\} + \sum\log(1-u),}
#' \deqn{U_c = n/c + a\sum\log w - (b-1)\sum r\log w,}
#' \deqn{U_k = \{nac - (ac+1)\sum w - c(b-1)\sum r(1-w)\}/k.}
#' \eqn{\log w} is formed as `-log1p(1/(k*x))`, which keeps full relative
#' precision for every \eqn{kx > 0}.
#'
#' @param x Complete sample (positive, finite).
#' @param th Natural-scale parameters: `c(a, b, c, k)`, or `c(b, c, k)` when
#'   `submodel = TRUE`.
#' @param submodel Logical; `TRUE` for the ED submodel (\eqn{a = 1}).
#' @return Score on the log scale, or zeros where it cannot be evaluated, so
#'   that a gradient-based optimiser sees a flat penalty region.
#' @noRd
.bd_score_eta <- function(x, th, submodel = FALSE) {
  np <- if (submodel) 3L else 4L
  if (length(th) != np || any(!is.finite(th)) || any(th <= 0))
    return(rep(0, np))
  if (submodel) {
    a <- 1; b <- th[1]; cc <- th[2]; k <- th[3]
  } else {
    a <- th[1]; b <- th[2]; cc <- th[3]; k <- th[4]
  }
  n   <- length(x)
  kt  <- k * x
  lw  <- -log1p(1 / kt)            # log w
  clw <- cc * lw                   # log u
  l1m <- log(-expm1(clw))          # log(1 - u)
  r   <- exp(clw - l1m)            # u / (1 - u)
  w   <- exp(lw)
  omw <- 1 / (1 + kt)              # 1 - w
  dab <- digamma(a + b)
  Ub <- n * (dab - digamma(b)) + sum(l1m)
  Uc <- n / cc + a * sum(lw) - (b - 1) * sum(r * lw)
  Uk <- (n * a * cc - (a * cc + 1) * sum(w) - cc * (b - 1) * sum(r * omw)) / k
  g <- if (submodel) {
    c(Ub, Uc, Uk) * c(b, cc, k)
  } else {
    Ua <- n * (dab - digamma(a)) + cc * sum(lw)
    c(Ua, Ub, Uc, Uk) * c(a, b, cc, k)
  }
  if (any(!is.finite(g))) rep(0, np) else as.numeric(g)
}

#' Starting Values of the Thesis Estimator
#'
#' Stage 1 of the thesis estimator fits the Exponentiated Danish submodel from
#' four starts, \eqn{(b, c) \in \{(1,1), (2,2), (0.5,3), (4,0.8)\}}, with
#' \eqn{k = 1/\mathrm{median}(x)}. Stage 2 fits the Beta-Danish from seven
#' starts placed along the \eqn{(a, c)} ridge through the Stage-1 optimum,
#' \eqn{(a, b, c, k) = (f, b_{ED}, c_{ED}/f, k_{ED})} for
#' \eqn{f \in \{1/30, 1/8, 1/3, 1, 3, 8, 30\}}. The start \eqn{f = 1} is the
#' ED optimum itself, so the Beta-Danish log-likelihood can never fall below
#' the ED one.
#'
#' @param k0 Data-driven rate, the reciprocal of the median time.
#' @param ed Named natural-scale ED optimum `c(b = , c = , k = )`.
#' @param lower,upper Named natural-scale bounds for `a`, `b`, `c`, `k`.
#' @return A list of named numeric vectors on the log scale.
#' @noRd
.bd_ed_starts <- function(k0, lower, upper) {
  pairs <- list(c(1, 1), c(2, 2), c(0.5, 3), c(4, 0.8))
  lim <- .bd_start_limits(lower[c("b", "c", "k")], upper[c("b", "c", "k")])
  lo <- lim$lo; hi <- lim$hi
  lapply(pairs, function(s) {
    v <- pmin(pmax(c(s[1], s[2], k0), lo), hi)
    stats::setNames(log(v), c("log_b", "log_c", "log_k"))
  })
}

#' Ridge Starts Through the ED Optimum (see .bd_ed_starts)
#' @noRd
.bd_ridge_starts <- function(ed, lower, upper,
                             factors = c(1/30, 1/8, 1/3, 1, 3, 8, 30)) {
  lim <- .bd_start_limits(lower[c("a", "b", "c", "k")], upper[c("a", "b", "c", "k")])
  lo <- lim$lo; hi <- lim$hi
  lapply(factors, function(f) {
    v <- c(f, ed[["b"]], ed[["c"]] / f, ed[["k"]])
    ## f = 1 is the ED optimum itself and is left untouched.
    if (f != 1) v <- pmin(pmax(v, lo), hi)
    stats::setNames(log(v), c("log_a", "log_b", "log_c", "log_k"))
  })
}

#' Where Starting Values May Sit
#'
#' Starts sit strictly inside the box, a factor of ten from each edge, as in
#' the thesis. For a box narrower than a factor of 100 that inner region is
#' empty, and the box itself is used instead.
#'
#' @noRd
.bd_start_limits <- function(lower, upper) {
  lo <- pmax(lower * 10, 1e-8)
  hi <- pmin(upper / 10, 1e8)
  bad <- lo > hi
  lo[bad] <- lower[bad]; hi[bad] <- upper[bad]
  list(lo = lo, hi = hi)
}

#' Parameter Box for the Univariate Fit
#'
#' The shape parameters `a`, `b` and `c` are confined to `box`, the thesis box
#' \eqn{[10^{-4}, 10^{4}]} by default. The rate `k` carries the unit of time,
#' so its bound is the same box widened to follow the data: it covers both
#' \eqn{[10^{-4}, 10^{4}]} and \eqn{k \cdot \mathrm{median}(t) \in} `box`.
#' Changing the unit of time therefore cannot push an otherwise interior fit
#' onto the boundary.
#'
#' @param time Observed times.
#' @param box Numeric vector `c(lower, upper)`.
#' @return A list with named numeric vectors `lower` and `upper`.
#' @noRd
.bd_param_box <- function(time, box = c(1e-4, 1e4)) {
  if (!is.numeric(box) || length(box) != 2L || any(!is.finite(box)) ||
      box[1] <= 0 || box[2] <= box[1])
    stop("'box' must be two positive numbers, c(lower, upper), with ",
         "lower < upper.", call. = FALSE)
  med <- stats::median(time[is.finite(time) & time > 0])
  if (!is.finite(med) || med <= 0) med <- 1
  lower <- c(a = box[1], b = box[1], c = box[1],
             k = min(box[1], box[1] / med))
  upper <- c(a = box[2], b = box[2], c = box[2],
             k = max(box[2], box[2] / med))
  list(lower = lower, upper = upper)
}

#' Standard Errors from a Covariance Matrix
#'
#' A variance that is not finite, is negative, or gives a standard error at or
#' below `1e-10` comes from a degenerate information matrix. Such a standard
#' error is reported as `NA`, never as `0`.
#'
#' @param V Covariance matrix.
#' @param nm Optional names for the result.
#' @return Numeric vector of standard errors.
#' @noRd
.bd_se_from_vcov <- function(V, nm = NULL) {
  d <- suppressWarnings(as.numeric(diag(as.matrix(V))))
  se <- rep(NA_real_, length(d))
  ok <- is.finite(d) & d > 0
  se[ok] <- sqrt(d[ok])
  se[!is.finite(se) | se <= 1e-10] <- NA_real_
  names(se) <- if (is.null(nm)) colnames(V) else nm
  se
}

#' Kaplan-Meier Steps
#'
#' Distinct event times and the Kaplan-Meier survivor estimate just after each
#' of them. Written out so that it can run on any monotone transform of the
#' times without refitting a `survfit` object.
#'
#' @noRd
.bd_km_steps <- function(time, status) {
  ut <- sort(unique(time[status == 1]))
  S  <- numeric(length(ut)); s <- 1
  for (i in seq_along(ut)) {
    nr <- sum(time >= ut[i]); dv <- sum(time == ut[i] & status == 1)
    if (nr > 0) s <- s * (1 - dv / nr)
    S[i] <- s
  }
  list(time = ut, surv = S)
}

#' Goodness of Fit on the Probability-Integral-Transform Scale
#'
#' The Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling statistics
#' with the empirical distribution function replaced by the Kaplan-Meier
#' estimator of the transformed times \eqn{u = F(t; \hat\theta)}, on which
#' the fitted distribution function is the identity:
#' \deqn{D = \sup_{0 \le u \le U} |\hat F(u) - u|,}
#' \deqn{W^2 = n\int_0^U (\hat F(u) - u)^2\,du,}
#' \deqn{A^2 = n\int_0^U \frac{(\hat F(u) - u)^2}{u(1-u)}\,du.}
#' The upper limit \eqn{U} is 1 when the largest observation is an event and
#' the transformed largest observation when it is censored. \eqn{\hat F} is
#' a step function, so both integrals are evaluated exactly, piece by piece.
#' Without censoring the three statistics are exactly the classical ones.
#'
#' Ported verbatim from `.gof_censored_pit()` of the thesis master script.
#'
#' @param time,status The sample (status 1 = event, 0 = right-censored).
#' @param Ffun Function of `t` returning the fitted distribution function.
#' @return Named numeric vector `c(D, W2, A2)`; `NA` when not computable.
#' @noRd
.bd_gof_pit <- function(time, status, Ffun) {
  na <- c(D = NA_real_, W2 = NA_real_, A2 = NA_real_)
  ok <- is.finite(time) & !is.na(status)
  time <- time[ok]; status <- as.integer(status[ok] > 0)
  n <- length(time)
  if (n < 2L || sum(status) < 1L) return(na)
  u <- tryCatch(as.numeric(Ffun(time)), error = function(e) rep(NA_real_, n))
  if (length(u) != n || any(!is.finite(u))) return(na)
  u  <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  km <- .bd_km_steps(u, status)
  uj <- km$time
  Fj <- 1 - km$surv
  U  <- max(u)
  m  <- length(uj)
  if (m > 0L && Fj[m] >= 1 - 1e-12) { Fj[m] <- 1; U <- 1 }
  left  <- c(0, uj)
  right <- c(uj, U)
  Fval  <- c(0, Fj)
  keep  <- right > left
  left <- left[keep]; right <- right[keep]; Fval <- Fval[keep]
  D  <- max(abs(Fval - left), abs(Fval - right))
  W2 <- n * sum(((Fval - left)^3 - (Fval - right)^3) / 3)
  G <- function(Fv, uu) {
    t1 <- ifelse(Fv == 0, 0, Fv^2 * log(uu))
    t2 <- ifelse(Fv == 1, 0, (1 - Fv)^2 * log1p(-uu))
    t1 - t2 - uu
  }
  lo_safe <- pmax(left, 1e-300)
  A2 <- n * sum(G(Fval, right) - ifelse(left <= 0, 0, G(Fval, lo_safe)))
  c(D = D, W2 = W2, A2 = A2)
}
