## Moment-based structural properties of the Beta-Danish distribution.
##
## Every integral below is taken on the finite Beta(a, b) scale. Under
## u = G(z) = {kz/(1+kz)}^c the variable u is exactly Beta(a, b) distributed and
## z = v / (k(1 - v)) with v = u^(1/c), so
##
##   E[g(Z)] = 1/B(a,b) * INT_0^1 g(z(u)) u^(a-1) (1-u)^(b-1) du.
##
## The interval is finite, the only singularity is the integrable endpoint
## u -> 1, and there is no series truncation.

#' Back-Transform from the Beta Scale to the Time Scale
#' @noRd
.bd_z_of_u <- function(u, c, k) {
  v <- u^(1 / c)
  v / (k * (1 - v))
}

#' Core Integrand: z(u)^r Times the Beta(a, b) Density
#' @noRd
.bd_mom_integrand <- function(u, r, a, b, c, k, lB) {
  z  <- .bd_z_of_u(u, c, k)
  ## r = 0 must contribute nothing: r * log(z) would be 0 * -Inf = NaN at the
  ## lower endpoint, where the integrand is simply the Beta(a, b) density.
  rz <- if (r == 0) 0 else r * log(z)
  lg <- rz + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
  out <- exp(lg)
  out[!is.finite(out)] <- 0
  out
}

#' Integrate the Core Integrand Between Two Beta-Scale Limits
#' @noRd
.bd_mom_int <- function(r, a, b, c, k, lo = 0, hi = 1,
                        rel.tol = 1e-10, subdivisions = 4000L) {
  lB <- lbeta(a, b)
  tryCatch(
    stats::integrate(.bd_mom_integrand, lower = lo, upper = hi,
                     r = r, a = a, b = b, c = c, k = k, lB = lB,
                     rel.tol = rel.tol, subdivisions = subdivisions,
                     stop.on.error = FALSE)$value,
    error = function(e) NA_real_)
}

#' Raw Moments of the Beta-Danish Distribution
#'
#' @param r Order or orders of the moment. Non-negative.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail; `E(Z^r)` is finite iff `b > r`.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param subdivisions Subdivision limit, passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `r`. Entries with `b <= r` are `Inf`.
#'
#' @details
#' The survival function is regularly varying with index \eqn{-b}, so
#' \eqn{E(Z^r)} is finite if and only if \eqn{b > r}. That condition is checked
#' rather than assumed: a request for a moment that does not exist returns
#' `Inf`, not a large finite number produced by a truncated sum.
#'
#' @seealso [bd_moment_summary()], [bd_incomplete_moment()]
#'
#' @export
#'
#' @examples
#' bd_moments(1:2, a = 1.5, b = 3, c = 2, k = 1)
#'
#' # The fourth moment does not exist when b = 3
#' bd_moments(4, a = 1.5, b = 3, c = 2, k = 1)
bd_moments <- function(r, a, b, c, k, rel.tol = 1e-10, subdivisions = 4000L) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  r <- as.numeric(r)
  if (any(is.na(r)) || any(r < 0))
    stop("'r' must be non-negative.", call. = FALSE)

  vapply(r, function(ri) {
    if (b <= ri) return(Inf)
    .bd_mom_int(ri, a, b, c, k, 0, 1, rel.tol, subdivisions)
  }, numeric(1))
}

#' Summary Moments: Mean, Variance, Skewness and Kurtosis
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param subdivisions Subdivision limit, passed to [stats::integrate()].
#'
#' @return A named numeric vector. Entries requiring a moment that does not
#'   exist are `NA`, with the governing condition given in the `condition`
#'   attribute.
#'
#' @details
#' Skewness needs \eqn{b > 3} and kurtosis \eqn{b > 4}. Both are frequently
#' unavailable for fitted values of \eqn{b}, which is a property of the family
#' rather than a limitation of the computation.
#'
#' @export
#'
#' @examples
#' bd_moment_summary(a = 1.5, b = 5, c = 2, k = 1)
#' bd_moment_summary(a = 1.5, b = 2.5, c = 2, k = 1)   # skew/kurt unavailable
bd_moment_summary <- function(a, b, c, k, rel.tol = 1e-10,
                              subdivisions = 4000L) {
  m <- bd_moments(1:4, a, b, c, k, rel.tol = rel.tol,
                  subdivisions = subdivisions)
  mu <- m[1]

  out <- c(mean = NA_real_, variance = NA_real_,
           sd = NA_real_, skewness = NA_real_, kurtosis = NA_real_)

  if (is.finite(mu)) out[["mean"]] <- mu
  if (all(is.finite(m[1:2]))) {
    v <- m[2] - mu^2
    out[["variance"]] <- v
    out[["sd"]] <- sqrt(max(v, 0))
  }
  if (all(is.finite(m[1:3])) && is.finite(out[["sd"]]) && out[["sd"]] > 0) {
    mu3 <- m[3] - 3 * mu * m[2] + 2 * mu^3
    out[["skewness"]] <- mu3 / out[["sd"]]^3
  }
  if (all(is.finite(m[1:4])) && is.finite(out[["sd"]]) && out[["sd"]] > 0) {
    mu4 <- m[4] - 4 * mu * m[3] + 6 * mu^2 * m[2] - 3 * mu^4
    out[["kurtosis"]] <- mu4 / out[["sd"]]^4
  }
  attr(out, "condition") <- "E(Z^r) is finite if and only if b > r"
  out
}

#' Incomplete Moments
#'
#' The lower incomplete moment \eqn{M_r(t) = \int_0^t z^r f(z)\,dz} and the
#' upper incomplete moment \eqn{\Upsilon_r(t) = \int_t^\infty z^r f(z)\,dz}.
#'
#' @param r Order of the moment.
#' @param t Vector of thresholds.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param lower Logical; `TRUE` (default) for \eqn{M_r(t)}, `FALSE` for
#'   \eqn{\Upsilon_r(t)}.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `t`.
#'
#' @details
#' The threshold is mapped to the Beta scale by \eqn{u_t = G(t)}, computed on
#' the log scale so it stays accurate when \eqn{G(t)} approaches one, and the
#' integral is taken over \eqn{(0, u_t)} or \eqn{(u_t, 1)}.
#'
#' @seealso [bd_mrl()], [bd_lorenz()], [bd_mean_deviation()]
#'
#' @export
#'
#' @examples
#' a <- 1.5; b <- 3; c <- 2; k <- 1
#' t <- 1
#' # The two halves sum to the complete moment
#' bd_incomplete_moment(1, t, a, b, c, k, lower = TRUE) +
#'   bd_incomplete_moment(1, t, a, b, c, k, lower = FALSE)
#' bd_moments(1, a, b, c, k)
bd_incomplete_moment <- function(r, t, a, b, c, k, lower = TRUE,
                                 rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (length(r) != 1L || is.na(r) || r < 0)
    stop("'r' must be a single non-negative number.", call. = FALSE)
  if (b <= r)
    warning("The complete moment of order ", r, " does not exist (b <= r); ",
            "the upper incomplete moment is infinite.", call. = FALSE)

  vapply(as.numeric(t), function(ti) {
    if (is.na(ti)) return(NA_real_)
    if (ti <= 0) return(if (lower) 0 else bd_moments(r, a, b, c, k))
    if (!is.finite(ti)) return(if (lower) bd_moments(r, a, b, c, k) else 0)
    ut <- exp(.bd_logG(k * ti, c))
    if (lower) .bd_mom_int(r, a, b, c, k, 0, ut, rel.tol)
    else       .bd_mom_int(r, a, b, c, k, ut, 1, rel.tol)
  }, numeric(1))
}

#' Conditional Moments
#'
#' \eqn{E(Z^r \mid Z > t)} or \eqn{E(Z^r \mid Z \le t)}.
#'
#' @param r Order of the moment.
#' @param t Vector of thresholds.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param upper Logical; `TRUE` (default) conditions on \eqn{Z > t}.
#' @return A numeric vector the length of `t`.
#' @export
#' @examples
#' bd_conditional_moment(1, t = c(0.5, 1, 2), a = 1.5, b = 3, c = 2, k = 1)
bd_conditional_moment <- function(r, t, a, b, c, k, upper = TRUE,
                                  rel.tol = 1e-10) {
  num <- bd_incomplete_moment(r, t, a, b, c, k, lower = !upper, rel.tol = rel.tol)
  den <- if (upper) sbetadanish(t, a, b, c, k) else pbetadanish(t, a, b, c, k)
  out <- num / den
  out[!is.finite(den) | den <= 0] <- NA_real_
  out
}

#' Mean Residual Life and Reversed Mean Residual Life
#'
#' @param t Vector of times.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `t`, or all `NA` when \eqn{b \le 1}.
#'
#' @details
#' Following the definitions in the underlying dissertation, with
#' \eqn{\Upsilon_1} the upper and \eqn{M_1} the lower incomplete first moment,
#' \deqn{m(t) = E(Z - t \mid Z > t) = \Upsilon_1(t)/S(t) - t,}
#' \deqn{\bar m(t) = E(t - Z \mid Z \le t) = t - M_1(t)/F(t).}
#' The reversed form is also called the mean inactivity time. Both require the
#' first moment to exist, hence \eqn{b > 1}.
#'
#' No monotonicity is asserted. Ageing classification requires a separate
#' argument and is not implied by these values; see [bd_hazard_shape()] for the
#' hazard-rate counterpart.
#'
#' @export
#'
#' @examples
#' bd_mrl(c(0.5, 1, 2, 5), a = 1.5, b = 3, c = 2, k = 1)
#' bd_rmrl(c(0.5, 1, 2, 5), a = 1.5, b = 3, c = 2, k = 1)
bd_mrl <- function(t, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The mean residual life needs a finite first moment, so b > 1. ",
            "NA returned.", call. = FALSE)
    return(rep(NA_real_, length(t)))
  }
  ups <- bd_incomplete_moment(1, t, a, b, c, k, lower = FALSE, rel.tol = rel.tol)
  S   <- sbetadanish(t, a, b, c, k)
  out <- ups / S - as.numeric(t)
  out[!is.finite(S) | S <= 0] <- NA_real_
  out
}

#' @rdname bd_mrl
#' @export
bd_rmrl <- function(t, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The reversed mean residual life needs a finite first moment, so ",
            "b > 1. NA returned.", call. = FALSE)
    return(rep(NA_real_, length(t)))
  }
  M1 <- bd_incomplete_moment(1, t, a, b, c, k, lower = TRUE, rel.tol = rel.tol)
  Fx <- pbetadanish(t, a, b, c, k)
  out <- as.numeric(t) - M1 / Fx
  out[!is.finite(Fx) | Fx <= 0] <- NA_real_
  out
}

#' Mean Deviations
#'
#' Mean deviation about the mean, \eqn{E|Z - \mu|}, or about the median.
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param about `"mean"` (default) or `"median"`.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single number, or `NA` when \eqn{b \le 1}.
#'
#' @details
#' \eqn{E|Z - \mu| = 2\mu F(\mu) - 2 M_1(\mu)} and
#' \eqn{E|Z - M| = \mu - 2 M_1(M)}, with \eqn{M} the median.
#'
#' @export
#'
#' @examples
#' bd_mean_deviation(a = 1.5, b = 3, c = 2, k = 1)
#' bd_mean_deviation(a = 1.5, b = 3, c = 2, k = 1, about = "median")
bd_mean_deviation <- function(a, b, c, k, about = c("mean", "median"),
                              rel.tol = 1e-10) {
  about <- match.arg(about)
  if (b <= 1) {
    warning("Mean deviations need a finite first moment, so b > 1. ",
            "NA returned.", call. = FALSE)
    return(NA_real_)
  }
  mu <- bd_moments(1, a, b, c, k, rel.tol = rel.tol)
  if (identical(about, "mean")) {
    2 * mu * pbetadanish(mu, a, b, c, k) -
      2 * bd_incomplete_moment(1, mu, a, b, c, k, rel.tol = rel.tol)
  } else {
    med <- qbetadanish(0.5, a, b, c, k)
    mu - 2 * bd_incomplete_moment(1, med, a, b, c, k, rel.tol = rel.tol)
  }
}

#' Lorenz and Bonferroni Curves
#'
#' @param p Vector of probabilities in \eqn{[0, 1]}.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A numeric vector the length of `p`.
#'
#' @details
#' \eqn{L(p) = M_1(Q(p))/\mu} and \eqn{B(p) = L(p)/p}, with \eqn{Q} the
#' quantile function. Both require \eqn{b > 1}.
#'
#' @export
#'
#' @examples
#' p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' bd_lorenz(p, a = 1.5, b = 3, c = 2, k = 1)
#' bd_bonferroni(p, a = 1.5, b = 3, c = 2, k = 1)
bd_lorenz <- function(p, a, b, c, k, rel.tol = 1e-10) {
  if (b <= 1) {
    warning("The Lorenz curve needs a finite mean, so b > 1. NA returned.",
            call. = FALSE)
    return(rep(NA_real_, length(p)))
  }
  mu <- bd_moments(1, a, b, c, k, rel.tol = rel.tol)
  q  <- qbetadanish(p, a, b, c, k)
  out <- bd_incomplete_moment(1, q, a, b, c, k, rel.tol = rel.tol) / mu
  out[is.na(p)] <- NA_real_
  out
}

#' @rdname bd_lorenz
#' @export
bd_bonferroni <- function(p, a, b, c, k, rel.tol = 1e-10) {
  L <- bd_lorenz(p, a, b, c, k, rel.tol = rel.tol)
  out <- L / as.numeric(p)
  out[!is.na(p) & p <= 0] <- NA_real_
  out
}

#' Probability Weighted Moments
#'
#' \eqn{\tau_{r,s,w} = E\{Z^r F(Z)^s (1 - F(Z))^w\}}.
#'
#' @param r Power on \eqn{Z}.
#' @param s Power on \eqn{F(Z)}.
#' @param w Power on \eqn{1 - F(Z)}.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single number.
#'
#' @details
#' Existence follows the same rule as the raw moments in the worst case,
#' \eqn{b > r}, and the weights can only reduce the tail contribution.
#'
#' @export
#'
#' @examples
#' # tau_{1,0,0} is the mean
#' bd_pwm(1, 0, 0, a = 1.5, b = 3, c = 2, k = 1)
#' bd_moments(1, a = 1.5, b = 3, c = 2, k = 1)
#'
#' bd_pwm(1, 1, 0, a = 1.5, b = 3, c = 2, k = 1)
bd_pwm <- function(r, s = 0, w = 0, a, b, c, k, rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (b <= r && s == 0 && w == 0) return(Inf)
  lB <- lbeta(a, b)

  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    Fz <- stats::pbeta(u, a, b)
    rz <- if (r == 0) 0 else r * log(z)
    lg <- rz + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
    out <- exp(lg) * Fz^s * (1 - Fz)^w
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}
