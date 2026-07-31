#' The Exponentiated Danish Distribution
#'
#' Density, distribution function, quantile function, survival function, hazard
#' function and random generation for the three-parameter Exponentiated Danish
#' distribution, the \eqn{a = 1} submodel of the Beta-Danish family.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param b,c,k Shape, shape and scale parameters.
#' @param log,log.p Logical; return values on the log scale.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P[X \le x]}.
#'
#' @return Vectors of the same form as the Beta-Danish equivalents.
#'
#' @details
#' These are thin wrappers that fix \eqn{a = 1}, provided because the ED
#' submodel is a named distribution in its own right in the underlying work and
#' writing `dbetadanish(x, 1, b, c, k)` obscures that.
#'
#' At \eqn{a = 1} the beta generator collapses and the distribution function
#' simplifies to \eqn{F(t) = 1 - \{1 - G(t)\}^{b}} with
#' \eqn{G(t) = \{kt/(1+kt)\}^{c}}. The upper tail still has index \eqn{-b}, so
#' the moment condition \eqn{b > r} is unchanged.
#'
#' @seealso [BetaDanish] for the four-parameter parent,
#'   [fit_betadanish()] with `submodel = TRUE` for estimation.
#'
#' @name ExponentiatedDanish
#'
#' @examples
#' ded(2, b = 3, c = 2, k = 1)
#' ped(2, b = 3, c = 2, k = 1)
#' qed(0.5, b = 3, c = 2, k = 1)
#' red(5, b = 3, c = 2, k = 1)
#'
#' # Identical to the a = 1 parent
#' all.equal(ded(2, 3, 2, 1), dbetadanish(2, 1, 3, 2, 1))
NULL

#' @rdname ExponentiatedDanish
#' @export
ded <- function(x, b, c, k, log = FALSE) {
  dbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}

#' @rdname ExponentiatedDanish
#' @export
ped <- function(q, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  pbetadanish(q, a = 1, b = b, c = c, k = k,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ExponentiatedDanish
#' @export
qed <- function(p, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  qbetadanish(p, a = 1, b = b, c = c, k = k,
              lower.tail = lower.tail, log.p = log.p)
}

#' @rdname ExponentiatedDanish
#' @export
red <- function(n, b, c, k) {
  rbetadanish(n, a = 1, b = b, c = c, k = k)
}

#' @rdname ExponentiatedDanish
#' @export
sed <- function(x, b, c, k, log = FALSE) {
  sbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}

#' @rdname ExponentiatedDanish
#' @export
hed <- function(x, b, c, k, log = FALSE) {
  hbetadanish(x, a = 1, b = b, c = c, k = k, log = log)
}
