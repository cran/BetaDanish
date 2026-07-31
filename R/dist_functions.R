#' The Beta-Danish Distribution
#'
#' Density, distribution function, quantile function, hazard function, survival
#' function and random generation for the four-parameter Beta-Danish
#' distribution.
#'
#' @param x,q Vector of quantiles (time points).
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param a Shape parameter (beta generator). Set `a = 1` for the
#'   three-parameter Exponentiated Danish (ED) submodel.
#' @param b Shape parameter (beta generator). Governs the upper tail: the
#'   survival function is regularly varying with index `-b`.
#' @param c Shape parameter (baseline shape).
#' @param k Scale parameter (baseline scale).
#' @param log,log.p Logical; if `TRUE`, densities/probabilities are returned on
#'   the log scale.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @details
#' With baseline \eqn{G(t) = \{kt/(1+kt)\}^{c}}, the Beta-Danish CDF is the
#' regularised incomplete beta function \eqn{F(t) = I_{G(t)}(a, b)} and the
#' density is
#' \deqn{f(t) = \frac{c\,k^{ca}\,t^{ca-1}}{B(a,b)\,(1+kt)^{ca+1}}
#'              \bigl\{1 - G(t)\bigr\}^{b-1}.}
#' The family accommodates decreasing, increasing, unimodal and bathtub-shaped
#' hazard rates.
#'
#' @section Numerical notes:
#' All three of the density, distribution and quantile functions are evaluated
#' so as to retain full relative precision in the upper tail, which is where
#' this family is distinctive and where naive implementations fail:
#'
#' * \eqn{\log\{1 - G(t)\}} is formed as
#'   `log(-expm1(c * -log1p(1/(k*t))))`. Writing \eqn{\log G} as
#'   \eqn{-c\,\log(1 + 1/kt)} avoids the cancellation in
#'   \eqn{\log(kt) - \log(1+kt)}, which loses all significant digits once
#'   \eqn{kt \gtrsim 10^{15}}.
#' * The survival function uses the beta mirror identity
#'   \eqn{1 - I_{y}(a,b) = I_{1-y}(b,a)}, so no probability near one is ever
#'   subtracted from one.
#' * `qbetadanish` obtains \eqn{1 - u} directly via
#'   \eqn{1 - q\beta(p; a, b) = q\beta(p; b, a)} with the tail flag reversed,
#'   so `1 - p` is never formed.
#'
#' Parameters recycle element-wise against `x`, `q` or `p` following the usual
#' convention for R distribution functions, so a per-observation scale (as
#' produced by an AFT link) may be supplied directly.
#'
#' As \eqn{t \to \infty}, \eqn{S(t) \propto t^{-b}}; consequently
#' \eqn{E(X^{r})} is finite if and only if \eqn{b > r}.
#'
#' @return
#' `dbetadanish` gives the density, `pbetadanish` the distribution function,
#' `qbetadanish` the quantile function, `sbetadanish` the survival function,
#' `hbetadanish` the hazard function, and `rbetadanish` generates random
#' deviates. Length is the maximum of the lengths of the numeric arguments.
#'
#' @references
#' Ahmad, B., & Danish, M. Y. (2025). Development and characterization of a
#' flexible three-parameter lifetime distribution: theoretical properties and
#' real-world applications. *Journal of Applied Mathematics, Statistics and
#' Informatics*, 21(1). \doi{10.2478/jamsi-2025-0010}
#'
#' @examples
#' dbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' pbetadanish(q = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' qbetadanish(p = 0.5, a = 1.5, b = 2, c = 3, k = 0.5)
#' hbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#' rbetadanish(n = 10, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # Per-observation scale, as used by the AFT link
#' dbetadanish(c(1, 2, 3), a = 1, b = 2, c = 1.5, k = c(0.5, 1, 2))
#'
#' # The survival tail is regularly varying with index -b
#' round(log(sbetadanish(c(1e10, 1e12), 1.5, 3, 2, 1)), 4)
#'
#' @name BetaDanish
NULL

## Recycle every numeric argument to the common length. Returns NULL if the
## common length is zero.
.bd_conform <- function(x, a, b, c, k) {
  lens <- c(length(x), length(a), length(b), length(c), length(k))
  ## R's distribution functions return a zero-length result when ANY
  ## argument is zero-length, rather than recycling up to length one.
  if (any(lens == 0L)) return(NULL)
  n <- max(lens)
  list(n = n,
       x = rep_len(as.numeric(x), n),
       a = rep_len(as.numeric(a), n),
       b = rep_len(as.numeric(b), n),
       c = rep_len(as.numeric(c), n),
       k = rep_len(as.numeric(k), n))
}

## Element-wise parameter validity.
.bd_bad <- function(p) {
  with(p, is.na(a) | is.na(b) | is.na(c) | is.na(k) |
          a <= 0 | b <= 0 | c <= 0 | k <= 0)
}

## log G(t) = -c * log1p(1/(k t)).  Full relative precision for every k*t > 0.
.bd_logG <- function(kt, c) -c * log1p(1 / kt)

## log{1 - G(t)}, stable for G arbitrarily close to either 0 or 1.
.bd_log1mG <- function(kt, c) log(-expm1(.bd_logG(kt, c)))

#' @rdname BetaDanish
#' @export
dbetadanish <- function(x, a, b, c, k, log = FALSE) {
  p <- .bd_conform(x, a, b, c, k)
  if (is.null(p)) return(numeric(0))
  out <- rep(-Inf, p$n)

  bad <- .bd_bad(p)
  if (any(bad)) {
    out[bad] <- NaN
    warning("Parameters a, b, c and k must be strictly positive; NaN returned ",
            "for ", sum(bad), " element(s).", call. = FALSE)
  }
  out[!bad & is.na(p$x)] <- NA_real_

  ok <- !bad & is.finite(p$x) & p$x > 0
  if (any(ok)) {
    xx <- p$x[ok]; aa <- p$a[ok]; bb <- p$b[ok]
    cc <- p$c[ok]; kk <- p$k[ok]
    kt <- kk * xx
    out[ok] <- log(cc) + cc * aa * log(kk) +
               (cc * aa - 1) * log(xx) -
               (cc * aa + 1) * log1p(kt) +
               (bb - 1) * .bd_log1mG(kt, cc) -
               lbeta(aa, bb)
  }
  if (log) out else exp(out)
}

#' @rdname BetaDanish
#' @export
pbetadanish <- function(q, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  p <- .bd_conform(q, a, b, c, k)
  if (is.null(p)) return(numeric(0))

  out <- rep(if (lower.tail) 0 else 1, p$n)             # value at q <= 0
  out[is.infinite(p$x) & p$x > 0] <- if (lower.tail) 1 else 0

  bad <- .bd_bad(p)
  out[bad] <- NaN
  out[!bad & is.na(p$x)] <- NA_real_

  ok <- !bad & is.finite(p$x) & p$x > 0
  if (any(ok)) {
    xx <- p$x[ok]; aa <- p$a[ok]; bb <- p$b[ok]
    cc <- p$c[ok]; kk <- p$k[ok]
    kt <- kk * xx
    out[ok] <- if (lower.tail) {
      stats::pbeta(exp(.bd_logG(kt, cc)), aa, bb,
                   lower.tail = TRUE, log.p = log.p)
    } else {
      ## Beta mirror: S(t) = 1 - I_G(a,b) = I_{1-G}(b,a). No cancellation.
      stats::pbeta(-expm1(.bd_logG(kt, cc)), bb, aa,
                   lower.tail = TRUE, log.p = log.p)
    }
  }
  if (log.p) {
    seeded <- !ok & !bad & !is.na(p$x)
    out[seeded] <- log(out[seeded])
  }
  out
}

#' @rdname BetaDanish
#' @export
qbetadanish <- function(p, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  cf <- .bd_conform(p, a, b, c, k)
  if (is.null(cf)) return(numeric(0))
  pp <- cf$x
  out <- rep(NA_real_, cf$n)

  bad <- .bd_bad(cf)
  out[bad] <- NaN

  ok <- !bad & !is.na(pp)
  if (any(ok)) {
    aa <- cf$a[ok]; bb <- cf$b[ok]; cc <- cf$c[ok]; kk <- cf$k[ok]
    ## w = 1 - u obtained directly; qbeta signals NaN for p outside its range.
    w  <- stats::qbeta(pp[ok], bb, aa, lower.tail = !lower.tail, log.p = log.p)
    lv <- log1p(-w) / cc                    # log v, v = u^(1/c)
    out[ok] <- exp(lv - log(-expm1(lv))) / kk
  }
  out
}

#' @rdname BetaDanish
#' @export
rbetadanish <- function(n, a, b, c, k) {
  if (length(n) > 1L) n <- length(n)
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0)
    stop("'n' must be a single non-negative number.", call. = FALSE)
  qbetadanish(stats::runif(n), a, b, c, k)
}

#' @rdname BetaDanish
#' @export
sbetadanish <- function(x, a, b, c, k, log = FALSE) {
  pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = log)
}

#' @rdname BetaDanish
#' @export
hbetadanish <- function(x, a, b, c, k, log = FALSE) {
  log_f <- dbetadanish(x, a, b, c, k, log = TRUE)
  log_s <- pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = TRUE)
  log_h <- log_f - log_s
  ## An exhausted survival with positive density is a divergent hazard, not a
  ## large finite one. 0/0 stays NaN rather than being silently invented.
  log_h[is.finite(log_f) & is.infinite(log_s) & log_s < 0] <- Inf
  if (log) log_h else exp(log_h)
}
