#' Kullback-Leibler Separation from the Exponentiated Danish Submodel
#'
#' Computes the Kullback-Leibler divergence from a Beta-Danish distribution to
#' its closest Exponentiated Danish (ED) submodel, and the sample size at which
#' the four-parameter model would be expected to win by a given AIC margin.
#'
#' @param a,b,c,k Parameters of the Beta-Danish distribution, each a single
#'   positive number.
#' @param tgrid Optional numeric vector of equally spaced points on the scale
#'   \eqn{t = \log z} for the trapezoid rule. `NULL` (default) builds a grid
#'   with spacing 0.01 that covers the distribution from its
#'   \eqn{10^{-15}} to its \eqn{1 - 10^{-15}} quantile. The thesis table used
#'   `seq(-80, 80, by = 0.01)`.
#' @param box Numeric vector `c(lower, upper)` bounding the ED parameters
#'   during the minimisation. Default `c(1e-4, 1e4)`, the thesis box.
#' @param aic_margin AIC margin for the sample-size statement. Default 10.
#'
#' @return A list with components
#'   \describe{
#'     \item{`kl`}{The minimised divergence \eqn{\mathrm{KL}(f_{BD}, f_{ED^*})}.}
#'     \item{`ed`}{Named vector `b`, `c`, `k` of the closest ED distribution.}
#'     \item{`n_aic`}{The sample size
#'       \eqn{n = (m + 2)/(2\,\mathrm{KL})}, with \eqn{m} the AIC margin,
#'       at which the expected log-likelihood advantage of the
#'       Beta-Danish over its closest submodel exceeds the AIC margin
#'       (\eqn{6/\mathrm{KL}} for a margin of 10).}
#'     \item{`mass`}{The integrated Beta-Danish density on the grid; it should
#'       be 1, and a warning is given if it is not.}
#'     \item{`at_box_edge`}{Logical; whether the closest ED distribution sits on
#'       the edge of `box`, in which case `kl` is an upper bound.}
#'     \item{`tgrid_range`}{The range of the integration grid on the log
#'       scale.}
#'   }
#'
#' @details
#' The divergence is
#' \deqn{\mathrm{KL}(f_{BD}, f_{ED^*}) = \min_{b', c', k'} \int_0^\infty
#'   f_{BD}(z) \log\frac{f_{BD}(z)}{f_{ED}(z; b', c', k')}\,dz,}
#' the divergence from the parent to its closest submodel. Since
#' \eqn{n\,\mathrm{KL}} is the expected log-likelihood advantage of the parent
#' in a sample of size \eqn{n}, and the parent pays 2 AIC units for its extra
#' parameter, a margin of \eqn{m} AIC units needs
#' \eqn{n = (m + 2)/(2\,\mathrm{KL})} observations.
#'
#' The integral is taken on the scale \eqn{t = \log z}, where both tails decay
#' exponentially, by the trapezoid rule. With \eqn{s = \log k + t} and
#' \eqn{\log w = -\log(1 + e^{-s})}, the log-density of \eqn{T = \log Z} is
#' written in the form, stable over the whole real line,
#' \deqn{\log f_T(t) = \log c + (ca + 1)\log w - s + (b - 1)\log(1 - w^{c})
#'   - \log B(a, b).}
#' The minimisation over the ED parameters uses 27 Nelder-Mead starts around
#' the tail-matched submodel (\eqn{b' = b}, \eqn{c' = ac}, \eqn{k' = k}) and a
#' final polishing run. This is the computation behind the Kullback-Leibler
#' table of the thesis, which replaced an estimate from a single simulated
#' sample.
#'
#' @references
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan.
#'
#' @seealso [compare_models()] for the likelihood ratio test of \eqn{a = 1}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' kl <- bd_kl_submodel(a = 3, b = 2, c = 2, k = 1)
#' kl$kl      # divergence to the closest ED distribution
#' kl$n_aic   # sample size for a 10-unit AIC margin
#' }
bd_kl_submodel <- function(a, b, c, k, tgrid = NULL, box = c(1e-4, 1e4),
                           aic_margin = 10) {
  vals <- c(a = a, b = b, c = c, k = k)
  if (length(a) != 1L || length(b) != 1L || length(c) != 1L ||
      length(k) != 1L || any(!is.finite(vals)) || any(vals <= 0))
    stop("'a', 'b', 'c' and 'k' must each be a single positive number.",
         call. = FALSE)
  if (!is.numeric(box) || length(box) != 2L || any(!is.finite(box)) ||
      box[1] <= 0 || box[2] <= box[1])
    stop("'box' must be two positive numbers, c(lower, upper), with ",
         "lower < upper.", call. = FALSE)
  if (!is.numeric(aic_margin) || length(aic_margin) != 1L ||
      !is.finite(aic_margin) || aic_margin < 0)
    stop("'aic_margin' must be a single non-negative number.", call. = FALSE)

  if (is.null(tgrid)) {
    tgrid <- .bd_kl_grid(a, b, c, k)
  } else {
    tgrid <- sort(as.numeric(tgrid[is.finite(tgrid)]))
    if (length(tgrid) < 3L)
      stop("'tgrid' needs at least three finite values.", call. = FALSE)
  }

  h  <- tgrid[2L] - tgrid[1L]
  wq <- rep(h, length(tgrid)); wq[c(1L, length(tgrid))] <- h / 2
  lp <- .bd_logfT(tgrid, a, b, c, k)
  pw <- wq * exp(lp)
  use <- is.finite(lp) & pw > 0
  mass <- sum(pw[use])
  if (abs(mass - 1) > 1e-6)
    warning(sprintf(paste0("The Beta-Danish density integrates to %.8f on ",
                           "the grid, not 1; widen 'tgrid' for an accurate ",
                           "divergence."), mass), call. = FALSE)

  obj <- function(eta) {
    th <- exp(eta)
    if (any(!is.finite(th)) || any(th < box[1]) || any(th > box[2])) return(1e6)
    lq <- .bd_logfT(tgrid[use], 1, th[1], th[2], th[3])
    v  <- sum(pw[use] * (lp[use] - lq))
    if (is.finite(v)) v else 1e6
  }

  ## 27 starts around the tail-matched submodel (b' = b, c' = a c, k' = k).
  zmed <- qbetadanish(0.5, a, b, c, k)
  best <- NULL
  for (fb in c(0.6, 1, 1.6)) for (fc in c(0.3, 1, 2.5)) for (fk in c(0.4, 1, 2.5)) {
    s <- c(b * fb, a * c * fc, if (fc != 1) fk / zmed else k * fk)
    s <- pmin(pmax(s, box[1] * 2), box[2] / 2)
    o <- tryCatch(stats::optim(log(s), obj, method = "Nelder-Mead",
                               control = list(maxit = 3000, reltol = 1e-14)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && o$value < 1e5 &&
        (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best))
    stop("The minimisation over the ED parameters failed from every start.",
         call. = FALSE)
  o <- tryCatch(stats::optim(best$par, obj, method = "Nelder-Mead",
                             control = list(maxit = 5000, reltol = 1e-15)),
                error = function(e) NULL)
  if (!is.null(o) && is.finite(o$value) && o$value < best$value) best <- o

  ed <- stats::setNames(exp(best$par), c("b", "c", "k"))
  kl <- max(best$value, 0)
  list(kl = kl,
       ed = ed,
       n_aic = if (kl > 1e-12) (aic_margin + 2) / (2 * kl) else Inf,
       mass = mass,
       at_box_edge = any(ed <= box[1] * 2.1 | ed >= box[2] / 2.1),
       tgrid_range = range(tgrid))
}

#' Log-Density of log Z for the Beta-Danish Distribution
#'
#' Stable over the whole real line: \eqn{\log w} is formed as
#' \eqn{-\log(1 + e^{-s})}, and where \eqn{w^c} rounds to one the term
#' \eqn{\log(1 - w^c)} is replaced by its limit
#' \eqn{\log c + \log(1 - w)}, with \eqn{\log(1 - w) = -\log(1 + e^{s})}.
#'
#' @noRd
.bd_logfT <- function(tg, a, b, cc, k) {
  s   <- log(k) + tg
  lw  <- -log1p(exp(-s))
  l1m <- log(-expm1(cc * lw))
  bad <- !is.finite(l1m)
  if (any(bad)) {
    sb <- s[bad]
    l1m[bad] <- log(cc) - ifelse(sb > 0, sb + log1p(exp(-sb)), log1p(exp(sb)))
  }
  log(cc) + (cc * a + 1) * lw - s + (b - 1) * l1m - lbeta(a, b)
}

#' Integration Grid on the log Scale
#'
#' Covers the \eqn{10^{-15}} and \eqn{1 - 10^{-15}} quantiles, padded by one
#' unit, with spacing 0.01 (coarser only if that would exceed 200000 points).
#' The quantiles are found on the log scale so that neither tail underflows.
#'
#' @noRd
.bd_kl_grid <- function(a, b, cc, k, p = 1e-15) {
  ## lower quantile: u = qbeta(p, a, b), leading-term inversion on underflow
  u <- stats::qbeta(p, a, b)
  lu <- if (is.finite(u) && u > 0) log(u) else (log(p) + log(a) + lbeta(a, b)) / a
  lv  <- lu / cc
  lo  <- lv - log(-expm1(lv)) - log(k)
  ## upper quantile: w = 1 - u = qbeta(p, b, a) by the beta mirror
  w  <- stats::qbeta(p, b, a)
  lwq <- if (is.finite(w) && w > 0) log(w) else (log(p) + log(b) + lbeta(a, b)) / b
  lv2 <- if (is.finite(w) && w > 0) log1p(-w) / cc else -exp(lwq) / cc
  l1v <- if (is.finite(w) && w > 0) log(-expm1(lv2)) else lwq - log(cc)
  hi  <- lv2 - l1v - log(k)
  if (!is.finite(lo)) lo <- -80
  if (!is.finite(hi)) hi <- 80
  lo <- max(lo - 1, -700 - log(k))
  hi <- min(hi + 1, 700 - log(k))
  if (hi <= lo) { lo <- -80; hi <- 80 }
  n <- min(floor((hi - lo) / 0.01) + 1, 2e5)
  seq(lo, hi, length.out = max(n, 3))
}
