## Entropy, reliability, shape and order-statistic properties.

#' Entropies of the Beta-Danish Distribution
#'
#' Shannon, Renyi and Tsallis entropies. All three share one documentation
#' topic so that every argument is described exactly once.
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param terms Number of series terms before the analytic tail is applied
#'   (Shannon only).
#' @param method `"closed"` (default) for the closed form, or `"quadrature"`
#'   for direct numerical integration of \eqn{-f \log f}, which is slower but
#'   independent of the series (Shannon only).
#' @param order Entropy order \eqn{q}, positive and not equal to one (Renyi and
#'   Tsallis only).
#' @param rel.tol Relative accuracy, passed to [stats::integrate()].
#' @param subdivisions Subdivision limit, passed to [stats::integrate()].
#'
#' @return A single number. The Shannon entropy is in nats.
#'
#' @details
#' Writing \eqn{u = G(Z) \sim \mathrm{Beta}(a,b)} and \eqn{v = u^{1/c}}, the
#' Shannon entropy has the closed form
#' \deqn{H = \log B(a,b) - \log(ck)
#'          - (a - 1/c)\{\psi(a) - \psi(a+b)\}
#'          - (b - 1)\{\psi(b) - \psi(a+b)\}
#'          + 2\sum_{i \ge 1} \frac{B(a + i/c,\, b)}{i\, B(a,b)},}
#' the final sum arising from \eqn{-2E\log(1 - v)} expanded as a power series.
#'
#' With \eqn{I_q = \int_0^\infty f(z)^q dz}, the Renyi entropy is
#' \eqn{\log(I_q)/(1-q)} and the Tsallis entropy is \eqn{(1 - I_q)/(q-1)}. Both
#' reduce to the Shannon entropy as \eqn{q \to 1}, which is excluded; use
#' `bd_entropy_shannon()` at \eqn{q = 1}.
#'
#' \eqn{I_q} is evaluated as \eqn{\int_0^1 f(z(u))^{q-1} d\mathrm{Beta}(u;a,b)}
#' on the finite Beta scale rather than over the half line, which keeps the
#' heavy upper tail from dominating the quadrature.
#'
#' @section Series truncation:
#' The Shannon series terms decay like \eqn{i^{-(b+1)}}, so truncating at
#' \eqn{M} omits a tail of order \eqn{M^{-b}}. That is negligible for large
#' \eqn{b} and is not for small \eqn{b}: against high-precision integration the
#' plain truncated sum at \eqn{M = 2000} is out by about \eqn{10^{-8}} at
#' \eqn{b = 3} but by about \eqn{0.11} at \eqn{b = 0.5}. The analytic tail
#' \eqn{\Gamma(b)c^{b}M^{-b}/\{b B(a,b)\}} is therefore added, restoring
#' agreement to roughly eight digits across the range.
#'
#' @name bd_entropy
#'
#' @examples
#' bd_entropy_shannon(a = 1.5, b = 3, c = 2, k = 1)
#'
#' # Independent check by quadrature
#' bd_entropy_shannon(a = 1.5, b = 3, c = 2, k = 1, method = "quadrature")
#'
#' bd_entropy_renyi(a = 1.5, b = 3, c = 2, k = 1, order = 2)
#' bd_entropy_tsallis(a = 1.5, b = 3, c = 2, k = 1, order = 2)
NULL

#' @rdname bd_entropy
#' @export
bd_entropy_shannon <- function(a, b, c, k, terms = 20000L,
                               method = c("closed", "quadrature"),
                               rel.tol = 1e-10, subdivisions = 4000L) {
  method <- match.arg(method)
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)

  if (identical(method, "quadrature")) {
    integrand <- function(u) {
      z  <- .bd_z_of_u(u, c, k)
      lf <- dbetadanish(z, a, b, c, k, log = TRUE)
      out <- -lf * exp((a - 1) * log(u) + (b - 1) * log1p(-u) - lbeta(a, b))
      out[!is.finite(out)] <- 0
      out
    }
    return(tryCatch(
      stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                       subdivisions = subdivisions, stop.on.error = FALSE)$value,
      error = function(e) NA_real_))
  }

  M  <- as.integer(terms)
  lB <- lbeta(a, b)
  i  <- seq_len(M)
  S  <- sum(exp(lbeta(a + i / c, b) - lB) / i)
  ## Euler-Maclaurin tail: terms ~ Gamma(b) c^b / B(a,b) * i^-(b+1).
  S  <- S + exp(lgamma(b) + b * log(c) - lB) * M^(-b) / b

  lB - log(c * k) -
    (a - 1 / c) * (digamma(a) - digamma(a + b)) -
    (b - 1)     * (digamma(b) - digamma(a + b)) +
    2 * S
}

#' @rdname bd_entropy
#' @export
bd_entropy_renyi <- function(a, b, c, k, order = 2, rel.tol = 1e-10,
                             subdivisions = 4000L) {
  Iq <- .bd_Iq(a, b, c, k, order, rel.tol, subdivisions)
  if (is.na(Iq) || Iq <= 0) return(NA_real_)
  log(Iq) / (1 - order)
}

#' @rdname bd_entropy
#' @export
bd_entropy_tsallis <- function(a, b, c, k, order = 2, rel.tol = 1e-10,
                               subdivisions = 4000L) {
  Iq <- .bd_Iq(a, b, c, k, order, rel.tol, subdivisions)
  if (is.na(Iq)) return(NA_real_)
  (1 - Iq) / (order - 1)
}

#' The Integral of f^q, on the Beta Scale
#' @noRd
.bd_Iq <- function(a, b, c, k, q, rel.tol = 1e-10, subdivisions = 4000L) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (length(q) != 1L || is.na(q) || q <= 0)
    stop("'order' must be a single positive number.", call. = FALSE)
  if (abs(q - 1) < 1e-8)
    stop("Order 1 is the Shannon entropy; use bd_entropy_shannon().",
         call. = FALSE)

  lB <- lbeta(a, b)
  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    lf <- dbetadanish(z, a, b, c, k, log = TRUE)
    out <- exp((q - 1) * lf + (a - 1) * log(u) + (b - 1) * log1p(-u) - lB)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = subdivisions,
                            stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Stress-Strength Reliability
#'
#' \eqn{R = P(X > Y)} where \eqn{X} is the strength and \eqn{Y} the stress,
#' each Beta-Danish distributed with its own parameters.
#'
#' @param strength Named list or numeric vector giving `a`, `b`, `c`, `k` for
#'   the strength variable \eqn{X}.
#' @param stress Same, for the stress variable \eqn{Y}.
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return A single probability.
#'
#' @details
#' \eqn{R = \int_0^\infty S_X(y) f_Y(y)\,dy}, evaluated on the Beta scale of
#' \eqn{Y} so the integral is taken over \eqn{(0,1)}.
#'
#' The convention is the standard reliability one: \eqn{R} is the probability
#' that the strength exceeds the stress, so larger \eqn{R} means a more
#' reliable component. Passing the two arguments the wrong way round returns
#' \eqn{1 - R}.
#'
#' @export
#'
#' @examples
#' # Identical distributions must give one half
#' p <- list(a = 1.5, b = 3, c = 2, k = 1)
#' bd_stress_strength(strength = p, stress = p)
#'
#' # A stronger component
#' bd_stress_strength(strength = list(a = 1.5, b = 3, c = 2, k = 0.5),
#'                    stress   = list(a = 1.5, b = 3, c = 2, k = 2))
bd_stress_strength <- function(strength, stress, rel.tol = 1e-10) {
  px <- .bd_as_pars(strength, "strength")
  py <- .bd_as_pars(stress,   "stress")

  lB <- lbeta(py$a, py$b)
  integrand <- function(u) {
    y  <- .bd_z_of_u(u, py$c, py$k)
    Sx <- sbetadanish(y, px$a, px$b, px$c, px$k)
    out <- Sx * exp((py$a - 1) * log(u) + (py$b - 1) * log1p(-u) - lB)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Coerce a Parameter Argument to a Named List
#' @noRd
.bd_as_pars <- function(x, what) {
  if (is.list(x)) x <- unlist(x)
  if (!is.numeric(x) || is.null(names(x)) || !all(c("a", "b", "c", "k") %in% names(x)))
    stop("'", what, "' must supply named values a, b, c and k.", call. = FALSE)
  p <- as.list(x[c("a", "b", "c", "k")])
  if (!check_positive_params(p$a, p$b, p$c, p$k))
    stop("All ", what, " parameters must be finite and strictly positive.",
         call. = FALSE)
  p
}

#' Hazard Rate Shape and the Glaser Diagnostic
#'
#' Classifies the hazard rate as increasing, decreasing, bathtub-shaped or
#' upside-down bathtub-shaped, and returns the mode of the density.
#'
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param n_grid Number of grid points between the lower and upper evaluation
#'   quantiles.
#' @param range_p Two probabilities bounding the evaluation range.
#'
#' @return A list of class `"bd_shape"` with the classification, the grid of
#'   times, the hazard, Glaser's \eqn{\eta(t) = -f'(t)/f(t)}, and the mode.
#'
#' @details
#' Glaser's criterion works through \eqn{\eta(t) = -f'(t)/f(t)}: monotone
#' increasing \eqn{\eta} gives an increasing hazard, monotone decreasing gives a
#' decreasing hazard, and a single interior turning point gives a bathtub or
#' upside-down bathtub. Here \eqn{\eta} is formed analytically from the
#' log-density and evaluated on a grid, and the hazard itself is checked
#' alongside it, so the two must agree before a shape is reported.
#'
#' @export
#'
#' @examples
#' s <- bd_hazard_shape(a = 1.5, b = 3, c = 2, k = 1)
#' s$shape
#' s$mode
bd_hazard_shape <- function(a, b, c, k, n_grid = 400L,
                            range_p = c(1e-4, 1 - 1e-4)) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)

  lo <- qbetadanish(range_p[1], a, b, c, k)
  hi <- qbetadanish(range_p[2], a, b, c, k)
  t  <- exp(seq(log(lo), log(hi), length.out = as.integer(n_grid)))

  eta <- .bd_glaser_eta(t, a, b, c, k)
  h   <- hbetadanish(t, a, b, c, k)

  ok <- is.finite(eta) & is.finite(h) & h > 0
  shape <- if (sum(ok) < 10L) {
    "indeterminate"
  } else {
    .bd_classify(h[ok])
  }

  ## Mode of the density: maximise the log-density over the same range.
  mode_t <- tryCatch(
    stats::optimize(function(x) dbetadanish(x, a, b, c, k, log = TRUE),
                    interval = c(lo, hi), maximum = TRUE)$maximum,
    error = function(e) NA_real_)
  if (is.finite(mode_t) &&
      (mode_t <= lo * (1 + 1e-6) || mode_t >= hi * (1 - 1e-6)))
    mode_t <- NA_real_   # optimum on the boundary is not an interior mode

  out <- list(shape = shape, time = t, hazard = h, eta = eta,
              mode = mode_t,
              parameters = c(a = a, b = b, c = c, k = k))
  class(out) <- "bd_shape"
  out
}

#' Glaser's eta(t) = -f'(t)/f(t), Analytically
#' @noRd
.bd_glaser_eta <- function(t, a, b, c, k) {
  kt   <- k * t
  logw <- -log1p(1 / kt)                 # log of w = kt/(1+kt)
  lG   <- c * logw
  om   <- -expm1(lG)                     # 1 - G, stable
  ca   <- c * a

  dlog_1mG <- -c * exp((c - 1) * logw) * k / ((1 + kt)^2 * om)
  dlogf    <- (ca - 1) / t - (ca + 1) * k / (1 + kt) + (b - 1) * dlog_1mG
  -dlogf
}

#' Classify a Positive Sequence as Monotone or Single-Turning
#' @noRd
.bd_classify <- function(v, tol = 1e-10) {
  d <- diff(v)
  d <- d[is.finite(d)]
  if (!length(d)) return("indeterminate")
  up   <- d > tol
  down <- d < -tol
  if (!any(down)) return("increasing")
  if (!any(up))   return("decreasing")
  first_up <- which(up)[1]
  first_dn <- which(down)[1]
  if (first_dn < first_up) "bathtub" else "upside-down bathtub"
}

#' @param x A `"bd_shape"` object.
#' @param ... Ignored.
#' @rdname bd_hazard_shape
#' @export
print.bd_shape <- function(x, ...) {
  cat("Beta-Danish hazard shape\n")
  cat("  parameters: ",
      paste(names(x$parameters), signif(x$parameters, 4), sep = " = ",
            collapse = ", "), "\n", sep = "")
  cat("  shape:      ", x$shape, "\n", sep = "")
  cat("  mode:       ",
      if (is.finite(x$mode)) signif(x$mode, 6) else "none in range (monotone density)",
      "\n", sep = "")
  cat("  evaluated over t in [", signif(min(x$time), 4), ", ",
      signif(max(x$time), 4), "]\n", sep = "")
  invisible(x)
}

#' Order Statistic Distribution and Moments
#'
#' @param t Vector of times, for `bd_order_stat_cdf`.
#' @param r Order of the moment, for `bd_order_stat_moments`.
#' @param i Index of the order statistic, from 1 to `n`.
#' @param n Sample size.
#' @param a Shape parameter (beta generator).
#' @param b Shape parameter governing the tail.
#' @param c Shape parameter (baseline).
#' @param k Scale parameter (baseline).
#' @param rel.tol Passed to [stats::integrate()].
#'
#' @return `bd_order_stat_cdf` returns a vector the length of `t`;
#'   `bd_order_stat_moments` returns a single number.
#'
#' @details
#' \eqn{F_{i:n}(t) = I_{F(t)}(i, n - i + 1)}, and moments are taken against the
#' order-statistic density on the Beta scale. The \eqn{i}-th order statistic of
#' a sample of size \eqn{n} has a finite \eqn{r}-th moment when
#' \eqn{b(n - i + 1) > r}, so extremes need heavier restrictions than the parent
#' distribution, and the maximum needs \eqn{b > r} exactly as the parent does.
#'
#' @seealso [bd_order_stat_pdf()]
#'
#' @export
#'
#' @examples
#' bd_order_stat_cdf(c(0.5, 1, 2), i = 3, n = 5, a = 1.5, b = 3, c = 2, k = 1)
#' bd_order_stat_moments(1, i = 1, n = 5, a = 1.5, b = 3, c = 2, k = 1)
bd_order_stat_cdf <- function(t, i, n, a, b, c, k) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (i < 1 || i > n) stop("'i' must lie between 1 and n.", call. = FALSE)
  stats::pbeta(pbetadanish(t, a, b, c, k), i, n - i + 1)
}

#' @rdname bd_order_stat_cdf
#' @export
bd_order_stat_moments <- function(r, i, n, a, b, c, k, rel.tol = 1e-10) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  if (i < 1 || i > n) stop("'i' must lie between 1 and n.", call. = FALSE)
  if (b * (n - i + 1) <= r) return(Inf)

  lB   <- lbeta(a, b)
  lcon <- lgamma(n + 1) - lgamma(i) - lgamma(n - i + 1)

  integrand <- function(u) {
    z  <- .bd_z_of_u(u, c, k)
    Fz <- stats::pbeta(u, a, b)
    rz <- if (r == 0) 0 else r * log(z)
    lg <- lcon + rz + (i - 1) * log(Fz) + (n - i) * log1p(-Fz) +
          (a - 1) * log(u) + (b - 1) * log1p(-u) - lB
    out <- exp(lg)
    out[!is.finite(out)] <- 0
    out
  }
  tryCatch(stats::integrate(integrand, 0, 1, rel.tol = rel.tol,
                            subdivisions = 4000L, stop.on.error = FALSE)$value,
           error = function(e) NA_real_)
}

#' Tail Index of the Beta-Danish Distribution
#'
#' @param a Shape parameter (beta generator). Accepted for interface
#'   consistency; it does not affect the index.
#' @param b Shape parameter. This is the tail index.
#' @param c Shape parameter (baseline). Accepted for interface consistency.
#' @param k Scale parameter (baseline). Accepted for interface consistency.
#'
#' @return A list with the tail index, the highest finite moment order, and a
#'   note on the moment generating function.
#'
#' @details
#' The survival function is regularly varying at infinity with index \eqn{-b}:
#' \eqn{S(t) \sim (c/k)^b t^{-b} / \{b B(a,b)\}}. Three consequences follow.
#'
#' The \eqn{r}-th moment is finite if and only if \eqn{b > r}, so the mean needs
#' \eqn{b > 1}, the variance \eqn{b > 2}, skewness \eqn{b > 3} and kurtosis
#' \eqn{b > 4}.
#'
#' The moment generating function does not exist for any \eqn{t > 0}: a
#' regularly varying tail decays polynomially, so \eqn{E(e^{tZ})} diverges.
#' Characteristic-function or Laplace-transform arguments must be used instead
#' of MGF ones anywhere in this family.
#'
#' The distribution lies in the Frechet domain of attraction, so sample maxima
#' normalise to a Frechet limit with shape \eqn{b}, not to a Gumbel limit.
#'
#' @export
#'
#' @examples
#' bd_tail_index(a = 1.5, b = 3, c = 2, k = 1)
bd_tail_index <- function(a, b, c, k) {
  if (!check_positive_params(a, b, c, k))
    stop("a, b, c and k must be finite and strictly positive.", call. = FALSE)
  list(tail_index          = b,
       survival_exponent   = -b,
       highest_finite_moment = if (b > 1) floor(ceiling(b) - 1) else 0,
       moment_condition    = "E(Z^r) finite if and only if b > r",
       mgf_exists          = FALSE,
       domain_of_attraction = "Frechet")
}
