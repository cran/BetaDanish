#' The Beta-Danish Distribution
#'
#' Density, distribution function, quantile function, hazard function, and
#' random generation for the four-parameter Beta-Danish distribution.
#'
#' @param x,q Vector of quantiles (time points).
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param a Shape parameter (beta generator). Set `a = 1` for the 3-parameter submodel.
#' @param b Shape parameter (beta generator / tail weight).
#' @param c Shape parameter (baseline shape).
#' @param k Scale parameter (baseline scale).
#' @param log,log.p Logical; if TRUE, probabilities/densities are given as log.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise P[X > x].
#'
#' @details
#' The Beta-Danish distribution is a highly flexible lifetime distribution capable
#' of modeling decreasing, increasing, unimodal, and bathtub-shaped hazard rates.
#'
#' @return
#' `dbetadanish` gives the density, `pbetadanish` gives the distribution function,
#' `qbetadanish` gives the quantile function, `hbetadanish` gives the hazard function,
#' and `rbetadanish` generates random deviates.
#'
#' @references
#' Ahmad, B., & Danish, M. Y. (2026). Development and Characterization of a
#' Flexible Three-Parameter Lifetime Distribution.
#'
#' @examples
#' # Density
#' dbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # CDF
#' pbetadanish(q = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # Hazard
#' hbetadanish(x = 2, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # Random generation
#' rbetadanish(n = 10, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' @name BetaDanish
NULL

#' @rdname BetaDanish
#' @export
dbetadanish <- function(x, a, b, c, k, log = FALSE) {
  if (!check_positive_params(a, b, c, k)) {
    warning("Parameters a, b, c, and k must be strictly positive.")
    return(rep(NaN, length(x)))
  }

  # Log-space calculations for numerical stability
  log_k <- log(k)
  log_x <- suppressWarnings(log(x))
  log_1_plus_kx <- log1p(k * x)

  # Baseline G(x) = (kx / (1 + kx))^c
  log_u <- log_k + log_x - log_1_plus_kx
  log_G <- c * log_u

  # Baseline density g(x) = c * k * (kx)^(c-1) / (1 + kx)^(c+1)
  log_g <- log(c) + log_k + (c - 1) * log_u - 2 * log_1_plus_kx

  # Beta-Danish PDF: f(x) = g(x) * G(x)^(a-1) * (1 - G(x))^(b-1) / B(a,b)
  G <- exp(log_G)
  G[G >= 1] <- 1 - 1e-16
  log_1_minus_G <- log1p(-G) # Stably computes log(1 - G)

  log_pdf <- log_g + (a - 1) * log_G + (b - 1) * log_1_minus_G - lbeta(a, b)

  # Boundary handling
  log_pdf[x <= 0] <- -Inf

  if (log) return(log_pdf)
  return(exp(log_pdf))
}

#' @rdname BetaDanish
#' @export
pbetadanish <- function(q, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  if (!check_positive_params(a, b, c, k)) return(rep(NaN, length(q)))

  q[q < 0] <- 0

  # Baseline G(q)
  u <- (k * q) / (1 + k * q)
  G_q <- u^c

  # F(x) = I_{G(x)}(a, b) = pbeta(G(x), a, b)
  p <- stats::pbeta(G_q, shape1 = a, shape2 = b, lower.tail = lower.tail, log.p = log.p)
  return(p)
}

#' @rdname BetaDanish
#' @export
qbetadanish <- function(p, a, b, c, k, lower.tail = TRUE, log.p = FALSE) {
  if (!check_positive_params(a, b, c, k)) return(rep(NaN, length(p)))

  # Invert the Beta generator
  y <- stats::qbeta(p, shape1 = a, shape2 = b, lower.tail = lower.tail, log.p = log.p)

  # Invert the baseline: G(x) = y  =>  x = y^(1/c) / (k * (1 - y^(1/c)))
  y_c <- y^(1/c)
  x <- y_c / (k * (1 - y_c))

  return(x)
}

#' @rdname BetaDanish
#' @export
rbetadanish <- function(n, a, b, c, k) {
  if (!check_positive_params(a, b, c, k)) stop("Parameters must be positive.")
  u <- stats::runif(n)
  qbetadanish(u, a, b, c, k)
}

#' @rdname BetaDanish
#' @export
hbetadanish <- function(x, a, b, c, k, log = FALSE) {
  log_pdf <- dbetadanish(x, a, b, c, k, log = TRUE)
  log_surv <- pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = TRUE)

  log_surv[log_surv == -Inf] <- -700
  log_haz <- log_pdf - log_surv

  if (log) return(log_haz)
  return(exp(log_haz))
}
