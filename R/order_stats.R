#' Density of the r-th Order Statistic
#'
#' Evaluates the probability density function of the r-th order statistic
#' from a sample of size n drawn from the Beta-Danish distribution.
#'
#' @param x Numeric vector of time points.
#' @param r Integer order (1 = minimum, n = maximum).
#' @param n Integer sample size.
#' @param a,b,c,k Positive parameters of the Beta-Danish distribution.
#' @param log Logical; if \code{TRUE} return the log-density.
#'
#' @return Numeric vector (or its log).
#'
#' @examples
#' tgrid <- seq(0.01, 5, length.out = 50)
#' bd_order_stat_pdf(tgrid, r = 5, n = 20,
#'                   a = 1.5, b = 2.5, c = 2, k = 1)
#'
#' @export
bd_order_stat_pdf <- function(x, r, n, a, b, c, k, log = FALSE) {
  if (!check_positive_params(a, b, c, k))
    stop("All parameters must be strictly positive.")
  if (!is.numeric(r) || !is.numeric(n) || r < 1 || n < r || r != round(r) || n != round(n))
    stop("r and n must be integers with 1 <= r <= n.")
  log_norm <- -lbeta(r, n - r + 1)
  log_f <- dbetadanish(x, a, b, c, k, log = TRUE)
  log_F <- pbetadanish(x, a, b, c, k, lower.tail = TRUE, log.p = TRUE)
  log_S <- pbetadanish(x, a, b, c, k, lower.tail = FALSE, log.p = TRUE)
  out <- log_norm + (r - 1) * log_F + (n - r) * log_S + log_f
  out[!is.finite(log_f) | x <= 0] <- -Inf
  if (log) out else exp(out)
}
