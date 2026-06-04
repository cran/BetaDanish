#' Shannon Entropy of the Beta-Danish Distribution
#'
#' Computes the differential Shannon entropy
#' \eqn{H(f) = -\int_0^\infty f(t) \log f(t)\, dt}
#' for the four-parameter Beta-Danish distribution by adaptive Gauss-Kronrod
#' quadrature on the log-pdf.
#'
#' @param a,b,c,k Positive parameters of the Beta-Danish distribution.
#' @param subdivisions,rel.tol Passed to \code{stats::integrate}.
#'
#' @return Scalar Shannon entropy (in nats); \code{NA_real_} on integration
#'   failure.
#'
#' @examples
#' bd_entropy_shannon(a = 1.5, b = 2.5, c = 2, k = 1)
#'
#' @export
bd_entropy_shannon <- function(a, b, c, k,
                               subdivisions = 2000, rel.tol = 1e-8) {
  if (!check_positive_params(a, b, c, k))
    stop("All parameters must be strictly positive.")
  integrand <- function(t) {
    lp <- dbetadanish(t, a, b, c, k, log = TRUE)
    out <- rep(0, length(t))
    finite <- is.finite(lp) & t > 0
    if (any(finite)) out[finite] <- exp(lp[finite]) * lp[finite]
    out
  }
  val <- tryCatch(
    stats::integrate(integrand, lower = 0, upper = Inf,
                     subdivisions = subdivisions, rel.tol = rel.tol)$value,
    error = function(e) NA_real_)
  if (is.na(val)) return(NA_real_)
  -val
}
