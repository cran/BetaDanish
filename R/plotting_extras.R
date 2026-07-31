#' Cox-Snell Residual Plot for AFT and Cure Fits
#'
#' Diagnostic Cox-Snell residual plot for a fitted AFT or cure model. Under a
#' correctly specified model the residuals behave like a unit exponential
#' sample, so the Kaplan-Meier estimate of their cumulative hazard should
#' follow the 45-degree line.
#'
#' @param x A fitted `"bd_aft"` or `"bd_cure"` object.
#' @param ... Further graphical parameters passed to `plot`.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 300
#' xcov  <- stats::rnorm(n)
#' k_i   <- exp(-0.5 - 0.3 * xcov)
#' t_sim <- rbetadanish(n, a = 1, b = 2, c = 1.5, k = k_i)
#' dat   <- data.frame(time = t_sim,
#'                     status = stats::rbinom(n, 1, 0.85),
#'                     x = xcov)
#' fit <- fit_bd_aft(survival::Surv(time, status) ~ x, data = dat, n_starts = 5)
#' plot(fit)
#' }
#'
#' @export
plot.bd_aft <- function(x, ...) {
  .bd_aft_or_cure_coxsnell(x, ...)
}

#' @rdname plot.bd_aft
#' @export
plot.bd_cure <- function(x, ...) {
  .bd_aft_or_cure_coxsnell(x, ...)
}

#' Natural-Scale Shape Parameters from an AFT or Cure Fit
#'
#' `fit_bd_aft()` and `fit_bd_cure()` optimise over `log_b` and `log_c`, so the
#' shape parameters must be looked up under those names and exponentiated.
#' Reading `coefficients["b"]` returns NA, which is what previously disabled
#' the Cox-Snell plots.
#'
#' @param fit A `"bd_aft"` or `"bd_cure"` object.
#' @return A list with numeric `b`, `c` and `delta`.
#' @noRd
.bd_shape_natural <- function(fit) {
  cf <- fit$coefficients
  pick <- function(nm) {
    if (paste0("log_", nm) %in% names(cf)) return(exp(as.numeric(cf[[paste0("log_", nm)]])))
    if (nm %in% names(cf))                 return(as.numeric(cf[[nm]]))
    NA_real_
  }
  delta_names <- grep("^delta_", names(cf), value = TRUE)
  list(b     = pick("b"),
       c     = pick("c"),
       delta = as.numeric(cf[delta_names]))
}

.bd_aft_or_cure_coxsnell <- function(x, ...) {
  d      <- x$data
  time   <- d$time
  status <- d$status
  X      <- d$X

  pars  <- .bd_shape_natural(x)
  b     <- pars$b
  cpar  <- pars$c
  delta <- pars$delta

  if (!is.finite(b) || !is.finite(cpar) || length(delta) == 0L ||
      any(!is.finite(delta))) {
    warning("Cox-Snell plot: the fit contains non-finite coefficients, so ",
            "residuals cannot be computed.", call. = FALSE)
    return(invisible(x))
  }
  if (length(delta) != ncol(X))
    stop("The coefficient layout does not match the design matrix.",
         call. = FALSE)

  k_i <- exp(as.numeric(X %*% delta))

  ## Cox-Snell residual is the fitted cumulative hazard, -log S(t).
  r <- -pbetadanish(time, a = 1, b = b, c = cpar, k = k_i,
                    lower.tail = FALSE, log.p = TRUE)

  keep <- is.finite(r) & r > 0
  if (!any(keep)) {
    warning("Cox-Snell plot: every residual is non-finite; nothing to plot.",
            call. = FALSE)
    return(invisible(x))
  }
  if (sum(keep) < length(r))
    warning(sprintf("Cox-Snell plot: dropped %d non-finite residual(s).",
                    sum(!keep)), call. = FALSE)

  r      <- r[keep]
  status <- status[keep]

  km_r <- survival::survfit(survival::Surv(r, status) ~ 1)
  H_r  <- -log(pmax(km_r$surv, 1e-12))

  graphics::plot(km_r$time, H_r, type = "s",
                 xlab = "Cox-Snell residual",
                 ylab = "Estimated cumulative hazard",
                 main = "Cox-Snell residuals", lwd = 2, ...)
  graphics::abline(0, 1, col = "red", lwd = 2, lty = 2)
  graphics::legend("topleft",
                   legend = c("KM cumulative hazard of residuals", "y = x"),
                   col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")
  invisible(x)
}
