## Diagnostic and exploratory plots.

#' Scaled Total Time on Test Plot
#'
#' Draws the scaled total time on test transform, a distribution-free way of
#' judging hazard shape before any model is fitted.
#'
#' @param time Numeric vector of observed times, or a fitted `"betadanish"`
#'   object, from which the times are taken.
#' @param status Optional event indicator. Censored observations are dropped,
#'   since the transform is defined for complete samples.
#' @param add Logical; add to an existing plot rather than starting a new one.
#' @param col,lwd Colour and line width for the curve.
#' @param main,xlab,ylab Labels.
#' @param ... Further graphical parameters.
#'
#' @return Invisibly, a data frame with the plotting coordinates `i_n` and
#'   `phi`, and an attribute `"shape"` giving the suggested hazard shape.
#'
#' @details
#' For an ordered sample \eqn{x_{(1)} \le \cdots \le x_{(n)}}, the scaled
#' transform at \eqn{i/n} is
#' \deqn{\phi(i/n) = \frac{\sum_{j=1}^{i} x_{(j)} + (n-i)x_{(i)}}
#'                        {\sum_{j=1}^{n} x_{(j)}}.}
#'
#' Read it against the diagonal. A curve entirely above the diagonal indicates
#' an increasing hazard, entirely below a decreasing one; a curve that starts
#' below and crosses above suggests a bathtub shape, and the reverse suggests a
#' unimodal one. A curve close to the diagonal indicates a constant hazard,
#' that is an exponential sample.
#'
#' This is a shape diagnostic, not a test. It is worth drawing before choosing
#' between the four-parameter model and its submodel, because it says which
#' hazard shapes the data can support without assuming any of them.
#'
#' @seealso [bd_hazard_shape()] for the fitted counterpart
#'
#' @export
#'
#' @examples
#' data(guinea_pig)
#' ttt <- bd_ttt_plot(guinea_pig$time)
#' attr(ttt, "shape")
bd_ttt_plot <- function(time, status = NULL, add = FALSE,
                        col = "steelblue", lwd = 2,
                        main = "Scaled total time on test",
                        xlab = "i / n", ylab = expression(phi(i/n)), ...) {

  if (inherits(time, "betadanish")) {
    status <- time$data$status
    time   <- time$data$time
  }
  time <- as.numeric(time)
  if (!is.null(status)) {
    keep <- as.numeric(status) == 1
    if (sum(keep) < 5L)
      stop("At least five uncensored observations are needed.", call. = FALSE)
    if (any(!keep))
      warning(sprintf(paste0("%d censored observation(s) dropped: the TTT ",
                             "transform is defined for complete samples."),
                      sum(!keep)), call. = FALSE)
    time <- time[keep]
  }
  time <- sort(time[is.finite(time) & time > 0])
  n <- length(time)
  if (n < 5L) stop("At least five positive times are needed.", call. = FALSE)

  i   <- seq_len(n)
  cs  <- cumsum(time)
  phi <- (cs + (n - i) * time) / cs[n]
  u   <- i / n

  shape <- .bd_ttt_shape(u, phi)

  if (!isTRUE(add)) {
    plot(c(0, u), c(0, phi), type = "l", col = col, lwd = lwd,
                   xlim = c(0, 1), ylim = c(0, 1),
                   main = main, xlab = xlab, ylab = ylab, ...)
    graphics::abline(0, 1, col = "grey50", lty = 2)
    graphics::legend("bottomright", bty = "n",
                     legend = c("TTT transform", "diagonal (constant hazard)"),
                     col = c(col, "grey50"), lty = c(1, 2), lwd = c(lwd, 1))
    graphics::mtext(paste("suggested hazard:", shape), side = 3, line = 0.2,
                    cex = 0.85)
  } else {
    graphics::lines(c(0, u), c(0, phi), col = col, lwd = lwd)
  }

  out <- data.frame(i_n = u, phi = phi)
  attr(out, "shape") <- shape
  invisible(out)
}

#' Classify a TTT Curve Against the Diagonal
#' @noRd
.bd_ttt_shape <- function(u, phi, tol = 0.02) {
  d <- phi - u
  above <- d > tol
  below <- d < -tol
  if (!any(above) && !any(below)) return("constant (exponential)")
  if (!any(below)) return("increasing")
  if (!any(above)) return("decreasing")
  ## Mixed: the side it starts on decides between bathtub and unimodal.
  first <- if (which(above)[1] < which(below)[1]) "above" else "below"
  if (identical(first, "below")) "bathtub" else "unimodal (upside-down bathtub)"
}

#' Plot a Profile Likelihood
#'
#' Draws the profile log-likelihood produced by [bd_profile_ci()], with the
#' critical threshold and the resulting interval marked.
#'
#' @param x An object of class `"bd_profile"`.
#' @param col,lwd Colour and line width for the profile curve.
#' @param main,xlab,ylab Labels. `NULL` uses sensible defaults.
#' @param ... Further graphical parameters.
#'
#' @return Invisibly returns `x`.
#'
#' @details
#' The horizontal line sits at \eqn{\ell_{\max} - \chi^2_{1,\alpha}/2}. Every
#' parameter value whose profile lies above it is inside the interval.
#'
#' When the curve does not fall back below that line before the right-hand edge
#' of the grid, there is no finite upper bound and the plot says so. Widening
#' the grid will not produce one; it will only confirm the flatness. That is the
#' situation the underlying dissertation records for the tail index on the
#' breaking-stress data.
#'
#' @seealso [bd_profile_ci()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(120, a = 1, b = 3, c = 2, k = 0.5, seed = 3)
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                       submodel = TRUE, n_starts = 1)
#' p <- bd_profile_ci(fit, "b", n_grid = 15L)
#' bd_profile_plot(p)
#' }
bd_profile_plot <- function(x, col = "steelblue", lwd = 2,
                            main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (!inherits(x, "bd_profile"))
    stop("'x' must be a bd_profile object from bd_profile_ci().", call. = FALSE)

  ok <- is.finite(x$profile)
  if (sum(ok) < 3L)
    stop("Too few finite profile values to plot.", call. = FALSE)

  thr <- x$logLik_max - x$critical
  if (is.null(main)) main <- paste("Profile log-likelihood for", x$parameter)
  if (is.null(xlab)) xlab <- x$parameter
  if (is.null(ylab)) ylab <- "profile log-likelihood"

  plot(x$grid[ok], x$profile[ok], type = "l", col = col, lwd = lwd,
                 main = main, xlab = xlab, ylab = ylab, ...)
  graphics::abline(h = thr, col = "red", lty = 2)
  graphics::abline(v = x$estimate, col = "grey40", lty = 3)

  if (is.finite(x$lower)) graphics::abline(v = x$lower, col = "red", lty = 3)
  if (is.finite(x$upper)) graphics::abline(v = x$upper, col = "red", lty = 3)

  lab <- sprintf("%s%% interval: [%s, %s]", format(100 * x$level),
                 signif(x$lower, 4),
                 if (is.infinite(x$upper)) "Inf" else signif(x$upper, 4))
  graphics::mtext(lab, side = 3, line = 0.2, cex = 0.85)

  if (isTRUE(x$open_above))
    graphics::legend("bottomright", bty = "n", text.col = "red3",
                     legend = "no finite upper bound: report a lower bound")

  invisible(x)
}

#' Diagnostic Plots for a Bayesian Fit
#'
#' Trace and posterior density plots for each parameter of a
#' [bayes_betadanish()] fit.
#'
#' @param x An object of class `"bd_bayes"`.
#' @param which Optional character vector of parameters to show. Defaults to
#'   all of them.
#' @param type `"both"` (default) for trace and density side by side, or
#'   `"trace"` or `"density"` alone.
#' @param col Colour for the traces and densities.
#' @param ... Further graphical parameters.
#'
#' @return Invisibly returns `x`.
#'
#' @details
#' Read the traces first. A well-mixed chain looks like noise around a stable
#' level, with no drift and no long excursions. Visible trend means the burn-in
#' was too short; a chain that sticks at one value for many iterations means
#' the proposal is too wide and almost every move is being rejected, which is
#' worth fixing with `tune` before interpreting anything.
#'
#' The graphical parameters are restored on exit, so the function leaves the
#' device as it found it.
#'
#' @seealso [bayes_betadanish()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("MCMCpack", quietly = TRUE) &&
#'     requireNamespace("coda", quietly = TRUE)) {
#'   dat <- simulate_bd_data(80, a = 1, b = 3, c = 2, k = 0.5, seed = 6)
#'   bfit <- bayes_betadanish(dat$time, dat$status, submodel = TRUE,
#'                            burnin = 200, mcmc = 800, seed = 1)
#'   plot(bfit)
#' }
#' }
plot.bd_bayes <- function(x, which = NULL, type = c("both", "trace", "density"),
                          col = "steelblue", ...) {
  type <- match.arg(type)
  dr <- as.matrix(x$draws)
  if (!is.matrix(dr) || !nrow(dr))
    stop("The fit contains no posterior draws.", call. = FALSE)

  nm <- colnames(dr)
  if (!is.null(which)) {
    miss <- setdiff(which, nm)
    if (length(miss))
      stop("Not in the posterior: ", paste(miss, collapse = ", "),
           ". Available: ", paste(nm, collapse = ", "), call. = FALSE)
    dr <- dr[, which, drop = FALSE]
    nm <- which
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  ncol_panel <- if (identical(type, "both")) 2L else 1L
  graphics::par(mfrow = c(length(nm), ncol_panel),
                mar = c(4, 4, 2, 1))

  it <- seq_len(nrow(dr))
  for (p in nm) {
    v <- dr[, p]
    if (type %in% c("both", "trace")) {
      plot(it, v, type = "l", col = col,
                     main = paste("Trace:", p),
                     xlab = "iteration", ylab = p, ...)
      graphics::abline(h = mean(v), col = "red", lty = 2)
    }
    if (type %in% c("both", "density")) {
      ## Plot the coordinates rather than the density object: relying on S3
      ## dispatch for plot.density is fragile inside a package namespace, and
      ## falling through to plot.default on a list is an error, not a warning.
      d <- stats::density(v)
      plot(d$x, d$y, type = "l", col = col, lwd = 2,
           main = paste("Posterior:", p), xlab = p, ylab = "density", ...)
      graphics::abline(v = mean(v), col = "red", lty = 2)
      if (!is.null(x$HPD) && p %in% rownames(x$HPD))
        graphics::abline(v = x$HPD[p, ], col = "grey40", lty = 3)
    }
  }
  invisible(x)
}
