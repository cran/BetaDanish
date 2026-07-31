#' Plot Diagnostics for Beta-Danish Fit
#'
#' Generates diagnostic plots for a fitted Beta-Danish model, including survival,
#' hazard, density, PP, and QQ plots.
#'
#' @param x A fitted `betadanish` object.
#' @param type Character string specifying the plot type: `"survival"`, `"hazard"`,
#'   `"density"`, `"pp"`, `"qq"`, or `"all"`.
#'
#' @param ... Additional arguments passed to the base `plot` function.
#' @return Invisibly returns the input betadanish object. Called mainly for its side effect of producing diagnostic plots.
#' @export
plot.betadanish <- function(x, type = c("survival", "hazard", "density", "pp", "qq", "all"), ...) {
  type <- match.arg(type)

  time <- x$data$time
  status <- x$data$status
  time_max <- max(time, na.rm = TRUE)
  seq_t <- seq(1e-5, time_max, length.out = 200)

  pars <- as.list(x$coefficients)
  if (x$submodel) pars$a <- 1

  # Helper to reset par if "all" is chosen
  old_par <- graphics::par(no.readonly = TRUE)
  if (type == "all") {
    graphics::par(mfrow = c(2, 3))
    on.exit(graphics::par(old_par))
  }

  plot_types <- if (type == "all") c("survival", "hazard", "density", "pp", "qq") else type

  for (pt in plot_types) {
    if (pt == "survival") {
      km <- survival::survfit(survival::Surv(time, status) ~ 1)
      graphics::plot(km, conf.int = FALSE, xlab = "Time", ylab = "Survival Probability",
                     main = "Survival Fit: Beta-Danish vs KM", ...)
      fit_probs <- pbetadanish(seq_t, pars$a, pars$b, pars$c, pars$k, lower.tail = FALSE)
      graphics::lines(seq_t, fit_probs, col = "red", lwd = 2)
      graphics::legend("topright", legend = c("Kaplan-Meier", "Beta-Danish"),
                       col = c("black", "red"), lty = 1, lwd = c(1, 2), bty = "n")

    } else if (pt == "hazard") {
      haz_vals <- hbetadanish(seq_t, pars$a, pars$b, pars$c, pars$k)
      graphics::plot(seq_t, haz_vals, type = "l", col = "blue", lwd = 2,
                     xlab = "Time", ylab = "Hazard Rate",
                     main = "Estimated Hazard Function", ...)
      graphics::grid()

    } else if (pt == "density") {
      # Density plot makes most sense for uncensored data, but we plot it anyway
      graphics::hist(time[status == 1], probability = TRUE, breaks = 20,
                     col = "lightgray", border = "white",
                     xlab = "Time", main = "Density Fit (Events Only)", ...)
      pdf_vals <- dbetadanish(seq_t, pars$a, pars$b, pars$c, pars$k)
      graphics::lines(seq_t, pdf_vals, col = "darkgreen", lwd = 2)

    } else if (pt == "pp") {
      # PP Plot (Empirical CDF vs Fitted CDF)
      km <- survival::survfit(survival::Surv(time, status) ~ 1)
      emp_cdf <- 1 - km$surv
      fit_cdf <- pbetadanish(km$time, pars$a, pars$b, pars$c, pars$k)
      graphics::plot(fit_cdf, emp_cdf, xlab = "Theoretical Probabilities",
                     ylab = "Empirical Probabilities", main = "P-P Plot", pch = 16, col = "blue", ...)
      graphics::abline(0, 1, col = "red", lwd = 2)

    } else if (pt == "qq") {
      # QQ Plot (Empirical Quantiles vs Fitted Quantiles)
      km <- survival::survfit(survival::Surv(time, status) ~ 1)
      emp_cdf <- 1 - km$surv
      # Avoid 1.0 for quantile function
      emp_cdf[emp_cdf >= 1] <- 0.9999
      fit_q <- qbetadanish(emp_cdf, pars$a, pars$b, pars$c, pars$k)
      graphics::plot(fit_q, km$time, xlab = "Theoretical Quantiles",
                     ylab = "Empirical Quantiles", main = "Q-Q Plot", pch = 16, col = "blue", ...)
      graphics::abline(0, 1, col = "red", lwd = 2)
    }
  }
  invisible(x)
}
