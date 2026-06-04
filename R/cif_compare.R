#' Compare Model-Based CIF to the Aalen-Johansen Estimator
#'
#' Computes the nonparametric Aalen-Johansen CIF (via \pkg{cmprsk}) for each
#' competing-risks cause, overlays it on the fitted Beta-Danish CIF, and
#' returns the Aalen-Johansen times/estimates, the fitted CIF values on a
#' common grid, and Gray's CIF-equality test where applicable.
#'
#' @param fit A fitted object of class \code{"bd_competing"}.
#' @param tmax Optional upper time for evaluation; default the 95th
#'   percentile of observed times.
#' @param n_grid Number of time points on the evaluation grid (default 160).
#' @param plot Logical; if \code{TRUE} (default) draws a panel of overlays.
#'
#' @return A list with elements \code{tgrid}, \code{cif_fit} (data frame
#'   long format), \code{cif_aj} (data frame long format) and optionally
#'   \code{gray_test}.
#'
#' @details Requires \pkg{cmprsk} (Suggests).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' T1 <- rbetadanish(200, 1.2, 1.5, 1.0, 0.4)
#' T2 <- rbetadanish(200, 1.0, 2.0, 1.0, 0.2)
#' C  <- stats::rexp(200, 0.05)
#' time  <- pmin(T1, T2, C)
#' cause <- ifelse(time == C, 0L, ifelse(T1 <= T2, 1L, 2L))
#' fit <- fit_bd_competing(time = time, cause = cause)
#' cif_compare(fit)   # requires cmprsk to be installed
#' }
#'
#' @export
cif_compare <- function(fit, tmax = NULL, n_grid = 160, plot = TRUE) {
  if (!inherits(fit, "bd_competing"))
    stop("fit must be a bd_competing object.")
  if (!requireNamespace("cmprsk", quietly = TRUE))
    stop("cif_compare() requires the 'cmprsk' package. Install with install.packages('cmprsk').")
  time  <- fit$data$time
  cause <- fit$data$cause
  if (is.null(tmax))
    tmax <- as.numeric(stats::quantile(time, 0.95, na.rm = TRUE))
  tgrid <- seq(0, tmax, length.out = n_grid)
  aj <- cmprsk::cuminc(ftime = time, fstatus = cause, cencode = 0)
  cause_curves <- aj[!grepl("Tests", names(aj))]

  cif_fit_long <- list()
  cif_aj_long  <- list()
  for (j_idx in seq_along(fit$causes)) {
    j <- fit$causes[j_idx]
    fitted_cif <- cif_betadanish(fit, tgrid, cause_idx = j)
    cif_fit_long[[j_idx]] <- data.frame(
      time = tgrid, cif = fitted_cif,
      cause = paste0("Cause_", j), source = "Beta-Danish fitted")
    aj_match <- grep(paste0(" ", j, "$"), names(cause_curves))
    if (length(aj_match) == 0L) aj_match <- j_idx
    aj_curve <- cause_curves[[aj_match[1]]]
    cif_aj_long[[j_idx]] <- data.frame(
      time = aj_curve$time, cif = aj_curve$est,
      cause = paste0("Cause_", j), source = "Aalen-Johansen")
  }
  cif_fit_df <- do.call(rbind, cif_fit_long)
  cif_aj_df  <- do.call(rbind, cif_aj_long)

  if (plot) {
    m <- length(fit$causes)
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    graphics::par(mfrow = c(1, m), mar = c(4.5, 4.5, 3, 1))
    for (j_idx in seq_along(fit$causes)) {
      j <- fit$causes[j_idx]
      fit_sub <- cif_fit_df[cif_fit_df$cause == paste0("Cause_", j), ]
      aj_sub  <- cif_aj_df[cif_aj_df$cause == paste0("Cause_", j), ]
      graphics::plot(fit_sub$time, fit_sub$cif, type = "l", lwd = 2,
                     col = "red", xlim = c(0, tmax), ylim = c(0, 1),
                     xlab = "Time", ylab = paste0("CIF, cause ", j),
                     main = paste0("Cause ", j))
      graphics::lines(aj_sub$time, aj_sub$cif, type = "s",
                      lwd = 2, lty = 2, col = "black")
      graphics::legend("bottomright",
                       legend = c("Beta-Danish fitted", "Aalen-Johansen"),
                       col = c("red", "black"), lwd = 2,
                       lty = c(1, 2), bty = "n")
    }
  }

  out <- list(tgrid = tgrid, cif_fit = cif_fit_df, cif_aj = cif_aj_df)
  if (!is.null(aj$Tests)) out$gray_test <- aj$Tests
  out
}
