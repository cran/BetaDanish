## Finite-sample simulation studies.
##
## Each runner draws `n_sim` samples at every sample size, refits, and reports
## bias, RMSE, the ratio of mean estimated standard error to empirical standard
## deviation, and Wald coverage. That last pair is the informative one: a ratio
## far below one means the standard errors understate the sampling variability,
## and coverage well below the nominal level follows from it.
##
## Fits that fail or that the degeneracy guard rejects are counted rather than
## silently dropped, because a high non-convergence rate is itself a finding
## about the parameter region being studied.

#' Simulate Competing Risks Data
#'
#' @param n Sample size.
#' @param pars List of parameter vectors, one per cause, each a named vector
#'   with `b`, `c`, `k` and optionally `a`. Defaults to two causes.
#' @param gammas Numeric vector of accelerated failure time coefficients, one
#'   per cause, acting on a single binary covariate. `NULL` for no covariate.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param seed Optional integer seed.
#'
#' @return A data frame with `time`, `cause` (0 for censored) and, when
#'   `gammas` is supplied, the covariate `x`.
#'
#' @seealso [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' d <- simulate_bd_competing_data(200, seed = 1)
#' table(d$cause)
simulate_bd_competing_data <- function(n,
                                       pars = list(c(a = 1, b = 2.5, c = 1.8, k = 0.04),
                                                   c(a = 1, b = 3.5, c = 1.2, k = 0.02)),
                                       gammas = NULL,
                                       censor_rate = 0.15,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- length(pars)
  if (m < 2L) stop("At least two causes are needed.", call. = FALSE)

  x <- if (is.null(gammas)) NULL else stats::rbinom(n, 1, 0.5)
  if (!is.null(gammas) && length(gammas) != m)
    stop("'gammas' must have one entry per cause.", call. = FALSE)

  Tm <- matrix(Inf, n, m)
  for (j in seq_len(m)) {
    p  <- pars[[j]]
    a  <- if ("a" %in% names(p)) p[["a"]] else 1
    t0 <- rbetadanish(n, a, p[["b"]], p[["c"]], p[["k"]])
    Tm[, j] <- if (is.null(gammas)) t0 else t0 * exp(x * gammas[j])
  }

  tmin  <- apply(Tm, 1, min)
  cause <- max.col(-Tm, ties.method = "first")

  cmax <- stats::quantile(tmin, 1 - censor_rate, na.rm = TRUE) * 2
  cens <- stats::runif(n, 0, cmax)
  cause[cens < tmin] <- 0L

  out <- data.frame(time = pmin(tmin, cens), cause = as.integer(cause))
  if (!is.null(x)) out$x <- x
  out
}

#' Summarise One Column of Simulation Replicates
#' @noRd
.bd_sim_summary <- function(est, se, truth, level = 0.95) {
  ok <- is.finite(est)
  z  <- stats::qnorm(1 - (1 - level) / 2)
  cov <- if (any(is.finite(se[ok]))) {
    mean(abs(est[ok] - truth) <= z * se[ok], na.rm = TRUE)
  } else NA_real_
  c(truth      = truth,
    mean       = mean(est[ok]),
    bias       = mean(est[ok]) - truth,
    rel_bias   = (mean(est[ok]) - truth) / truth,
    rmse       = sqrt(mean((est[ok] - truth)^2)),
    emp_sd     = stats::sd(est[ok]),
    mean_se    = mean(se[ok], na.rm = TRUE),
    se_ratio   = mean(se[ok], na.rm = TRUE) / stats::sd(est[ok]),
    coverage   = cov)
}

#' Finite-Sample Simulation Study for the Univariate Model
#'
#' Draws samples from a known Beta-Danish or Exponentiated Danish distribution,
#' refits each one, and reports the finite-sample behaviour of the maximum
#' likelihood estimator.
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param truth Named numeric vector of true parameters. Must contain `b`, `c`
#'   and `k`, and `a` unless `submodel` is `TRUE`.
#' @param submodel Logical; fit the three-parameter ED submodel.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param n_starts Random starts per fit, on top of the deterministic grid.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"` with a `results` data frame,
#'   one row per sample size and parameter.
#'
#' @details
#' Read `se_ratio` and `coverage` together. A ratio near one with coverage near
#' the nominal level means the asymptotic standard errors are trustworthy at
#' that sample size. A ratio well below one means they are optimistic, and
#' coverage will fall short accordingly.
#'
#' `n_fail` counts replicates where no admissible optimum was found. A high
#' rate is not a nuisance to be suppressed: it says the parameter region being
#' studied is hard to estimate at that sample size.
#'
#' @seealso [bd_simulation_cure()], [bd_simulation_competing()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Deliberately tiny so the example runs in seconds. A real study would use
#' # n_sim in the hundreds; see the vignette.
#' s <- bd_simulation_study(n = 60, n_sim = 3,
#'                          truth = c(b = 3, c = 2, k = 0.5),
#'                          submodel = TRUE, n_starts = 1,
#'                          seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_study <- function(n = c(50, 100, 200), n_sim = 200,
                                truth = c(a = 1.5, b = 3, c = 2, k = 0.5),
                                submodel = FALSE, censor_rate = 0,
                                n_starts = 3, level = 0.95,
                                seed = NULL, quiet = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  need <- if (submodel) c("b", "c", "k") else c("a", "b", "c", "k")
  if (!all(need %in% names(truth)))
    stop("'truth' must name: ", paste(need, collapse = ", "), call. = FALSE)
  a_true <- if (submodel) 1 else truth[["a"]]

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(need),
                  dimnames = list(NULL, need))
    se  <- est
    for (r in seq_len(n_sim)) {
      t0 <- rbetadanish(nn, a_true, truth[["b"]], truth[["c"]], truth[["k"]])
      status <- rep(1L, nn)
      if (censor_rate > 0) {
        cmax <- stats::quantile(t0, 1 - censor_rate) * 2
        cc   <- stats::runif(nn, 0, cmax)
        status <- as.integer(t0 <= cc)
        t0 <- pmin(t0, cc)
      }
      dat <- data.frame(time = t0, status = status)
      fit <- tryCatch(suppressWarnings(
        fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
                       submodel = submodel, n_starts = n_starts,
                       check_identifiability = FALSE)),
        error = function(e) NULL)
      if (is.null(fit)) next
      est[r, ] <- fit$coefficients[need]
      se[r, ]  <- sqrt(pmax(diag(fit$vcov)[need], 0))
    }

    n_ok <- sum(stats::complete.cases(est))
    for (p in need) {
      s <- .bd_sim_summary(est[, p], se[, p], truth[[p]], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = p, as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, truth = truth,
                            submodel = submodel, censor_rate = censor_rate,
                            level = level),
              type = "univariate")
  class(out) <- "bd_simulation"
  out
}

#' Generate Mixture-Cure Data
#'
#' A cured subject never fails, so it is censored at whatever censoring time it
#' draws. Susceptible subjects fail according to the ED distribution. Written
#' here rather than reusing `simulate_bd_cure_data()`, whose interface requires
#' explicit regression coefficients and design matrices that an intercept-only
#' study does not need.
#'
#' @noRd
.bd_sim_cure_data <- function(n, b, c, k, cure_fraction, censor_rate) {
  cured <- stats::rbinom(n, 1, cure_fraction) == 1L
  t <- rep(Inf, n)
  ns <- sum(!cured)
  if (ns) t[!cured] <- rbetadanish(ns, 1, b, c, k)

  fin <- t[is.finite(t)]
  cmax <- if (length(fin)) {
    as.numeric(stats::quantile(fin, max(1 - censor_rate, 0.05))) * 3
  } else 1
  cens <- stats::runif(n, 0, cmax)

  data.frame(time = pmin(t, cens), status = as.integer(t <= cens))
}

#' Finite-Sample Simulation Study for the Cure Model
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param truth Named vector with `b`, `c` and `k` for the susceptible
#'   population.
#' @param cure_fraction True proportion of long-term survivors.
#' @param censor_rate Approximate proportion of censoring among susceptibles.
#' @param type `"mixture"` or `"promotion"`.
#' @param n_starts Random starts per fit.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"`.
#'
#' @seealso [bd_simulation_study()], [fit_bd_cure()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' s <- bd_simulation_cure(n = 80, n_sim = 2, n_starts = 1,
#'                         seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_cure <- function(n = c(100, 200), n_sim = 100,
                               truth = c(b = 2, c = 1.5, k = 0.5),
                               cure_fraction = 0.3, censor_rate = 0.2,
                               type = c("mixture", "promotion"),
                               n_starts = 3, level = 0.95,
                               seed = NULL, quiet = FALSE) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)
  need <- c("b", "c")

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(need), dimnames = list(NULL, need))
    se  <- est
    for (r in seq_len(n_sim)) {
      dat <- .bd_sim_cure_data(nn, truth[["b"]], truth[["c"]], truth[["k"]],
                               cure_fraction, censor_rate)
      fit <- tryCatch(suppressWarnings(
        fit_bd_cure(formula_aft = survival::Surv(time, status) ~ 1,
                    formula_cure = ~ 1, data = dat, type = type,
                    n_starts = n_starts)),
        error = function(e) NULL)
      if (is.null(fit)) next
      cf <- fit$coefficients
      v  <- sqrt(pmax(diag(fit$vcov), 0))
      for (p in need) {
        nm <- paste0("log_", p)
        if (nm %in% names(cf)) {
          est[r, p] <- exp(cf[[nm]])
          se[r, p]  <- exp(cf[[nm]]) * v[which(names(cf) == nm)]
        }
      }
    }
    n_ok <- sum(stats::complete.cases(est))
    for (p in need) {
      s <- .bd_sim_summary(est[, p], se[, p], truth[[p]], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = p, as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, truth = truth,
                            cure_fraction = cure_fraction, type = type,
                            level = level),
              type = paste0("cure (", type, ")"))
  class(out) <- "bd_simulation"
  out
}

#' Finite-Sample Simulation Study for Competing Risks
#'
#' @param n Vector of sample sizes.
#' @param n_sim Number of replicates at each sample size.
#' @param pars List of true parameter vectors, one per cause.
#' @param gammas Optional accelerated failure time coefficients, one per cause.
#' @param censor_rate Approximate proportion of right-censored observations.
#' @param submodel Logical; fit Exponentiated Danish cause-specific kernels.
#' @param n_starts Starting points per fit.
#' @param level Nominal coverage level.
#' @param seed Optional integer seed.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_simulation"`, with `parameter` naming the
#'   cause and quantity, for example `c1:b`.
#'
#' @seealso [bd_simulation_study()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' s <- bd_simulation_competing(n = 120, n_sim = 2, n_starts = 1,
#'                              seed = 1, quiet = TRUE)
#' s
#' }
bd_simulation_competing <- function(n = c(200, 400), n_sim = 100,
                                    pars = list(c(b = 2.5, c = 1.8, k = 0.04),
                                                c(b = 3.5, c = 1.2, k = 0.02)),
                                    gammas = NULL, censor_rate = 0.15,
                                    submodel = TRUE, n_starts = 2,
                                    level = 0.95, seed = NULL, quiet = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  m <- length(pars)
  base <- c("b", "c", "k")
  targets <- unlist(lapply(seq_len(m), function(j) paste0("c", j, ":", base)))
  truths  <- unlist(lapply(seq_len(m), function(j) unname(pars[[j]][base])))
  if (!is.null(gammas)) {
    targets <- c(targets, paste0("c", seq_len(m), ":gamma_x"))
    truths  <- c(truths, gammas)
  }

  rows <- list()
  for (nn in n) {
    say("n = ", nn, ": ", n_sim, " replicate(s)")
    est <- matrix(NA_real_, n_sim, length(targets),
                  dimnames = list(NULL, targets))
    se  <- est
    for (r in seq_len(n_sim)) {
      d <- tryCatch(simulate_bd_competing_data(nn, pars = pars,
                                               gammas = gammas,
                                               censor_rate = censor_rate),
                    error = function(e) NULL)
      if (is.null(d) || length(unique(d$cause[d$cause > 0])) < m) next
      fit <- tryCatch(suppressWarnings(
        fit_bd_competing(d$time, d$cause,
                         covariates = if (is.null(gammas)) NULL else ~ x,
                         data = d, submodel = submodel, n_starts = n_starts)),
        error = function(e) NULL)
      if (is.null(fit)) next
      cm <- fit$coefficients; sm <- fit$se
      for (j in seq_len(m)) {
        for (p in base) {
          est[r, paste0("c", j, ":", p)] <- cm[j, p]
          se[r,  paste0("c", j, ":", p)] <- sm[j, p]
        }
        if (!is.null(gammas) && "gamma_x" %in% colnames(cm)) {
          est[r, paste0("c", j, ":gamma_x")] <- cm[j, "gamma_x"]
          se[r,  paste0("c", j, ":gamma_x")] <- sm[j, "gamma_x"]
        }
      }
    }
    n_ok <- sum(stats::complete.cases(est))
    for (i in seq_along(targets)) {
      s <- .bd_sim_summary(est[, targets[i]], se[, targets[i]], truths[i], level)
      rows[[length(rows) + 1L]] <- data.frame(
        n = nn, parameter = targets[i], as.list(s),
        n_ok = n_ok, n_fail = n_sim - n_ok,
        row.names = NULL, check.names = FALSE)
    }
  }

  out <- list(results = do.call(rbind, rows),
              design = list(n = n, n_sim = n_sim, pars = pars,
                            gammas = gammas, submodel = submodel,
                            level = level),
              type = "competing risks")
  class(out) <- "bd_simulation"
  out
}

#' @param x A `"bd_simulation"` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @rdname bd_simulation_study
#' @export
print.bd_simulation <- function(x, digits = 4, ...) {
  cat("\nBeta-Danish simulation study: ", x$type, "\n", sep = "")
  cat("  replicates: ", x$design$n_sim,
      "   sample sizes: ", paste(x$design$n, collapse = ", "), "\n", sep = "")
  cat("  nominal coverage: ", format(100 * x$design$level), "%\n\n", sep = "")

  res <- x$results
  show <- c("n", "parameter", "truth", "mean", "bias", "rmse",
            "se_ratio", "coverage", "n_fail")
  print(format(res[, intersect(show, names(res))], digits = digits),
        row.names = FALSE)

  if (any(res$n_fail > 0))
    cat("\n  Some replicates found no admissible optimum. A high count means\n",
        "  the parameter region is hard to estimate at that sample size.\n",
        sep = "")
  bad <- res$se_ratio < 0.9 & is.finite(res$se_ratio)
  if (any(bad))
    cat("\n  se_ratio below 0.9 for ", sum(bad), " row(s): the asymptotic\n",
        "  standard errors understate the sampling variability there, and\n",
        "  coverage falls short of nominal as a result.\n", sep = "")
  cat("\n")
  invisible(x)
}

#' @param object A `"bd_simulation"` object.
#' @rdname bd_simulation_study
#' @export
summary.bd_simulation <- function(object, ...) object$results
