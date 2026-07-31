## Internal helpers for the CSV pipeline. Every model fit goes through
## .bd_try(), so a single failure is recorded and the run continues instead of
## aborting and losing the analyses that did succeed.

#' Run an Expression, Recording Rather Than Propagating Failure
#'
#' Warnings are captured rather than discarded. The identifiability diagnostics
#' from [fit_betadanish()] arrive as warnings, and they are among the most
#' important things a user of this pipeline needs to see, so they are collected
#' and surfaced on the result instead of being muffled into silence.
#'
#' @noRd
.bd_try <- function(expr, label, store) {
  ws <- character(0)
  res <- withCallingHandlers(
    tryCatch(expr,
             error = function(e) structure(list(message = conditionMessage(e)),
                                           class = "bd_failed")),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(ws))
    store$warnings <- c(store$warnings,
                        stats::setNames(ws, rep(label, length(ws))))
  if (inherits(res, "bd_failed")) {
    store$failures <- c(store$failures,
                        stats::setNames(res$message, label))
    return(NULL)
  }
  res
}

#' Parameter Estimate Table for a betadanish Fit
#' @noRd
.bd_est_table <- function(fit, label) {
  est <- fit$coefficients
  se  <- sqrt(pmax(diag(fit$vcov), 0))
  z   <- stats::qnorm(0.975)
  data.frame(
    model      = label,
    parameter  = names(est),
    estimate   = as.numeric(est),
    std_error  = as.numeric(se),
    lower_95   = as.numeric(est - z * se),
    upper_95   = as.numeric(est + z * se),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

#' Information-Criteria Table Across Fits
#' @noRd
.bd_ic_table <- function(fits) {
  rows <- lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    if (is.null(f) || is.null(f$logLik)) return(NULL)
    data.frame(model   = nm,
               n       = .bd_or(f$nobs, NA_integer_),
               npar    = .bd_or(f$npar, length(f$coefficients)),
               logLik  = as.numeric(f$logLik),
               AIC     = tryCatch(stats::AIC(f), error = function(e) NA_real_),
               BIC     = tryCatch(stats::BIC(f), error = function(e) NA_real_),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Build a Survival Formula from Retained Covariate Names
#' @noRd
.bd_surv_formula <- function(covs) {
  rhs <- if (length(covs)) paste(covs, collapse = " + ") else "1"
  stats::as.formula(paste0("survival::Surv(time, status) ~ ", rhs),
                    env = parent.frame())
}

#' Write One PNG, Closing the Device Whatever Happens
#'
#' The device is closed before the file is tested, so a half-written PNG is
#' never reported as a success. Graphical parameters need no saving here: they
#' are per-device, and the device is discarded.
#'
#' @noRd
.bd_png <- function(path, expr, width = 1600, height = 1200, res = 180) {
  opened <- FALSE
  ok <- tryCatch({
    grDevices::png(path, width = width, height = height, res = res)
    opened <- TRUE
    force(expr)
    TRUE
  }, error = function(e) FALSE)
  if (opened) try(grDevices::dev.off(), silent = TRUE)
  if (isTRUE(ok) && file.exists(path)) path else character(0)
}

#' Analyse Survival Data Straight from a CSV File
#'
#' A single entry point that reads a delimited or Excel file, fits the
#' requested Beta-Danish model or models, assembles the results into tidy
#' tables, and optionally writes those tables and the diagnostic figures to a
#' directory.
#'
#' @param file Path to a `.csv`, `.txt`, `.tsv`, `.xls` or `.xlsx` file.
#' @param time_col,status_col,cause_col Column names, passed to
#'   [read_survival_data()]. `NULL` means guess.
#' @param covariates Covariate columns to use. `NULL` means none for
#'   `analysis = "univariate"`, and every non-response column for the
#'   regression analyses.
#' @param analysis Which analysis to run: `"univariate"` (default), `"aft"`,
#'   `"cure"` or `"competing"`.
#' @param model For `analysis = "univariate"`, whether to fit the
#'   four-parameter Beta-Danish (`"BD"`), the three-parameter Exponentiated
#'   Danish submodel (`"ED"`), or `"both"` (default) and compare them.
#' @param compare Logical; benchmark against standard lifetime distributions
#'   via [compare_distributions()]. Requires the `flexsurv` package; skipped
#'   with a recorded note if it is absent. Default `TRUE`.
#' @param bayes Logical; also run [bayes_betadanish()]. Requires `MCMCpack`.
#'   Default `FALSE`, since it is slow.
#' @param cure_formula Right-hand-side formula for the cure fraction when
#'   `analysis = "cure"`, for example `~ group`. Defaults to the retained
#'   covariates.
#' @param cure_type `"mixture"` (default) or `"promotion"`.
#' @param output_dir Directory to write results to. `NULL` (default) writes
#'   nothing. The directory is created if it does not exist.
#' @param n_starts Number of random starting points for each fit.
#' @param seed Optional integer seed, for reproducible multi-start fitting.
#' @param quiet Logical; suppress progress messages.
#'
#' @return An object of class `"bd_analysis"`: a list with `call`, `analysis`,
#'   `data`, `data_report`, `fits`, `tables`, `extras`, `warnings`, `failures`,
#'   `output_dir` and `files`. It has `print`, `summary` and `plot` methods.
#'
#' @details
#' Nothing is written to disk unless `output_dir` is supplied. When it is, the
#' function writes one CSV per result table and a PNG per diagnostic figure,
#' and records the paths in the `files` component.
#'
#' Each fit is attempted independently. If one fails, the message is recorded
#' in `failures` and the remaining analyses still run, so a difficult
#' four-parameter fit does not cost you the submodel results alongside it.
#'
#' Warnings raised during fitting are captured in `warnings` rather than
#' printed as they occur. This is where the identifiability diagnostics appear,
#' so it is worth reading: a four-parameter fit that converges cleanly but sits
#' on the `b = 1` ridge will say so there.
#'
#' For `analysis = "univariate"` with `model = "both"`, a likelihood ratio test
#' of the submodel is included. Read it alongside those diagnostics: see the
#' Identifiability section of [fit_betadanish()].
#'
#' @seealso [read_survival_data()], [bd_csv_template()], [fit_betadanish()],
#'   [fit_bd_aft()], [fit_bd_cure()], [fit_bd_competing()]
#'
#' @export
#'
#' @examples
#' f <- system.file("extdata", "censored_sample.csv", package = "BetaDanish")
#'
#' # Submodel only, no benchmarking: fast enough for an example
#' res <- bd_analyze_csv(f, analysis = "univariate", model = "ED",
#'                       compare = FALSE, n_starts = 3, seed = 1, quiet = TRUE)
#' res
#' res$tables$information_criteria
#'
#' \donttest{
#' # Both models, with tables and figures written to a temporary directory
#' out <- file.path(tempdir(), "bd_results")
#' res2 <- bd_analyze_csv(f, model = "both", compare = FALSE,
#'                        output_dir = out, n_starts = 5, seed = 1)
#' basename(res2$files)
#'
#' # AFT regression, covariates taken from the file
#' g <- system.file("extdata", "covariate_sample.csv", package = "BetaDanish")
#' res3 <- bd_analyze_csv(g, analysis = "aft", covariates = c("age", "thickness"),
#'                        n_starts = 5, seed = 1)
#' summary(res3)
#' }
bd_analyze_csv <- function(file,
                           time_col = NULL, status_col = NULL,
                           covariates = NULL, cause_col = NULL,
                           analysis = c("univariate", "aft", "cure", "competing"),
                           model = c("both", "BD", "ED"),
                           compare = TRUE,
                           bayes = FALSE,
                           cure_formula = NULL,
                           cure_type = c("mixture", "promotion"),
                           output_dir = NULL,
                           n_starts = 10,
                           seed = NULL,
                           quiet = FALSE) {

  analysis  <- match.arg(analysis)
  model     <- match.arg(model)
  cure_type <- match.arg(cure_type)
  if (!is.null(seed)) set.seed(seed)
  say <- function(...) if (!isTRUE(quiet)) message(...)

  store <- new.env(parent = emptyenv())
  store$failures <- character(0)
  store$warnings <- character(0)

  ## ---- read -----------------------------------------------------------------
  covar_arg <- covariates
  if (is.null(covar_arg) && analysis %in% c("aft", "cure")) covar_arg <- "all"

  dat <- read_survival_data(file, time_col = time_col, status_col = status_col,
                            covar_cols = covar_arg, cause_col = cause_col,
                            quiet = quiet)
  rep  <- attr(dat, "bd_data_report")
  covs <- rep$covariates

  say("Read ", rep$rows_kept, " row(s); ", rep$n_events, " event(s); ",
      sprintf("%.1f%%", 100 * rep$censoring_prop), " censored.")

  out <- list(call = match.call(), analysis = analysis, model = model,
              data = dat, data_report = rep,
              fits = list(), tables = list(), extras = list(),
              failures = character(0), output_dir = NULL, files = character(0))

  ## ---- fit ------------------------------------------------------------------
  if (analysis == "univariate") {
    fml <- .bd_surv_formula(character(0))

    if (model %in% c("BD", "both")) {
      say("Fitting the four-parameter Beta-Danish model.")
      out$fits$BD <- .bd_try(
        fit_betadanish(fml, data = dat, submodel = FALSE, n_starts = n_starts),
        "fit_betadanish (BD)", store)
    }
    if (model %in% c("ED", "both")) {
      say("Fitting the three-parameter Exponentiated Danish submodel.")
      out$fits$ED <- .bd_try(
        fit_betadanish(fml, data = dat, submodel = TRUE, n_starts = n_starts),
        "fit_betadanish (ED)", store)
    }

    if (!is.null(out$fits$BD) && !is.null(out$fits$ED)) {
      lrt <- .bd_try(compare_models(out$fits$BD, out$fits$ED),
                     "compare_models", store)
      if (!is.null(lrt)) out$tables$likelihood_ratio_test <- as.data.frame(lrt)
    }

  } else if (analysis == "aft") {
    if (!length(covs))
      say("No covariates retained; the AFT fit reduces to an intercept-only ",
          "model, equivalent to the ED submodel with k = exp(intercept).")
    fml <- .bd_surv_formula(covs)
    say("Fitting the AFT model: ", deparse(fml[[3]]))
    out$fits$AFT <- .bd_try(
      fit_bd_aft(fml, data = dat, n_starts = n_starts), "fit_bd_aft", store)

  } else if (analysis == "cure") {
    if (is.null(cure_formula))
      cure_formula <- stats::as.formula(
        paste0("~ ", if (length(covs)) paste(covs, collapse = " + ") else "1"))
    say("Fitting the ", cure_type, " cure model; cure fraction ~ ",
        deparse(cure_formula[[2]]))
    out$fits$CURE <- .bd_try(
      fit_bd_cure(formula_aft  = .bd_surv_formula(character(0)),
                  formula_cure = cure_formula,
                  data = dat, type = cure_type, n_starts = n_starts),
      "fit_bd_cure", store)

  } else if (analysis == "competing") {
    if (!"cause" %in% names(dat))
      stop("analysis = \"competing\" needs a cause column. Pass cause_col, ",
           "and code it 0 for censored with 1, 2, ... for the causes.",
           call. = FALSE)
    say("Fitting cause-specific models for ",
        length(setdiff(unique(dat$cause), 0)), " cause(s).")
    out$fits$CR <- .bd_try(
      fit_bd_competing(time = dat$time, cause = dat$cause, n_starts = n_starts),
      "fit_bd_competing", store)

    if (!is.null(out$fits$CR)) {
      if (requireNamespace("cmprsk", quietly = TRUE)) {
        cc <- .bd_try(cif_compare(out$fits$CR, plot = FALSE), "cif_compare", store)
        if (!is.null(cc)) out$extras$cif <- cc
      } else {
        store$failures <- c(store$failures,
                            c(cif_compare = "the 'cmprsk' package is not installed"))
      }
    }
  }

  ## ---- estimate and IC tables ----------------------------------------------
  uni <- Filter(function(f) inherits(f, "betadanish"), out$fits)
  if (length(uni)) {
    out$tables$estimates <- do.call(
      rbind, Map(.bd_est_table, uni, names(uni)))
    out$tables$information_criteria <- .bd_ic_table(uni)

    gof <- lapply(uni, function(f) .bd_try(gof_betadanish(f), "gof_betadanish", store))
    gof <- Filter(Negate(is.null), gof)
    if (length(gof)) {
      out$tables$goodness_of_fit <- do.call(rbind, Map(function(g, nm)
        data.frame(model = nm, as.list(g$IC), KS = g$KS_Statistic,
                   row.names = NULL, check.names = FALSE),
        gof, names(gof)))
    }
  }

  for (nm in intersect(names(out$fits), c("AFT", "CURE"))) {
    f <- out$fits[[nm]]
    if (is.null(f)) next
    s <- .bd_try(summary(f), paste0("summary (", nm, ")"), store)
    if (!is.null(s)) out$extras[[tolower(nm)]] <- s
  }

  ## ---- optional benchmarking ------------------------------------------------
  primary <- if (!is.null(out$fits$BD)) out$fits$BD else out$fits$ED
  if (isTRUE(compare) && inherits(primary, "betadanish")) {
    if (requireNamespace("flexsurv", quietly = TRUE)) {
      say("Benchmarking against standard lifetime distributions.")
      cmp <- .bd_try(compare_distributions(primary), "compare_distributions", store)
      if (!is.null(cmp)) out$tables$distribution_benchmark <- as.data.frame(cmp)
    } else {
      store$failures <- c(store$failures,
                          c(compare_distributions = "the 'flexsurv' package is not installed"))
    }
  }

  ## ---- optional Bayesian run ------------------------------------------------
  if (isTRUE(bayes)) {
    if (requireNamespace("MCMCpack", quietly = TRUE)) {
      say("Running the Bayesian sampler; this takes a while.")
      out$extras$bayes <- .bd_try(
        bayes_betadanish(dat$time, dat$status,
                         submodel = !identical(model, "BD"), seed = seed),
        "bayes_betadanish", store)
    } else {
      store$failures <- c(store$failures,
                          c(bayes_betadanish = "the 'MCMCpack' package is not installed"))
    }
  }

  ## ---- data report table ----------------------------------------------------
  out$tables$data_report <- data.frame(
    file           = rep$file,
    time_col       = .bd_or(rep$time_col, NA_character_),
    status_col     = .bd_or(rep$status_col, NA_character_),
    rows_read      = rep$rows_read,
    rows_kept      = rep$rows_kept,
    n_events       = rep$n_events,
    censoring_prop = rep$censoring_prop,
    covariates     = paste(rep$covariates, collapse = "; "),
    grid_step      = .bd_or(rep$grid_step, NA_real_),
    row.names = NULL, stringsAsFactors = FALSE)

  out$failures <- store$failures
  out$warnings <- store$warnings
  class(out) <- "bd_analysis"

  ## ---- optional output ------------------------------------------------------
  if (!is.null(output_dir)) out <- .bd_write_outputs(out, output_dir, quiet)

  if (length(out$warnings) && !isTRUE(quiet))
    message("Completed with ", length(out$warnings),
            " model warning(s); see $warnings. Identifiability diagnostics ",
            "appear here.")
  if (length(out$failures) && !isTRUE(quiet))
    message("Completed with ", length(out$failures),
            " recorded failure(s); see $failures.")

  out
}

#' Write Tables and Figures for a bd_analysis Object
#' @noRd
.bd_write_outputs <- function(obj, output_dir, quiet = FALSE) {
  say <- function(...) if (!isTRUE(quiet)) message(...)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir))
    stop("Could not create '", output_dir, "'.", call. = FALSE)

  written <- character(0)

  for (nm in names(obj$tables)) {
    tb <- obj$tables[[nm]]
    if (is.null(tb) || !nrow(as.data.frame(tb))) next
    p <- file.path(output_dir, paste0(nm, ".csv"))
    ok <- tryCatch({ utils::write.csv(tb, p, row.names = FALSE); TRUE },
                   error = function(e) FALSE)
    if (ok) written <- c(written, p)
  }

  prim <- if (!is.null(obj$fits$BD)) obj$fits$BD else obj$fits$ED
  if (inherits(prim, "betadanish")) {
    for (ty in c("survival", "hazard", "density", "pp", "qq")) {
      written <- c(written,
                   .bd_png(file.path(output_dir, paste0(ty, ".png")),
                           plot(prim, type = ty)))
    }
  }
  for (nm in c("AFT", "CURE")) {
    f <- obj$fits[[nm]]
    if (is.null(f)) next
    written <- c(written,
                 .bd_png(file.path(output_dir, paste0("coxsnell_", tolower(nm), ".png")),
                         plot(f)))
  }
  if (!is.null(obj$fits$CR) && requireNamespace("cmprsk", quietly = TRUE)) {
    written <- c(written,
                 .bd_png(file.path(output_dir, "cif.png"),
                         cif_compare(obj$fits$CR, plot = TRUE)))
  }

  say("Wrote ", length(written), " file(s) to ", normalizePath(output_dir))
  obj$output_dir <- normalizePath(output_dir, mustWork = FALSE)
  obj$files      <- written
  obj
}

#' @param x A `"bd_analysis"` object.
#' @param ... Ignored.
#' @rdname bd_analyze_csv
#' @export
print.bd_analysis <- function(x, ...) {
  cat("\nBetaDanish CSV analysis\n")
  cat("-----------------------\n")
  cat("Analysis:      ", x$analysis, "\n", sep = "")
  cat("Observations:  ", x$data_report$rows_kept,
      "  (events: ", x$data_report$n_events,
      ", censored: ", sprintf("%.1f%%", 100 * x$data_report$censoring_prop),
      ")\n", sep = "")
  if (length(x$data_report$covariates))
    cat("Covariates:    ", paste(x$data_report$covariates, collapse = ", "),
        "\n", sep = "")
  cat("Models fitted: ",
      if (length(x$fits)) paste(names(Filter(Negate(is.null), x$fits)),
                                collapse = ", ") else "none", "\n", sep = "")

  ic <- x$tables$information_criteria
  if (!is.null(ic)) { cat("\n"); print(ic, row.names = FALSE) }

  if (length(x$files))
    cat("\nOutput:        ", length(x$files), " file(s) in ", x$output_dir,
        "\n", sep = "")
  if (length(x$warnings))
    cat("\nWarnings:      ", length(x$warnings),
        " (see $warnings)\n", sep = "")
  if (length(x$failures))
    cat("\nFailures:      ", length(x$failures),
        " (see $failures)\n", sep = "")
  cat("\n")
  invisible(x)
}

#' @param object A `"bd_analysis"` object.
#' @rdname bd_analyze_csv
#' @export
summary.bd_analysis <- function(object, ...) {
  structure(object, class = c("summary.bd_analysis", "bd_analysis"))
}

#' @rdname bd_analyze_csv
#' @export
print.summary.bd_analysis <- function(x, ...) {
  y <- x; class(y) <- "bd_analysis"
  print(y)

  for (nm in names(x$tables)) {
    if (identical(nm, "information_criteria")) next
    tb <- x$tables[[nm]]
    if (is.null(tb)) next
    cat("-- ", gsub("_", " ", nm), " --\n", sep = "")
    print(tb, row.names = FALSE)
    cat("\n")
  }
  for (nm in names(x$extras)) {
    cat("-- ", nm, " --\n", sep = "")
    print(x$extras[[nm]])
    cat("\n")
  }
  if (length(x$warnings)) {
    cat("-- model warnings --\n")
    for (i in seq_along(x$warnings))
      cat("  ", names(x$warnings)[i], ": ", x$warnings[i], "\n", sep = "")
    cat("\n")
  }
  if (length(x$failures)) {
    cat("-- failures --\n")
    for (i in seq_along(x$failures))
      cat("  ", names(x$failures)[i], ": ", x$failures[i], "\n", sep = "")
    cat("\n")
  }
  invisible(x)
}

#' @param type Plot type passed to [plot.betadanish()].
#' @rdname bd_analyze_csv
#' @export
plot.bd_analysis <- function(x, type = "survival", ...) {
  prim <- if (!is.null(x$fits$BD)) x$fits$BD else x$fits$ED
  if (inherits(prim, "betadanish")) return(plot(prim, type = type, ...))
  for (nm in c("AFT", "CURE")) {
    if (!is.null(x$fits[[nm]])) return(plot(x$fits[[nm]], ...))
  }
  if (!is.null(x$fits$CR)) return(cif_compare(x$fits$CR, plot = TRUE))
  warning("Nothing to plot: no model was fitted successfully.", call. = FALSE)
  invisible(x)
}
