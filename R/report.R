#' Comprehensive Beta-Danish Analysis Pipeline
#'
#' Runs a complete end-to-end analysis: reads data, fits the 4-parameter and
#' 3-parameter models, compares them, benchmarks against standard distributions,
#' and generates diagnostic plots.
#'
#' @param file Path to the CSV or Excel file containing the data.
#' @param time_col Name of the time column.
#'
#' @param status_col Name of the status column (optional).
#' @return Invisibly returns a list containing the fitted full model and submodel objects. The function is mainly called for its side effects of printing an analysis report and producing diagnostic plots.
#' @export
analyze_betadanish <- function(file, time_col, status_col = NULL) {
  cat("======================================================\n")
  cat("       BETA-DANISH COMPREHENSIVE ANALYSIS REPORT      \n")
  cat("======================================================\n\n")

  # 1. Read Data
  cat("[1] Reading Data...\n")
  dat <- read_survival_data(file, time_col, status_col)
  cat("    Successfully loaded", nrow(dat), "observations.\n\n")

  # 2. Fit Models
  cat("[2] Fitting Models...\n")
  form <- survival::Surv(time, status) ~ 1

  cat("    Fitting Full 4-Parameter Model...\n")
  fit_full <- fit_betadanish(form, data = dat, submodel = FALSE)

  cat("    Fitting 3-Parameter Submodel (a=1)...\n")
  fit_sub <- fit_betadanish(form, data = dat, submodel = TRUE)

  # 3. Print Summaries
  cat("\n[3] Model Summaries:\n")
  print(summary(fit_full))

  # 4. Compare Nesting
  cat("\n[4] Nested Model Comparison:\n")
  compare_models(fit_full, fit_sub)

  # 5. Compare Distributions
  if (requireNamespace("flexsurv", quietly = TRUE)) {
    cat("\n[5] Benchmarking Against Standard Distributions:\n")
    comp_tab <- compare_distributions(fit_full)
    print(comp_tab, row.names = FALSE)
  } else {
    cat("\n[5] Benchmarking skipped (install 'flexsurv' to enable).\n")
  }

  # 6. Plotting
  cat("\n[6] Generating Diagnostic Plots...\n")
  plot(fit_full, type = "all")

  cat("\n======================================================\n")
  cat("                  ANALYSIS COMPLETE                   \n")
  cat("======================================================\n")

  invisible(list(full_model = fit_full, sub_model = fit_sub))
}
