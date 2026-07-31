#' @keywords internal
#' @aliases BetaDanish-package
#'
#' @details
#' The `BetaDanish` package provides a comprehensive suite of tools for survival
#' and reliability analysis using the Beta-Danish distribution. It includes core
#' mathematical functions (density, distribution, quantile, random generation,
#' and hazard), robust Maximum Likelihood Estimation (MLE) for both complete and
#' right-censored data, and advanced modeling capabilities including Accelerated
#' Failure Time (AFT) regression, cure-rate models, and competing risks analysis.
#'
#' The three-parameter Exponentiated Danish (ED) submodel is introduced in the
#' published article (Ahmad & Danish, 2025). The full four-parameter Beta-Danish
#' distribution is developed in the accompanying doctoral thesis (Ahmad, 2026).
#'
#' @references
#' Ahmad, B., & Danish, M. Y. (2025). Development and characterization of a
#' flexible three-parameter lifetime distribution: theoretical properties and
#' real-world applications. *Journal of Applied Mathematics, Statistics and
#' Informatics*, 21(1). \doi{10.2478/jamsi-2025-0010}
#'
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan. Supervised by
#' Dr. Muhammad Yameen Danish.
#'
#' @importFrom maxLik maxLik
#' @importFrom survival Surv survfit is.Surv
#' @importFrom stats AIC BIC logLik optim pbeta qbeta rbeta runif pnorm pchisq model.extract model.frame model.response delete.response terms na.omit na.pass median var rnorm rbinom rpois rexp uniroot
#' @importFrom graphics lines plot legend grid abline hist
#' @importFrom utils head tail packageVersion
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "time", "status", "group", "cause", "cif", "type", "Model", "S",
    "Param", "Bias", "Estimate", "Lower", "Upper", "post_cure"
  ))
}
