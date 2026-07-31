#' Bladder Cancer Remission Times
#'
#' Remission times (in months) for 128 bladder cancer patients. This is a
#' complete (uncensored) sample widely used in lifetime distribution literature
#' to demonstrate decreasing or right-skewed hazard rates.
#'
#' @format A data frame with 128 rows and 2 columns:
#' \describe{
#'   \item{time}{Remission time in months}
#'   \item{status}{Event indicator (1 = event occurred)}
#' }
#' @source Lee, E. T., & Wang, J. W. (2003). Statistical Methods for Survival Data Analysis (3rd ed.). Wiley.
#' @examples
#' data(remission)
#' \donttest{
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = remission)
#' summary(fit)
#' }
"remission"

#' Breaking Stress of Carbon Fibres
#'
#' Breaking stress (in GPa) of 100 carbon fibre specimens. This dataset exhibits
#' a unimodal (increasing-then-decreasing) hazard pattern that classical
#' distributions like the Weibull cannot adequately capture.
#'
#' @format A data frame with 100 rows and 2 columns:
#' \describe{
#'   \item{time}{Breaking stress in GPa}
#'   \item{status}{Event indicator (1 = event occurred)}
#' }
#' @source Nichols, M. D., & Padgett, W. J. (2006). A bootstrap control chart for Weibull percentiles. Quality and Reliability Engineering International, 22(2), 141-151.
#' @examples
#' data(carbon_fibres)
#' \donttest{
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = carbon_fibres)
#' }
"carbon_fibres"

#' Bone Marrow Transplant Survival
#'
#' Survival times for 91 patients with refractory acute lymphoblastic leukemia
#' who received either an allogeneic or autologous bone marrow transplant.
#' This dataset includes right-censoring and a treatment covariate, making it
#' ideal for demonstrating cure-rate models and AFT regression.
#'
#' @format A data frame with 91 rows and 3 columns:
#' \describe{
#'   \item{time}{Survival time in days}
#'   \item{status}{Event indicator (1 = death/relapse, 0 = censored)}
#'   \item{group}{Treatment group (0 = Allogeneic, 1 = Autologous)}
#' }
#' @source Klein, J. P., & Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and Truncated Data (2nd ed.). Springer.
#' @examples
#' data(transplant)
#' \donttest{
#' # Fit a model with a covariate
#' fit <- fit_bd_aft(survival::Surv(time, status) ~ group, data = transplant)
#' }
"transplant"

#' Aarset Device Failure Times
#'
#' Times to failure of 50 devices, exhibiting a classic bathtub-shaped hazard rate.
#' This is a standard benchmark dataset in reliability engineering.
#'
#' @format A data frame with 50 rows and 2 columns:
#' \describe{
#'   \item{time}{Failure time}
#'   \item{status}{Event indicator (1 = event occurred)}
#' }
#' @source Aarset, M. V. (1987). How to Identify a Bathtub Hazard Rate. IEEE Transactions on Reliability, R-36(1), 106-108.
#' @examples
#' data(aarset)
#' \donttest{
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = aarset)
#' plot(fit, type = "hazard")
#' }
"aarset"

#' Acute Myelogenous Leukemia Survival
#'
#' Survival times (in weeks) for 23 patients with acute myelogenous leukemia.
#' A classic, small dataset perfect for fast testing of censored data workflows.
#'
#' @format A data frame with 23 rows and 3 columns:
#' \describe{
#'   \item{time}{Survival time in weeks}
#'   \item{status}{Event indicator (1 = event, 0 = censored)}
#'   \item{group}{Treatment group (Maintained vs Non-maintained)}
#' }
#' @source Miller, R. G. (1997). Survival Analysis. Wiley.
#' @examples
#' data(leukemia)
#' \donttest{
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = leukemia)
#' }
"leukemia"

#' Malignant Melanoma Survival After Surgery
#'
#' Survival times for 205 patients with malignant melanoma after surgery.
#' This rich clinical dataset includes multiple covariates and heavy censoring.
#'
#' @format A data frame with 205 rows and 7 columns:
#' \describe{
#'   \item{time}{Survival time in days}
#'   \item{status}{Event indicator (1 = died from melanoma, 0 = alive, 2 = other death)}
#'   \item{thickness}{Tumor thickness in mm}
#'   \item{sex}{Patient sex (1 = male, 0 = female)}
#'   \item{age}{Patient age in years}
#'   \item{ulcer}{Ulceration indicator (1 = present, 0 = absent)}
#'   \item{year}{Year of operation}
#' }
#' @source Andersen, P. K., Borgan, O., Gill, R. D., & Keiding, N. (1993). Statistical Models Based on Counting Processes. Springer.
#' @examples
#' data(melanoma)
#' \donttest{
#' # Treat status 1 as event, others as censored
#' melanoma$event <- ifelse(melanoma$status == 1, 1, 0)
#' fit <- fit_betadanish(survival::Surv(time, event) ~ age + thickness, data = melanoma)
#' }
"melanoma"
#' Guinea Pig Survival Times
#'
#' Survival times, in days, of 72 guinea pigs injected with virulent
#' tubercle bacilli, taken from the principal regimen of the study of
#' Bjerkedal (1960). A complete sample: every animal was observed to death,
#' so `status` is 1 throughout.
#'
#' @format A data frame with 72 rows and 2 columns:
#' \describe{
#'   \item{time}{Survival time in days}
#'   \item{status}{Event indicator, 1 for all observations (no censoring)}
#' }
#'
#' @details
#' These data are the reference case for a well-identified four-parameter
#' fit. The scaled total-time-on-test transform is unimodal and the upper
#' tail is determinate rather than heavy, and the maximum-likelihood fit of
#' the full Beta-Danish model attains a finite interior optimum: every
#' estimate lies strictly inside the parameter space, with
#' \eqn{\hat b = 3.64} (Wald standard error 1.20), so that
#' \eqn{(\hat b - 1)/\mathrm{SE} = 2.20} places \eqn{b} more than two
#' standard errors clear of the \eqn{b = 1} identifiability ridge.
#'
#' By contrast, on `remission` and `carbon_fibres` the four-parameter fit
#' sits on the flat \eqn{(a, c)} direction of the likelihood with
#' \eqn{\hat a < 1} and is only weakly identified. Use this dataset when
#' you want to see the parent model behaving well; see the Identifiability
#' section of [fit_betadanish()] for what to watch for elsewhere.
#'
#' @source Bjerkedal, T. (1960). Acquisition of resistance in guinea pigs
#'   infected with different doses of virulent tubercle bacilli.
#'   *American Journal of Epidemiology*, 72(1), 130-148.
#'   \doi{10.1093/oxfordjournals.aje.a120129}
#'
#' @examples
#' data(guinea_pig)
#' summary(guinea_pig$time)
#' \donttest{
#' fit_full <- fit_betadanish(survival::Surv(time, status) ~ 1,
#'                            data = guinea_pig)
#' fit_sub  <- fit_betadanish(survival::Surv(time, status) ~ 1,
#'                            data = guinea_pig, submodel = TRUE)
#' compare_models(fit_full, fit_sub)
#' }
"guinea_pig"
