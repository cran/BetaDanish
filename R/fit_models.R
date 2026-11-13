#' Fit the Beta-Danish Distribution to Survival Data
#'
#' Fits the Beta-Danish distribution by maximum likelihood. Complete and
#' right-censored samples are both supported, via a `survival::Surv` response.
#'
#' @param formula A formula whose left-hand side is a `Surv` object. Use
#'   `~ 1` for a model without covariates.
#' @param data A data frame containing the variables in `formula`.
#' @param submodel Logical; if `TRUE`, fits the three-parameter Exponentiated
#'   Danish (ED) submodel by fixing `a = 1`.
#' @param n_starts Integer; number of random starting points added to the
#'   thesis starts and the deterministic grid. Default 10.
#' @param method Character; optimisation method passed to `maxLik::maxLik`.
#' @param check_identifiability Logical; if `TRUE` (default), issue warnings
#'   when the fit lands in a region where the parameters are weakly identified.
#' @param penalty Non-negative ridge penalty applied on the log-parameter
#'   scale. `0` (default) is ordinary maximum likelihood. A small positive
#'   value stabilises the fit when the likelihood is nearly flat, at the cost
#'   of bias toward `penalty_center`.
#' @param penalty_center Numeric vector on the log-parameter scale toward which
#'   the penalty shrinks, or `NULL` (default) to shrink toward the best
#'   unpenalised start.
#' @param grouped Logical; if `TRUE`, use the grouped (interval) likelihood
#'   appropriate to times recorded on a coarse grid. Default `FALSE`.
#' @param delta Recording increment for `grouped = TRUE`. `NULL` (default)
#'   infers it from the spacing of the observed times.
#' @param box Numeric vector `c(lower, upper)` confining the shape parameters
#'   `a`, `b` and `c`. Default `c(1e-4, 1e4)`, the box used in the thesis.
#'   The rate `k` gets the same box widened to follow the unit of time; see
#'   the Parameter box section.
#'
#' @return An object of S3 class `"betadanish"` with components including
#'   `coefficients`, `logLik`, `vcov`, `npar`, `nobs`, `AIC`, `BIC`,
#'   `convergence`, `diagnostics`, `at_boundary` (logical), `edge_par` (the
#'   parameters at the box edge), `flag` (`"ok"` or a description),
#'   `bounds` (the natural-scale box) and, for the four-parameter model,
#'   `submodel_logLik`, the log-likelihood of the Stage-1 ED fit.
#'
#' @details
#' Optimisation is carried out on log-transformed parameters; estimates and the
#' variance-covariance matrix are returned on the natural scale, the latter via
#' the delta method.
#'
#' @section Estimator:
#' The search follows the two-stage estimator of the thesis.
#'
#' * **Stage 1** fits the Exponentiated Danish submodel from four starts,
#'   \eqn{(b, c) \in \{(1,1), (2,2), (0.5,3), (4,0.8)\}} with
#'   \eqn{k = 1/\mathrm{median}(t)}.
#' * **Stage 2** (four-parameter model only) starts from seven points along
#'   the \eqn{(a, c)} ridge through the Stage-1 optimum,
#'   \eqn{(f, b_{ED}, c_{ED}/f, k_{ED})} for
#'   \eqn{f \in \{1/30, 1/8, 1/3, 1, 3, 8, 30\}}. The start \eqn{f = 1} is
#'   the ED optimum itself, so the Beta-Danish log-likelihood can never fall
#'   below the ED one: nesting holds by construction.
#'
#' The deterministic grid and the `n_starts` random starts of earlier versions
#' are kept as further candidates, drawn from the random-number stream in the
#' same order as before. For a complete, ungrouped sample the analytic score
#' (thesis Appendix B.1) is used as the gradient; censored and grouped
#' likelihoods are differentiated numerically. No start value depends on the
#' parameters that generated the data.
#'
#' @section Parameter box:
#' The shape parameters are confined to `box`, \eqn{[10^{-4}, 10^{4}]} by
#' default. The rate `k` is measured in the unit of time, so its bound covers
#' both \eqn{[10^{-4}, 10^{4}]} and \eqn{k \cdot \mathrm{median}(t) \in}
#' `box`; rescaling the times cannot push an interior fit onto the boundary.
#'
#' A fit within 2\% of the box edge (on the log scale) is flagged
#' `at_boundary = TRUE`. Its log-likelihood is a **supremum**, not an attained
#' maximum, so `AIC`, `BIC`, the standard errors and the goodness-of-fit
#' statistics of [gof_betadanish()] are all `NA`. An interior fit with
#' \eqn{\hat a < 0.01} or \eqn{\hat c > 100} is flagged as weakly identified
#' on the \eqn{(a, c)} ridge but is still a proper maximum and keeps its
#' information criteria.
#'
#' @section Grouped data:
#' Survival times are often recorded on a grid -- whole days, whole months --
#' and the point-density likelihood is not appropriate for them. It treats a
#' rounded value as an exact observation, which can overstate the information
#' in the sample and make uncertainty estimates too small. With `grouped = TRUE` an
#' event recorded at \eqn{t} contributes
#' \eqn{\log\{F(t + \delta/2) - F(t - \delta/2)\}} instead of
#' \eqn{\log f(t)}; censored observations are unchanged. The cell probability
#' falls back to \eqn{f(t)\delta} only where the difference of two nearly
#' equal distribution values has cancelled to zero.
#'
#' `read_survival_data()` reports an inferred `grid_step` for exactly this
#' purpose, and this function warns when the times look grid-recorded but
#' `grouped = FALSE`.
#'
#' The grouped likelihood is a harder surface to optimise than the
#' point-density one. On the same data and the same starting grid, the
#' optimiser can stop where the observed information is not positive
#' definite, and the delta-method variances then come back non-positive.
#' Check `fit$diagnostics$vcov_singular` before trusting a standard error
#' from a grouped fit; leaving `check_identifiability = TRUE` will warn
#' about it automatically.
#'
#' @section Penalised fitting:
#' `penalty > 0` adds \eqn{\lambda \sum (\theta - \mu)^2} on the
#' log-parameter scale. This is worth reaching for when the likelihood is flat
#' along the \eqn{(a, c)} direction and the unpenalised optimiser wanders, but
#' it is a deliberate bias: the estimates are shrunk toward `penalty_center`.
#'
#' The reported `logLik` is always the **unpenalised** log-likelihood evaluated
#' at the penalised estimate, so that AIC, BIC and likelihood ratio tests
#' remain comparable across fits. The objective actually maximised is stored
#' separately as `penalised_logLik`. Treat the degrees of freedom as nominal:
#' shrinkage reduces the effective number of parameters, so information
#' criteria are conservative under penalisation.
#'
#' @section Identifiability:
#' The four-parameter model is not uniformly well identified, and a converged
#' fit is not by itself evidence that it is. Two regions warrant care.
#'
#' * **The \eqn{b = 1} ridge.** At \eqn{b = 1} the beta generator collapses and
#'   the model is non-identifiable. A fit with \eqn{\hat b} within about two
#'   standard errors of one lies close to that ridge; the likelihood is nearly
#'   flat along it, so the individual estimates carry little information even
#'   though the fitted survival curve may look excellent.
#' * **Lower-tail \eqn{(a, c)} confounding.** Near the lower tail, \eqn{a} and
#'   \eqn{c} enter almost exclusively through the product \eqn{ca}, so the
#'   expected Fisher information is close to singular in that direction. A
#'   fitted correlation between \eqn{\hat a} and \eqn{\hat c} above about 0.95
#'   in absolute value indicates that only the product is being estimated.
#'
#' In either case the ED submodel (`submodel = TRUE`) is usually the honest
#' report, and a likelihood ratio test via [compare_models()] will normally
#' fail to reject it. Set `check_identifiability = FALSE` to silence the
#' warnings once you have satisfied yourself that they are understood.
#'
#' @seealso [compare_models()], [gof_betadanish()], [compare_distributions()]
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim_time   <- rbetadanish(150, a = 1.5, b = 3, c = 2, k = 0.5)
#' sim_status <- rbinom(150, 1, 0.85)
#' dat <- data.frame(time = sim_time, status = sim_status)
#'
#' fit <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat)
#' summary(fit)
#'
#' fit_sub <- fit_betadanish(survival::Surv(time, status) ~ 1, data = dat,
#'                           submodel = TRUE)
#' compare_models(fit, fit_sub)
#' }
fit_betadanish <- function(formula, data, submodel = FALSE, n_starts = 10,
                           method = "BFGS", check_identifiability = TRUE,
                           penalty = 0, penalty_center = NULL,
                           grouped = FALSE, delta = NULL,
                           box = c(1e-4, 1e4)) {

  surv_data <- extract_surv_data(formula, data)
  time   <- surv_data$time
  status <- surv_data$status

  if (!is.numeric(penalty) || length(penalty) != 1L || is.na(penalty) || penalty < 0)
    stop("'penalty' must be a single non-negative number.", call. = FALSE)

  grid_step <- .bd_grid_step(time)
  if (isTRUE(grouped)) {
    if (is.null(delta)) delta <- grid_step
    if (!is.finite(delta) || delta <= 0)
      stop("grouped = TRUE needs a recording increment. The spacing of the ",
           "times did not imply one, so supply 'delta' explicitly.",
           call. = FALSE)
  } else if (is.finite(grid_step) && isTRUE(check_identifiability)) {
    warning(sprintf(paste0("The times look recorded on a grid of %s. The ",
                           "point-density likelihood treats them as exact, ",
                           "which understates the standard errors. Consider ",
                           "grouped = TRUE."), format(grid_step)),
            call. = FALSE)
  }

  ## ---- parameter box (thesis: [1e-4, 1e4]) ---------------------------------
  bounds <- .bd_param_box(time, box)
  lo_eta <- stats::setNames(log(bounds$lower), paste0("log_", names(bounds$lower)))
  hi_eta <- stats::setNames(log(bounds$upper), paste0("log_", names(bounds$upper)))
  box_tol <- 0.02
  clamp <- function(pars) {
    nm <- names(pars)
    stats::setNames(pmin(pmax(as.numeric(pars), lo_eta[nm]), hi_eta[nm]), nm)
  }

  ev <- status == 1

  ## Unpenalised log-likelihood, exact or grouped, on the natural scale.
  loglik_nat <- function(a_par, b_par, c_par, k_par) {
    lp <- if (isTRUE(grouped)) {
      .bd_log_cell(time[ev], delta, a_par, b_par, c_par, k_par)
    } else {
      dbetadanish(time[ev], a_par, b_par, c_par, k_par, log = TRUE)
    }
    v <- sum(lp)
    if (any(!ev))
      v <- v + sum(pbetadanish(time[!ev], a_par, b_par, c_par, k_par,
                               lower.tail = FALSE, log.p = TRUE))
    v
  }

  ## Log-scale parameters -> natural scale, inside the box. A vector without
  ## log_a is the ED submodel (a = 1); this serves Stage 1 of the full model.
  nat_of <- function(pars) {
    e <- exp(clamp(pars))
    c(a = if ("log_a" %in% names(e)) e[["log_a"]] else 1,
      b = e[["log_b"]], c = e[["log_c"]], k = e[["log_k"]])
  }

  loglik_fun <- function(pars) {
    th <- nat_of(pars)
    suppressWarnings(loglik_nat(th[["a"]], th[["b"]], th[["c"]], th[["k"]]))
  }

  pen_center <- penalty_center

  ll_fun <- function(pars) {
    loglik <- loglik_fun(pars)
    if (!is.finite(loglik)) return(-1e10)
    if (penalty > 0 && !is.null(pen_center)) {
      mu <- pen_center[names(pars)]
      if (!anyNA(mu)) loglik <- loglik - penalty * sum((clamp(pars) - mu)^2)
    }
    if (!is.finite(loglik)) return(-1e10)
    loglik
  }

  ## Analytic score (thesis Appendix B.1): complete, ungrouped samples only.
  ## Outside the box the objective is flat, so the gradient is zero there.
  use_score <- !isTRUE(grouped) && all(ev)
  score_fun <- function(pars, penalised) {
    nm <- names(pars)
    if (!is.finite(loglik_fun(pars))) return(stats::setNames(rep(0, length(pars)), nm))
    cl  <- clamp(pars)
    sub <- !("log_a" %in% nm)
    ord <- if (sub) c("log_b", "log_c", "log_k") else c("log_a", "log_b", "log_c", "log_k")
    g <- stats::setNames(.bd_score_eta(time, exp(as.numeric(cl[ord])), submodel = sub), ord)
    g <- g[nm]
    if (penalised && penalty > 0 && !is.null(pen_center)) {
      mu <- pen_center[nm]
      if (!anyNA(mu)) g <- g - 2 * penalty * (cl - mu)
    }
    inside <- as.numeric(cl == as.numeric(pars))
    stats::setNames(as.numeric(g) * inside, nm)
  }
  grad_unpen <- if (use_score) function(p) score_fun(p, FALSE) else NULL
  grad_pen   <- if (use_score) function(p) score_fun(p, TRUE)  else NULL
  ll_unpen_finite <- function(p) {
    v <- loglik_fun(p); if (is.finite(v)) v else -1e10
  }

  ## ---- starting values -------------------------------------------------------
  avg_t <- mean(time[status == 1], na.rm = TRUE)
  if (is.na(avg_t) || avg_t <= 0) avg_t <- mean(time, na.rm = TRUE)
  k_base <- 1 / avg_t

  ## The deterministic grid and the random starts of version 0.3.0, drawn in
  ## the same order, so the candidate set only ever grows.
  start_list <- .bd_default_starts(submodel, k_base)
  for (i in seq_len(n_starts)) {
    core <- c(log_b = log(stats::runif(1, 0.5, 5)),
              log_c = log(stats::runif(1, 0.5, 5)),
              log_k = log(k_base * stats::runif(1, 0.5, 2)))
    start_list[[length(start_list) + 1L]] <- if (submodel) core else
      c(log_a = log(stats::runif(1, 0.5, 5)), core)
  }

  ## Admissibility: finite, and no runaway log-likelihood. The shape bound of
  ## earlier versions is replaced by the box itself.
  accept <- .bd_make_accept(n = length(time), lim_shape = Inf, lim_scale = Inf)

  ## Thesis Stage 1 (ED) and, for the full model, Stage 2 ridge starts.
  med_t <- stats::median(time[ev])
  if (!is.finite(med_t) || med_t <= 0) med_t <- stats::median(time)
  k0 <- 1 / max(med_t, 1e-8)
  ed_starts <- .bd_ed_starts(k0, bounds$lower, bounds$upper)

  ed_stage <- NULL
  if (submodel) {
    thesis_starts <- ed_starts
  } else {
    ed_stage <- optim_multistart(ll_unpen_finite, ed_starts, method = method,
                                 accept = accept, grad = grad_unpen)
    thesis_starts <- list()
    if (!is.null(ed_stage)) {
      e3 <- exp(clamp(ed_stage$estimate))
      thesis_starts <- .bd_ridge_starts(
        c(b = e3[["log_b"]], c = e3[["log_c"]], k = e3[["log_k"]]),
        bounds$lower, bounds$upper)
    }
  }
  start_list <- c(thesis_starts, start_list)

  ## With a penalty and no explicit centre, shrink toward the unpenalised
  ## optimum rather than toward an arbitrary point.
  if (penalty > 0 && is.null(pen_center)) {
    pilot <- optim_multistart(ll_unpen_finite, start_list, method = method,
                              accept = accept, grad = grad_unpen)
    if (!is.null(pilot)) pen_center <- clamp(pilot$estimate)
  }

  fit <- optim_multistart(ll_fun, start_list, method = method, accept = accept,
                          grad = grad_pen)
  if (is.null(fit))
    stop("No admissible optimum was found. Every start either failed or ",
         "returned a runaway log-likelihood. Try method = \"NM\", or fit the ",
         "ED submodel with submodel = TRUE.", call. = FALSE)

  ## Polish the selected optimum to a tight tolerance (the thesis engine used
  ## reltol 1e-10). The default relative tolerance of maxLik, about 1.5e-8,
  ## can leave the log-likelihood some 1e-5 short of the maximum on a sample
  ## of a few hundred, which a profile likelihood can then exceed.
  pol <- tryCatch(
    maxLik::maxLik(logLik = ll_fun, grad = grad_pen, start = fit$estimate,
                   method = method,
                   control = list(reltol = 1e-12, iterlim = 2000)),
    error = function(e) NULL)
  if (!is.null(pol) && is.finite(pol$maximum) && pol$maximum >= fit$maximum &&
      isTRUE(accept(pol))) {
    for (an in c("bd_starts_ok", "bd_starts_rejected", "bd_loglik_spread"))
      attr(pol, an) <- attr(fit, an)
    fit <- pol
  }

  est_log <- clamp(fit$estimate)
  est_nat <- exp(est_log)
  names(est_nat) <- sub("^log_", "", names(est_log))
  npar <- length(est_nat)
  nobs <- length(time)

  ## ---- boundary screen (thesis PART R) ---------------------------------------
  nm_log <- names(est_log)
  span <- hi_eta[nm_log] - lo_eta[nm_log]
  edge <- (est_log - lo_eta[nm_log]) / span < box_tol |
          (hi_eta[nm_log] - est_log) / span < box_tol
  names(edge) <- names(est_nat)
  at_boundary <- any(edge)

  if (at_boundary) {
    ## A supremum on the edge of the box has no information-based uncertainty.
    vcov_nat <- matrix(NA_real_, npar, npar)
  } else {
    vcov_log <- tryCatch(solve(-fit$hessian),
                         error = function(e) matrix(NA_real_, npar, npar))
    J        <- diag(est_nat, nrow = npar)
    vcov_nat <- J %*% vcov_log %*% t(J)
  }
  rownames(vcov_nat) <- colnames(vcov_nat) <- names(est_nat)

  flag <- character(0)
  if (at_boundary)
    flag <- c(flag, paste0("at box edge (", paste(names(est_nat)[edge], collapse = ","), ")"))
  if (!submodel && (est_nat[["a"]] < 1e-2 || est_nat[["c"]] > 1e2))
    flag <- c(flag, "weakly identified: (a,c) ridge")
  flag <- if (length(flag)) paste(flag, collapse = "; ") else "ok"

  ## The objective may be penalised; the reported log-likelihood never is, so
  ## that AIC, BIC and likelihood ratio tests stay comparable across fits.
  ll_unpen <- loglik_fun(est_log)
  if (!is.finite(ll_unpen)) ll_unpen <- fit$maximum

  out <- list(
    coefficients = est_nat,
    logLik       = ll_unpen,
    penalised_logLik = if (penalty > 0) fit$maximum else NA_real_,
    penalty      = penalty,
    grouped      = isTRUE(grouped),
    delta        = if (isTRUE(grouped)) delta else NA_real_,
    vcov         = vcov_nat,
    npar         = npar,
    nobs         = nobs,
    nevent       = sum(status == 1),
    AIC          = if (at_boundary) NA_real_ else 2 * npar - 2 * ll_unpen,
    BIC          = if (at_boundary) NA_real_ else npar * log(nobs) - 2 * ll_unpen,
    at_boundary  = at_boundary,
    edge_par     = paste(names(est_nat)[edge], collapse = ","),
    flag         = flag,
    box          = box,
    bounds       = bounds,
    submodel_logLik = if (is.null(ed_stage)) NA_real_ else
                        loglik_fun(clamp(ed_stage$estimate)),
    convergence  = fit$code,
    message      = fit$message,
    starts_ok        = .bd_or(attr(fit, "bd_starts_ok"), NA_integer_),
    starts_rejected  = .bd_or(attr(fit, "bd_starts_rejected"), NA_integer_),
    loglik_spread    = .bd_or(attr(fit, "bd_loglik_spread"), NA_real_),
    submodel     = submodel,
    data         = list(time = time, status = status),
    formula      = formula,
    call         = match.call()
  )
  out$diagnostics <- .bd_fit_diagnostics(out)
  class(out) <- "betadanish"

  if (isTRUE(check_identifiability)) .bd_warn_diagnostics(out$diagnostics)

  out
}

#' Log Cell Probability for Grouped Times
#'
#' \eqn{\log\{F(t + \delta/2) - F(t - \delta/2)\}}, falling back to
#' \eqn{\log f(t) + \log\delta} only where the difference of two nearly equal
#' distribution values has cancelled to zero. Capped at `0`, since a
#' probability cannot exceed one.
#'
#' @noRd
.bd_log_cell <- function(t, delta, a, b, c, k) {
  L  <- pmax(t - delta / 2, 0)
  U  <- t + delta / 2
  pc <- pbetadanish(U, a, b, c, k) - pbetadanish(L, a, b, c, k)

  out <- numeric(length(t))
  bad <- !is.finite(pc) | pc <= 0
  if (any(!bad)) out[!bad] <- pmin(log(pc[!bad]), 0)
  if (any(bad))
    out[bad] <- pmin(dbetadanish(t[bad], a, b, c, k, log = TRUE) + log(delta), 0)
  out
}

#' Assemble Convergence and Identifiability Diagnostics
#' @noRd
.bd_fit_diagnostics <- function(fit) {
  se <- .bd_se_from_vcov(fit$vcov, names(fit$coefficients))
  boundary <- isTRUE(fit$at_boundary)

  d <- list(
    converged        = isTRUE(fit$convergence %in% c(0L, 1L, 2L)),
    convergence_code = fit$convergence,
    vcov_singular    = !boundary &&
                       (anyNA(fit$vcov) || any(!is.finite(diag(fit$vcov))) ||
                        any(diag(fit$vcov) <= 0)),
    at_boundary      = boundary,
    edge_par         = .bd_or(fit$edge_par, ""),
    ac_ridge         = isTRUE(grepl("(a,c) ridge", .bd_or(fit$flag, ""), fixed = TRUE)),
    near_b_ridge     = NA,
    b_distance_se    = NA_real_,
    ac_correlation   = NA_real_,
    starts_ok        = .bd_or(fit$starts_ok, NA_integer_),
    starts_rejected  = .bd_or(fit$starts_rejected, NA_integer_),
    loglik_spread    = .bd_or(fit$loglik_spread, NA_real_),
    loglik_per_obs   = if (is.null(fit$nobs) || !isTRUE(fit$nobs > 0)) NA_real_
                       else as.numeric(fit$logLik) / fit$nobs
  )

  if ("b" %in% names(fit$coefficients) && is.finite(se[["b"]]) && se[["b"]] > 0) {
    d$b_distance_se <- abs(fit$coefficients[["b"]] - 1) / se[["b"]]
    d$near_b_ridge  <- d$b_distance_se < 2
  }

  if (!fit$submodel && all(c("a", "c") %in% rownames(fit$vcov))) {
    vaa <- fit$vcov["a", "a"]; vcc <- fit$vcov["c", "c"]
    vac <- fit$vcov["a", "c"]
    if (is.finite(vaa) && is.finite(vcc) && is.finite(vac) && vaa > 0 && vcc > 0)
      d$ac_correlation <- vac / sqrt(vaa * vcc)
  }
  d
}

#' Emit Identifiability Warnings
#'
#' Every field is normalised to a length-one value of the expected type before
#' any condition is evaluated. A fitted object saved before a later diagnostic
#' existed carries no `starts_rejected` or `loglik_spread`, and must degrade
#' quietly rather than error. Guarding each condition separately is not enough:
#' `isTRUE()` absorbs `NA` and zero-length results, but `abs(NULL)` raises a
#' non-numeric argument error before `isTRUE()` ever sees it.
#'
#' @param d A diagnostics list, possibly incomplete.
#' @noRd
.bd_warn_diagnostics <- function(d) {
  if (!is.list(d)) return(invisible(NULL))

  num1 <- function(nm) {
    v <- d[[nm]]
    if (is.null(v) || length(v) != 1L || !is.numeric(v)) NA_real_ else as.numeric(v)
  }
  flag <- function(nm) isTRUE(d[[nm]])

  code       <- if (is.null(d$convergence_code) ||
                    length(d$convergence_code) != 1L) NA else d$convergence_code
  b_dist     <- num1("b_distance_se")
  ac_cor     <- num1("ac_correlation")
  n_rejected <- num1("starts_rejected")
  spread     <- num1("loglik_spread")

  if (!flag("converged"))
    warning("The optimiser reported code ", code,
            "; treat the estimates as provisional and increase n_starts.",
            call. = FALSE)

  if (flag("at_boundary")) {
    ep <- d$edge_par
    ep <- if (is.character(ep) && length(ep) == 1L && nzchar(ep)) ep else "a parameter"
    warning(sprintf(paste0("The fit reached the edge of the parameter box ",
                           "(%s). Its log-likelihood is a supremum, not an ",
                           "attained maximum, so AIC, BIC, standard errors ",
                           "and goodness-of-fit statistics are not reported. ",
                           "For the four-parameter model the ED submodel ",
                           "(submodel = TRUE) is usually the honest report."),
                    ep), call. = FALSE)
  }

  if (flag("vcov_singular"))
    warning("The observed information matrix is singular or not positive ",
            "definite, so standard errors are unreliable. This usually means ",
            "the likelihood is flat in at least one direction.", call. = FALSE)

  if (flag("near_b_ridge"))
    warning(sprintf(paste0("b-hat is only %.2f standard errors from 1, close to ",
                           "the b = 1 non-identifiability ridge. Consider the ",
                           "ED submodel (submodel = TRUE). See the ",
                           "Identifiability section of ?fit_betadanish."),
                    b_dist), call. = FALSE)

  if (isTRUE(abs(ac_cor) > 0.95))
    warning(sprintf(paste0("The fitted correlation between a-hat and c-hat is ",
                           "%.3f, so effectively only the product c*a is ",
                           "identified. Individual estimates of a and c should ",
                           "not be interpreted."),
                    ac_cor), call. = FALSE)

  if (isTRUE(n_rejected > 0))
    warning(sprintf(paste0("%.0f starting point(s) reached a degenerate ridge ",
                           "and were discarded. The reported fit is the best ",
                           "admissible optimum. If this is most of the grid, ",
                           "the four-parameter model is a poor choice for ",
                           "these data."), n_rejected), call. = FALSE)

  if (isTRUE(spread > 2))
    warning(sprintf(paste0("Accepted optima span %.2f log-likelihood units, ",
                           "so the surface has several local maxima and the ",
                           "reported fit may not be global. Increase ",
                           "n_starts."), spread), call. = FALSE)

  invisible(NULL)
}
