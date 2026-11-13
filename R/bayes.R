#' Bayesian Estimation for the Beta-Danish Distribution
#'
#' Samples from the posterior of the Beta-Danish or Exponentiated Danish
#' parameters with vague \eqn{\Gamma(0.01, 0.01)} priors on the positive
#' parameters. The default sampler is the multi-chain adaptive Metropolis
#' algorithm of the thesis, written in base R, with a convergence check on
#' every parameter.
#'
#' @param time Numeric vector of observed times.
#' @param status Numeric vector of event indicators (1 = event, 0 = right-censored).
#' @param submodel Logical; \code{TRUE} for the 3-parameter ED submodel,
#'   \code{FALSE} for the 4-parameter full model.
#' @param burnin Burn-in iterations (per chain for the adaptive sampler).
#' @param mcmc Post-burnin iterations (per chain for the adaptive sampler).
#' @param tune Random-walk tuning parameter, used only by
#'   \code{method = "metrop1R"}.
#' @param theta_init Optional starting values on the log scale. For the
#'   adaptive sampler the first chain starts here and the others at random
#'   perturbations of it; by default the maximum-likelihood estimate from
#'   [fit_betadanish()] is used.
#' @param seed Optional integer seed.
#' @param verbose Integer; passed to MCMCmetrop1R (0 = silent). Ignored by the
#'   adaptive sampler.
#' @param method \code{"adaptive"} (default) for the multi-chain adaptive
#'   Metropolis sampler, or \code{"metrop1R"} for the single-chain sampler of
#'   \pkg{MCMCpack} used by earlier versions.
#' @param n_chains Number of chains for the adaptive sampler. Default 4.
#' @param thin Thinning interval for the adaptive sampler. Default 1.
#' @param target_acc Target acceptance rate of the scale adaptation. Default
#'   0.234.
#' @param rhat_max,ess_min Convergence gate: every parameter must have a
#'   potential scale reduction factor at most \code{rhat_max} (default 1.05)
#'   and an effective sample size at least \code{ess_min} (default 400).
#'
#' @return An object of class \code{"bd_bayes"}. For the adaptive sampler its
#'   components are \code{draws} (the retained natural-scale draws of all
#'   chains stacked, a \pkg{coda} \code{mcmc} object when \pkg{coda} is
#'   installed), \code{chains} (a list with one natural-scale matrix per
#'   chain), \code{summary} (a data frame of posterior means, standard
#'   deviations, quantiles and 95\% HPD limits; for the full model it also
#'   has a row for the identified composite \code{a_times_c}), \code{HPD} (a
#'   matrix of 95\% HPD intervals), \code{diagnostics} (R-hat, effective
#'   sample size and acceptance rate per parameter), \code{converged},
#'   \code{acceptance}, \code{settings}, \code{method}, \code{submodel} and
#'   \code{call}. For \code{method = "metrop1R"} the object is as in earlier
#'   versions: \code{draws} (mcmc object), \code{summary}, \code{HPD},
#'   \code{submodel}, \code{call} and \code{method}.
#'
#' @details
#' **Adaptive sampler.** Each chain is a random-walk Metropolis sampler on the
#' log-parameter scale. During the first 85\% of the burn-in the global scale
#' of the proposal is adapted by a Robbins-Monro step towards the target
#' acceptance rate, every 200 iterations with gain
#' \eqn{\min\{0.6, (200/i)^{0.6}\}}; once enough iterations have accumulated
#' the proposal covariance is replaced by the learned covariance of the chain,
#' scaled by \eqn{2.38^2/d}, and the global scale is re-based to 1 once and
#' keeps adapting on top of it. Adaptation stops inside the burn-in, so every
#' retained draw comes from a fixed kernel. The parameters are confined to
#' \eqn{[10^{-7}, 10^{7}]} on the natural scale, which keeps a chain from
#' drifting off along the \eqn{(a, c)} ridge under the vague priors.
#'
#' Convergence is judged on the log scale by the potential scale reduction
#' factor and the effective sample size, computed with \pkg{coda} when it is
#' installed and by base-R equivalents otherwise. A fit whose chains fail the
#' gate is returned with \code{converged = FALSE} and a warning; increase
#' \code{burnin} and \code{mcmc} before using it. For the four-parameter model
#' the posterior of \eqn{ac} is usually far better determined than those of
#' \eqn{a} and \eqn{c}, which is why it is summarised separately.
#'
#' **MCMCpack sampler.** \code{method = "metrop1R"} calls
#' \code{MCMCpack::MCMCmetrop1R()} with the fixed tuning constant \code{tune};
#' it requires \pkg{MCMCpack} and \pkg{coda} (Suggests). The thesis found its
#' acceptance rate for the four-parameter model to be as low as 4--9\%, which
#' is why it is no longer the default.
#'
#' @references
#' Ahmad, B. (2026). *Modeling Diverse Survival Patterns: The Development and
#' Characterization of a New Four-Parameter Lifetime Distribution* (Ph.D.
#' thesis). Allama Iqbal Open University, Islamabad, Pakistan.
#'
#' @seealso [plot.bd_bayes()]
#'
#' @examples
#' \donttest{
#' dat <- simulate_bd_data(100, a = 1, b = 3, c = 2, k = 0.5, seed = 1)
#' fit <- bayes_betadanish(dat$time, dat$status, submodel = TRUE,
#'                         burnin = 1000, mcmc = 2000, seed = 1)
#' fit$summary
#' fit$diagnostics
#' }
#'
#' @export
bayes_betadanish <- function(time, status = NULL, submodel = TRUE,
                             burnin = 5000, mcmc = 15000, tune = 0.5,
                             theta_init = NULL, seed = NULL, verbose = 0,
                             method = c("adaptive", "metrop1R"),
                             n_chains = 4L, thin = 1L, target_acc = 0.234,
                             rhat_max = 1.05, ess_min = 400) {
  method <- match.arg(method)
  time <- as.numeric(time)
  if (is.null(status)) status <- rep(1, length(time))
  status <- as.integer(status)
  if (length(status) != length(time))
    stop("time and status must have the same length.")

  nm  <- if (submodel) c("b", "c", "k") else c("a", "b", "c", "k")
  dgamma_log <- function(x, shape, rate)
    stats::dgamma(x, shape = shape, rate = rate, log = TRUE)

  ## Log posterior on the log scale: likelihood, Gamma(0.01, 0.01) priors and
  ## the Jacobian of the log transform.
  loglik <- function(th) {
    a <- if (submodel) 1 else th[["a"]]
    log_f <- dbetadanish(time, a = a, b = th[["b"]], c = th[["c"]], k = th[["k"]],
                         log = TRUE)
    log_S <- pbetadanish(time, a = a, b = th[["b"]], c = th[["c"]], k = th[["k"]],
                         lower.tail = FALSE, log.p = TRUE)
    sum(log_f[status == 1]) + sum(log_S[status == 0])
  }

  if (method == "metrop1R") {
    if (!requireNamespace("MCMCpack", quietly = TRUE))
      stop("method = \"metrop1R\" requires MCMCpack. Install with install.packages('MCMCpack').")
    if (!requireNamespace("coda", quietly = TRUE))
      stop("method = \"metrop1R\" requires coda. Install with install.packages('coda').")
    logpost_m1r <- function(theta) {
      th <- stats::setNames(exp(theta), nm)
      if (any(!is.finite(th)) || any(th <= 0)) return(-Inf)
      ll <- suppressWarnings(loglik(th))
      lp <- sum(dgamma_log(th, 0.01, 0.01))
      ll + lp + sum(theta)
    }
    if (is.null(theta_init))
      theta_init <- if (submodel) log(c(2, 2, 1)) else log(c(1.5, 2, 2, 1))
    names(theta_init) <- paste0("log_", nm)
    if (!is.null(seed)) set.seed(seed)
    post <- MCMCpack::MCMCmetrop1R(logpost_m1r, theta.init = theta_init,
                                   burnin = burnin, mcmc = mcmc,
                                   tune = tune, verbose = verbose)
    draws_nat <- exp(as.matrix(post))
    colnames(draws_nat) <- nm
    draws_mcmc <- coda::as.mcmc(draws_nat)
    out <- list(
      draws    = draws_mcmc,
      summary  = summary(draws_mcmc),
      HPD      = coda::HPDinterval(draws_mcmc),
      submodel = submodel,
      method   = method,
      call     = match.call()
    )
    class(out) <- "bd_bayes"
    return(out)
  }

  ## ---- adaptive multi-chain Metropolis (thesis PART R) ------------------------
  n_chains <- as.integer(n_chains); thin <- as.integer(thin)
  burnin <- as.integer(burnin); mcmc <- as.integer(mcmc)
  if (n_chains < 1L || thin < 1L || burnin < 10L || mcmc < 10L)
    stop("'n_chains' and 'thin' must be at least 1, and 'burnin' and 'mcmc' ",
         "at least 10.", call. = FALSE)

  lo <- log(1e-7); hi <- log(1e7)
  logpost <- function(eta) {
    if (any(!is.finite(eta)) || any(eta < lo) || any(eta > hi)) return(-Inf)
    th <- stats::setNames(exp(eta), nm)
    pr <- sum(dgamma_log(th, 0.01, 0.01))
    if (!is.finite(pr)) return(-Inf)
    v <- pr + suppressWarnings(loglik(th)) + sum(eta)
    if (is.finite(v)) v else -Inf
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(theta_init)) {
    mle <- tryCatch(suppressWarnings(
      fit_betadanish(survival::Surv(time, status) ~ 1,
                     data = data.frame(time = time, status = status),
                     submodel = submodel, n_starts = 0,
                     check_identifiability = FALSE)),
      error = function(e) NULL)
    theta_init <- if (!is.null(mle)) log(mle$coefficients[nm]) else
      if (submodel) log(c(2, 2, 1 / stats::median(time))) else
        log(c(1.5, 2, 2, 1 / stats::median(time)))
  }
  if (length(theta_init) != length(nm))
    stop("'theta_init' must have length ", length(nm), " (log scale: ",
         paste(nm, collapse = ", "), ").", call. = FALSE)
  theta_init <- stats::setNames(as.numeric(theta_init), nm)
  if (!is.finite(logpost(theta_init)))
    stop("The log posterior is not finite at the starting values; supply ",
         "'theta_init' on the log scale.", call. = FALSE)

  cfg <- list(n_iter = burnin + mcmc, burn = burnin, thin = thin,
              target_acc = target_acc, adapt_upto = 0.85)
  adapt_end <- as.integer(floor(cfg$burn * cfg$adapt_upto))
  cfg$window    <- as.integer(min(200L, max(10L, adapt_end %/% 10L)))
  cfg$cov_start <- as.integer(min(3000L, max(cfg$window, adapt_end %/% 2L)))

  chains <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    st <- if (i == 1L) theta_init else theta_init + stats::rnorm(length(nm), 0, 1)
    names(st) <- nm
    if (!is.finite(logpost(st))) st <- theta_init
    chains[[i]] <- .bd_adaptive_chain(logpost, st, nm, cfg)
  }

  diag_tab <- .bd_chain_diagnostics(chains, nm)
  converged <- all(is.finite(diag_tab$ESS)) && all(diag_tab$ESS >= ess_min) &&
    (n_chains < 2L || (all(is.finite(diag_tab$Rhat)) && all(diag_tab$Rhat <= rhat_max)))
  if (!converged)
    warning(sprintf(paste0("The chains do not pass the convergence gate ",
                           "(R-hat <= %.2f and effective sample size >= %g ",
                           "for every parameter). Increase 'burnin' and ",
                           "'mcmc' before using the posterior summaries."),
                    rhat_max, ess_min), call. = FALSE)

  chains_nat <- lapply(chains, function(M) {
    E <- exp(M); colnames(E) <- nm; attr(E, "acceptance_rate") <- NULL; E
  })
  draws_nat <- do.call(rbind, chains_nat)
  summ <- .bd_posterior_summary(draws_nat, submodel)
  hpd  <- as.matrix(summ[summ$parameter %in% nm, c("HPD_lower", "HPD_upper")])
  dimnames(hpd) <- list(nm, c("lower", "upper"))

  out <- list(
    draws       = if (requireNamespace("coda", quietly = TRUE))
                    coda::as.mcmc(draws_nat) else draws_nat,
    chains      = chains_nat,
    summary     = summ,
    HPD         = hpd,
    diagnostics = diag_tab,
    converged   = converged,
    acceptance  = vapply(chains, function(M) attr(M, "acceptance_rate"), numeric(1)),
    settings    = list(n_chains = n_chains, burnin = burnin, mcmc = mcmc,
                       thin = thin, target_acc = target_acc,
                       adapt_end = adapt_end, window = cfg$window,
                       cov_start = cfg$cov_start, rhat_max = rhat_max,
                       ess_min = ess_min),
    submodel    = submodel,
    method      = method,
    call        = match.call()
  )
  class(out) <- "bd_bayes"
  out
}

#' One Chain of the Adaptive Metropolis Sampler
#'
#' Ported from `.r_chain()` of the thesis master script. Random-walk
#' Metropolis on the log scale; Robbins-Monro adaptation of the global scale
#' towards `target_acc` every `window` iterations; the learned covariance
#' replaces the spherical proposal from `cov_start` on, with the global scale
#' re-based to 1 once (not at every window) and then adapted further.
#' Adaptation stops at `adapt_upto * burn`, inside the burn-in.
#'
#' @return Matrix of retained log-scale draws with attribute
#'   `acceptance_rate`.
#' @noRd
.bd_adaptive_chain <- function(logpost, start, nm, cfg) {
  d <- length(start); cur <- start; names(cur) <- nm
  cur_lp <- logpost(cur)
  if (!is.finite(cur_lp)) stop("log posterior not finite at the start")
  lambda <- 2.4 / sqrt(d)
  R <- diag(0.05, d)
  mu <- cur; Sr <- diag(1e-4, d)
  acc <- 0L; win_acc <- 0L; win <- cfg$window; cov_on <- FALSE
  adapt_end <- as.integer(floor(cfg$burn * cfg$adapt_upto))
  keep_at <- seq.int(cfg$burn + 1L, cfg$n_iter, by = cfg$thin)
  keep <- matrix(NA_real_, length(keep_at), d); colnames(keep) <- nm; j <- 0L
  for (i in seq_len(cfg$n_iter)) {
    prop <- cur + lambda * as.vector(crossprod(R, stats::rnorm(d)))
    lp   <- logpost(prop)
    if (is.finite(lp) && log(stats::runif(1)) < lp - cur_lp) {
      cur <- prop; cur_lp <- lp; acc <- acc + 1L; win_acc <- win_acc + 1L
    }
    if (i <= adapt_end) {
      dl <- cur - mu; mu <- mu + dl / i
      Sr <- Sr + (tcrossprod(dl, cur - mu) - Sr) / i
      if (i %% win == 0L) {
        lambda <- lambda * exp(min(0.6, (win / i)^0.6) * (win_acc / win - cfg$target_acc))
        lambda <- min(max(lambda, 1e-6), 50); win_acc <- 0L
        if (i >= cfg$cov_start) {
          Sc <- (2.38^2 / d) * (Sr + diag(1e-8, d))
          Rc <- tryCatch(chol(Sc), error = function(e) NULL)
          if (!is.null(Rc) && all(is.finite(Rc))) {
            R <- Rc
            if (!cov_on) { lambda <- 1; cov_on <- TRUE }
          }
        }
      }
    }
    if (i > cfg$burn && ((i - cfg$burn - 1L) %% cfg$thin == 0L)) {
      j <- j + 1L; if (j <= nrow(keep)) keep[j, ] <- cur
    }
  }
  keep <- keep[seq_len(min(j, nrow(keep))), , drop = FALSE]
  attr(keep, "acceptance_rate") <- acc / cfg$n_iter
  keep
}

#' Convergence Diagnostics for a List of Chains
#'
#' Potential scale reduction factor and effective sample size per parameter on
#' the log scale: \pkg{coda} when installed, base-R equivalents otherwise
#' (Gelman-Rubin and batch means). Ported from `.r_diag()` of the thesis
#' master script.
#'
#' @noRd
.bd_chain_diagnostics <- function(chains, nm) {
  np <- length(nm)
  rh <- rep(NA_real_, np); es <- rep(NA_real_, np)
  if (requireNamespace("coda", quietly = TRUE)) {
    ml <- tryCatch(coda::as.mcmc.list(lapply(chains, function(M) {
      M2 <- matrix(as.numeric(M), ncol = np, dimnames = list(NULL, nm))
      coda::as.mcmc(M2)
    })), error = function(e) NULL)
    if (!is.null(ml)) {
      if (length(chains) >= 2L) {
        v <- tryCatch(coda::gelman.diag(ml, autoburnin = FALSE,
                                        multivariate = FALSE)$psrf[, 1],
                      error = function(e) NULL)
        if (!is.null(v) && length(v) == np) rh <- as.numeric(v)
      }
      v <- tryCatch(as.numeric(coda::effectiveSize(ml)), error = function(e) NULL)
      if (!is.null(v) && length(v) == np) es <- v
    }
  }
  if (all(is.na(rh)) && length(chains) >= 2L)
    rh <- vapply(seq_len(np), function(j) .bd_rhat(chains, j), numeric(1))
  if (all(is.na(es)))
    es <- vapply(seq_len(np), function(j)
      sum(vapply(chains, function(M) .bd_ess(M[, j]), numeric(1))), numeric(1))
  data.frame(parameter = nm, Rhat = rh, ESS = es,
             acceptance_rate = mean(vapply(chains, function(M)
               attr(M, "acceptance_rate"), numeric(1))),
             stringsAsFactors = FALSE)
}

#' @noRd
.bd_rhat <- function(chains, j) {
  xs <- lapply(chains, function(M) M[, j]); m <- length(xs); n <- min(lengths(xs))
  if (m < 2L || n < 10L) return(NA_real_)
  xs <- lapply(xs, function(v) v[seq_len(n)])
  W <- mean(vapply(xs, stats::var, numeric(1)))
  B <- n * stats::var(vapply(xs, mean, numeric(1)))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt((((n - 1) / n) * W + B / n) / W)
}

#' @noRd
.bd_ess <- function(v) {
  n <- length(v); if (n < 20L) return(NA_real_)
  b <- max(2L, floor(sqrt(n))); nb <- floor(n / b); if (nb < 2L) return(NA_real_)
  bm <- vapply(seq_len(nb), function(i) mean(v[((i - 1) * b + 1):(i * b)]), numeric(1))
  s2 <- stats::var(v); sb <- stats::var(bm)
  if (!is.finite(s2) || s2 <= 0 || !is.finite(sb) || sb <= 0) return(NA_real_)
  as.numeric(n * s2 / (b * sb))
}

#' Highest Posterior Density Interval from Draws
#' @noRd
.bd_hpd <- function(v, lev = 0.95) {
  v <- sort(v[is.finite(v)]); n <- length(v)
  if (n < 10L) return(c(NA_real_, NA_real_))
  m <- max(1L, floor(lev * n)); if (m >= n) return(c(v[1], v[n]))
  w <- v[(m + 1L):n] - v[1L:(n - m)]; i <- which.min(w); c(v[i], v[i + m])
}

#' Posterior Summary Table
#' @noRd
.bd_posterior_summary <- function(dr, submodel) {
  one <- function(v, lab) {
    v <- v[is.finite(v)]
    q <- stats::quantile(v, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
    h <- .bd_hpd(v)
    data.frame(parameter = lab, mean = mean(v), sd = stats::sd(v),
               q2.5 = q[1], q25 = q[2], median = q[3], q75 = q[4],
               q97.5 = q[5], HPD_lower = h[1], HPD_upper = h[2],
               stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(colnames(dr), function(p) one(dr[, p], p)))
  if (!submodel) out <- rbind(out, one(dr[, "a"] * dr[, "c"], "a_times_c"))
  rownames(out) <- NULL
  out
}

#' @export
print.bd_bayes <- function(x, ...) {
  cat("\nBeta-Danish Bayesian Fit (",
      if (x$submodel) "3-Parameter ED Submodel"
      else "Full 4-Parameter Beta-Danish", ")\n", sep = "")
  if (identical(x$method, "adaptive")) {
    s <- x$settings
    cat("Adaptive Metropolis: ", s$n_chains, " chain(s), ", s$burnin,
        " burn-in + ", s$mcmc, " iterations each, thin ", s$thin, "\n\n", sep = "")
    tab <- x$summary
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = 5)
    print(tab, row.names = FALSE)
    cat("\nConvergence diagnostics (log scale):\n")
    dg <- x$diagnostics
    dg$Rhat <- round(dg$Rhat, 4); dg$ESS <- round(dg$ESS, 1)
    dg$acceptance_rate <- round(dg$acceptance_rate, 3)
    print(dg, row.names = FALSE)
    cat("\nConverged (R-hat <= ", s$rhat_max, ", ESS >= ", s$ess_min, "): ",
        if (isTRUE(x$converged)) "yes" else "NO", "\n", sep = "")
  } else {
    cat("Random-walk Metropolis via MCMCpack::MCMCmetrop1R\n\n")
    print(x$summary)
    cat("\n95% HPD intervals:\n")
    print(x$HPD)
  }
  invisible(x)
}
