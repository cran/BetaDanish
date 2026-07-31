#' Simulate Data from the Beta-Danish Distribution
#'
#' Generates synthetic survival data from the Beta-Danish distribution, with
#' optional right-censoring.
#'
#' @param n Integer; number of observations to simulate.
#' @param a,b,c,k Numeric; parameters of the Beta-Danish distribution.
#' @param censor_rate Numeric; rate parameter for the exponential censoring
#'   distribution. If `0` (default), no censoring is applied.
#' @param seed Integer; optional seed for reproducibility.
#'
#' @return A data frame with columns `time` and `status`.
#' @export
#'
#' @examples
#' # Simulate complete data
#' dat <- simulate_bd_data(n = 100, a = 1.5, b = 2, c = 3, k = 0.5)
#'
#' # Simulate censored data
#' dat_cens <- simulate_bd_data(n = 100, a = 1.5, b = 2, c = 3, k = 0.5, censor_rate = 0.1)
simulate_bd_data <- function(n, a, b, c, k, censor_rate = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate true event times
  true_times <- rbetadanish(n, a, b, c, k)

  if (censor_rate > 0) {
    # Generate censoring times from Exponential distribution
    censor_times <- stats::rexp(n, rate = censor_rate)

    # Observed time is the minimum of true time and censoring time
    time <- pmin(true_times, censor_times)
    status <- ifelse(true_times <= censor_times, 1, 0)
  } else {
    time <- true_times
    status <- rep(1, n)
  }

  data.frame(time = time, status = status)
}

#' Simulate Beta-Danish Cure Data
#'
#' Generates synthetic survival data from a Beta-Danish mixture or promotion-time
#' cure model, incorporating covariates.
#'
#' @param n Integer; number of observations.
#' @param type Character; `"mixture"` or `"promotion"`.
#' @param a,b,c Numeric; baseline shape parameters.
#' @param delta Numeric vector; coefficients for the latency scale `k`.
#' @param gamma Numeric vector; coefficients for the incidence/cure component.
#' @param X Matrix; design matrix for latency (must match length of `delta`).
#' @param Z Matrix; design matrix for incidence (must match length of `gamma`).
#' @param target_censor Numeric; target proportion of censoring to calibrate the
#'   exponential censoring rate. Default is 0.3.
#' @param seed Integer; optional seed.
#'
#' @return A list containing the simulated `data` (time, status), the `cured`
#'   indicator, and the true parameters.
#' @export
simulate_bd_cure_data <- function(n, type = c("mixture", "promotion"), a = 1, b = 2, c = 1.5,
                                  delta, gamma, X, Z, target_censor = 0.3, seed = NULL) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  # Latency scale
  eta_k <- as.numeric(X %*% delta)
  k_i <- exp(eta_k)

  T_true <- rep(Inf, n)
  cured <- rep(FALSE, n)

  if (type == "mixture") {
    # Susceptible probability
    eta_pi <- as.numeric(Z %*% gamma)
    pi_i <- 1 / (1 + exp(-eta_pi))

    cured <- stats::rbinom(n, 1, prob = 1 - pi_i) == 1

    # Generate times for susceptible
    n_sus <- sum(!cured)
    if (n_sus > 0) {
      T_true[!cured] <- rbetadanish(n_sus, a, b, c, k_i[!cured])
    }

  } else {
    # Promotion-time
    eta_th <- as.numeric(Z %*% gamma)
    theta_i <- exp(eta_th)

    # Number of clonogenic cells
    N_cells <- stats::rpois(n, lambda = theta_i)
    cured <- (N_cells == 0)

    for (i in which(!cured)) {
      # Time is the minimum of N_cells latent times
      latent_times <- rbetadanish(N_cells[i], a, b, c, k_i[i])
      T_true[i] <- min(latent_times)
    }
  }

  # Calibrate censoring to hit target_censor
  finite_T <- T_true[is.finite(T_true)]
  if (length(finite_T) > 10) {
    obj_fun <- function(rate) {
      C <- stats::rexp(length(finite_T), rate)
      mean(C < finite_T) - target_censor
    }
    opt_rate <- tryCatch(stats::uniroot(obj_fun, c(1e-6, 10))$root, error = function(e) 0.05)
  } else {
    opt_rate <- 0.05
  }

  C_times <- stats::rexp(n, rate = opt_rate)
  time <- pmin(T_true, C_times)
  status <- ifelse(T_true <= C_times, 1, 0)

  list(
    data = data.frame(time = time, status = status),
    cured = cured,
    true_params = list(a=a, b=b, c=c, delta=delta, gamma=gamma)
  )
}
