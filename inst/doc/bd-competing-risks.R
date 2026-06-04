## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  warning = FALSE,
  message = FALSE
)
have_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
knitr::opts_chunk$set(eval = have_cmprsk)

## ----data---------------------------------------------------------------------
library(BetaDanish)
set.seed(2026)
n  <- 400
T1 <- rbetadanish(n, a = 1.2, b = 1.5, c = 1.0, k = 0.4)
T2 <- rbetadanish(n, a = 1.0, b = 2.0, c = 1.0, k = 0.2)
C  <- stats::rexp(n, 0.05)
time  <- pmin(T1, T2, C)
cause <- ifelse(time == C, 0L, ifelse(T1 <= T2, 1L, 2L))
table(cause)

## ----fit----------------------------------------------------------------------
fit <- fit_bd_competing(time = time, cause = cause)
print(fit)

## ----cif----------------------------------------------------------------------
res <- cif_compare(fit, plot = TRUE)

## ----gray---------------------------------------------------------------------
if (!is.null(res$gray_test)) {
  print(res$gray_test)
} else {
  cat("Gray's test was not produced.\n")
}

