#' The Kumaraswamy distribution in the median parametrization
#'
#' Density, distribution function, quantile function and random generation
#' for the Kumaraswamy distribution parametrized by its median and a shape
#' parameter (Mitnik and Baek, 2013; Bayer et al., 2017), on a bounded support
#' `(lower, upper)`.
#'
#' For the standard support `(0, 1)` the density is
#' \deqn{f(y \mid \mu, \varphi) = \varphi \delta y^{\varphi - 1}
#'   (1 - y^{\varphi})^{\delta - 1}, \qquad
#'   \delta = \frac{\log(0.5)}{\log(1 - \mu^{\varphi})},}
#' with median \eqn{\mu \in (0, 1)} and shape \eqn{\varphi > 0}. The
#' distribution function is \eqn{F(y) = 1 - (1 - y^{\varphi})^{\delta}} and
#' the quantile function is
#' \eqn{F^{-1}(p) = \{1 - (1 - p)^{1/\delta}\}^{1/\varphi}}. By construction
#' \eqn{F(\mu) = 0.5}.
#'
#' For a general support, `x`, `q` **and** `median` are all expressed on the
#' original scale `(lower, upper)`: the variable is mapped to the unit
#' interval by \eqn{(y - a)/(b - a)}, and the density is divided by
#' \eqn{b - a}. The shape parameter is scale free.
#'
#' The support is the open interval: the density is zero at and outside the
#' bounds, and the distribution function is 0 below `lower` and 1 above
#' `upper`. Computations use `log1p()`/`expm1()` so that tail probabilities
#' and quantiles remain accurate near the bounds.
#'
#' Arguments are recycled to the length of the longest one. Invalid parameters
#' (`median` outside `(lower, upper)`, non-positive `shape`, or
#' `lower >= upper`) give `NaN` with a warning; missing values propagate.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations. If `length(n) > 1`, the length is taken
#'   to be the number required.
#' @param median Median of the distribution, on the same scale as `x`.
#' @param shape Shape parameter (\eqn{\varphi > 0}).
#' @param lower,upper Bounds of the support.
#' @param log,log.p Logical; if `TRUE`, probabilities/densities are given on
#'   the log scale.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#'
#' @return `dkuma()` gives the density, `pkuma()` the distribution function,
#'   `qkuma()` the quantile function and `rkuma()` random deviates. The length
#'   of the result is the length of the longest argument (for `rkuma()`, `n`).
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. *Journal of Hydrology*, 46, 79--88.
#'
#' Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
#' median-dispersion re-parameterizations for regression modeling and
#' simulation-based estimation. *Statistical Papers*, 54, 177--192.
#'
#' Bayer, F. M., Bayer, D. M. and Pumi, G. (2017). Kumaraswamy autoregressive
#' moving average models for double bounded environmental data. *Journal of
#' Hydrology*, 555, 385--396.
#'
#' @examples
#' dkuma(c(0.2, 0.5, 0.8), median = 0.5, shape = 3)
#' pkuma(0.5, median = 0.5, shape = 3) # 0.5 by construction
#' qkuma(c(0.1, 0.5, 0.9), median = 0.7, shape = 10)
#'
#' # NDVI-like support (-1, 1): median and quantiles on the original scale
#' qkuma(0.5, median = 0.4, shape = 5, lower = -1, upper = 1)
#'
#' set.seed(1)
#' rkuma(5, median = 0.3, shape = 2)
#' @name kumaraswamy
NULL

kuma_delta <- function(m, shape) {
  log(0.5) / log1p(-m^shape)
}

kuma_invalid <- function(m, shape, lower, upper) {
  !is.na(m) & !is.na(shape) & !is.na(lower) & !is.na(upper) &
    (lower >= upper | m <= 0 | m >= 1 | shape <= 0 | !is.finite(shape))
}

#' @rdname kumaraswamy
#' @export
dkuma <- function(x, median, shape, lower = 0, upper = 1, log = FALSE) {
  a <- recycle_args(x, median, shape, lower, upper)
  x <- a[[1]]
  lower <- a[[4]]
  width <- a[[5]] - lower
  z <- (x - lower) / width
  m <- (a[[2]] - lower) / width
  shape <- a[[3]]

  out <- rep(NA_real_, length(z))
  invalid <- kuma_invalid(m, shape, lower, a[[5]])
  ok <- !is.na(z) & !is.na(m) & !is.na(shape) & !invalid
  inside <- ok & z > 0 & z < 1
  delta <- kuma_delta(m[inside], shape[inside])
  out[inside] <- log(shape[inside]) + log(delta) +
    (shape[inside] - 1) * log(z[inside]) +
    (delta - 1) * kuma_log1m_zphi(z[inside], shape[inside]) -
    log(width[inside])
  out[ok & !inside] <- -Inf
  out[invalid] <- NaN
  if (any(invalid)) warning("NaNs produced", call. = FALSE)
  if (log) out else exp(out)
}

# log(1 - z^shape), accurate for z close to 0 and to 1.
kuma_log1m_zphi <- function(z, shape) {
  log1mexp(shape * log(z))
}

# Log of the upper tail probability, log P(X > x), on the unit scale.
kuma_log_survival <- function(z, m, shape) {
  kuma_delta(m, shape) * kuma_log1m_zphi(z, shape)
}

#' @rdname kumaraswamy
#' @export
pkuma <- function(q, median, shape, lower = 0, upper = 1,
                  lower.tail = TRUE, log.p = FALSE) {
  a <- recycle_args(q, median, shape, lower, upper)
  lower <- a[[4]]
  width <- a[[5]] - lower
  z <- (a[[1]] - lower) / width
  m <- (a[[2]] - lower) / width
  shape <- a[[3]]

  log_surv <- rep(NA_real_, length(z))
  invalid <- kuma_invalid(m, shape, lower, a[[5]])
  ok <- !is.na(z) & !is.na(m) & !is.na(shape) & !invalid
  inside <- ok & z > 0 & z < 1
  log_surv[inside] <- kuma_log_survival(z[inside], m[inside], shape[inside])
  log_surv[ok & z <= 0] <- 0
  log_surv[ok & z >= 1] <- -Inf

  out <- if (lower.tail) log1mexp(log_surv) else log_surv
  out[invalid] <- NaN
  if (any(invalid)) warning("NaNs produced", call. = FALSE)
  if (log.p) out else exp(out)
}

#' @rdname kumaraswamy
#' @export
qkuma <- function(p, median, shape, lower = 0, upper = 1,
                  lower.tail = TRUE, log.p = FALSE) {
  a <- recycle_args(p, median, shape, lower, upper)
  p <- a[[1]]
  lower <- a[[4]]
  width <- a[[5]] - lower
  m <- (a[[2]] - lower) / width
  shape <- a[[3]]

  bad_p <- !is.na(p) & (if (log.p) p > 0 else (p < 0 | p > 1))
  p_ok <- ifelse(bad_p, NA_real_, p)
  # log of the upper tail probability 1 - F(x)
  log_surv <- if (log.p) {
    if (lower.tail) log1mexp(p_ok) else p_ok
  } else {
    if (lower.tail) log1p(-p_ok) else log(p_ok)
  }
  invalid <- kuma_invalid(m, shape, lower, a[[5]]) | bad_p
  ok <- !is.na(log_surv) & !is.na(m) & !is.na(shape) & !invalid

  out <- rep(NA_real_, length(p))
  delta <- kuma_delta(m[ok], shape[ok])
  z <- (-expm1(log_surv[ok] / delta))^(1 / shape[ok])
  out[ok] <- lower[ok] + width[ok] * z
  out[invalid] <- NaN
  if (any(invalid)) warning("NaNs produced", call. = FALSE)
  out
}

#' @rdname kumaraswamy
#' @export
rkuma <- function(n, median, shape, lower = 0, upper = 1) {
  if (length(n) > 1L) n <- length(n)
  if (!is.numeric(n) || is.na(n) || n < 0) {
    stop("`n` must be a non-negative number.", call. = FALSE)
  }
  n <- floor(n)
  if (n == 0) return(numeric(0))
  if (length(median) == 0L || length(shape) == 0L) {
    warning("NAs produced", call. = FALSE)
    return(rep(NA_real_, n))
  }
  u <- stats::runif(n)
  qkuma(u, rep_len(median, n), rep_len(shape, n), lower, upper)
}
