#' Diagnostics of a 2-D KARMA fit
#'
#' Summarizes the quantile residuals over the effective region, which should
#' be approximately independent standard normal when the model is adequate
#' (Dunn and Smyth, 1996): their mean, standard deviation, skewness and
#' kurtosis (3 under normality), the proportion outside \eqn{\pm 3}, and their
#' sample spatial autocorrelations at the first-order causal lags. Also
#' reports convergence and the reliability of the information matrix.
#'
#' The autocorrelation at lag \eqn{(i,j)} is the Pearson correlation between
#' \eqn{r[n,m]} and \eqn{r[n-i,m-j]} over all pairs available in the
#' effective region. Values close to zero suggest that the spatial dependence
#' has been captured; no formal test is implied.
#'
#' @param object A `karma2d_fit` object.
#' @param lags Data frame with columns `row_lag` and `col_lag`.
#' @return An object of class `karma2d_diagnostics`.
#' @examples
#' y <- karma2d_simulate(30, 30, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 21)
#' karma2d_diagnostics(karma2d_fit(y, order = c(1, 0)))
#' @export
karma2d_diagnostics <- function(object, lags = karma2d_lags(1)) {
  stopifnot(inherits(object, "karma2d_fit"))
  r <- object$residuals
  rows <- object$effective_region$rows
  cols <- object$effective_region$cols
  r_eff <- r[rows, cols, drop = FALSE]
  v <- as.vector(r_eff)
  v <- v[is.finite(v)]
  m <- mean(v)
  s <- stats::sd(v)
  central <- v - m
  moments <- c(
    mean = m, sd = s,
    skewness = mean(central^3) / mean(central^2)^1.5,
    kurtosis = mean(central^4) / mean(central^2)^2,
    prop_outside_3 = mean(abs(v) > 3)
  )
  acf <- vapply(seq_len(nrow(lags)), function(s) {
    i <- lags$row_lag[s]
    j <- lags$col_lag[s]
    nr <- nrow(r_eff)
    nc <- ncol(r_eff)
    if (i >= nr || j >= nc) return(NA_real_)
    a <- r_eff[(i + 1):nr, (j + 1):nc]
    b <- r_eff[1:(nr - i), 1:(nc - j)]
    stats::cor(as.vector(a), as.vector(b), use = "complete.obs")
  }, numeric(1))
  structure(
    list(residual_moments = moments,
         spatial_acf = data.frame(row_lag = lags$row_lag,
                                  col_lag = lags$col_lag, acf = acf),
         convergence = object$convergence,
         optim_message = object$optim$message,
         inference = object$inference),
    class = "karma2d_diagnostics"
  )
}

#' @export
print.karma2d_diagnostics <- function(x, digits = 4L, ...) {
  cat("2-D KARMA fit diagnostics (quantile residuals, effective region)\n\n")
  print(round(x$residual_moments, digits))
  cat("\nSpatial autocorrelation of residuals:\n")
  print(x$spatial_acf, digits = digits, row.names = FALSE)
  cat(sprintf("\nBFGS convergence code: %d\n", x$convergence))
  cat(sprintf("Information matrix (%s): %s\n", x$inference$type,
              x$inference$status))
  invisible(x)
}
