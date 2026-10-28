# Quantile residuals r = Phi^{-1}(F(y | mu, shape)) (Dunn and Smyth, 1996),
# computed from whichever tail is smaller so that extreme residuals do not
# overflow to +/- Inf. y and mu on the unit scale; NA in mu gives NA.
quantile_residuals <- function(y, mu, shape) {
  out <- y
  out[] <- NA_real_
  ok <- !is.na(mu) & !is.na(y)
  log_lower <- pkuma(y[ok], mu[ok], shape, log.p = TRUE)
  log_upper <- pkuma(y[ok], mu[ok], shape, lower.tail = FALSE, log.p = TRUE)
  out[ok] <- ifelse(log_lower < log_upper,
                    stats::qnorm(log_lower, log.p = TRUE),
                    -stats::qnorm(log_upper, log.p = TRUE))
  out
}

#' Residuals of a 2-D KARMA fit
#'
#' @param object A `karma2d_fit` object.
#' @param type `"quantile"` (default): quantile residuals
#'   \eqn{r[n,m] = \Phi^{-1}(F(y[n,m] \mid \hat\mu[n,m], \hat\varphi))};
#'   `"response"`: \eqn{y - \hat\mu} on the original scale;
#'   `"link"`: \eqn{g(y) - \hat\eta}, the errors entering the MA part.
#' @param ... Unused.
#' @return A matrix of the same dimension as the data, with `NA` in the first
#'   \eqn{w = \max(p, q)} rows and columns.
#' @examples
#' y <- karma2d_simulate(25, 25, alpha = 1.1, ma = c(0.35, 0.28, 0.18),
#'                       shape = 12.4, seed = 3)
#' fit <- karma2d_fit(y, order = c(0, 1))
#' r <- residuals(fit)
#' summary(as.vector(r))
#' @export
residuals.karma2d_fit <- function(object, type = c("quantile", "response",
                                                   "link"), ...) {
  type <- match.arg(type)
  switch(type,
    quantile = object$residuals,
    response = object$data - object$fitted_median,
    link = {
      link_obj <- make_karma2d_link(object$link)
      sup <- object$support
      unit <- function(v) (v - sup[["lower"]]) / (sup[["upper"]] - sup[["lower"]])
      link_obj$linkfun(unit(object$data)) -
        link_obj$linkfun(unit(object$fitted_median))
    }
  )
}

#' Residuals of a fitted 2-D KARMA model on a new image
#'
#' Applies the fitted model to `newdata` (for instance, the full scene after
#' fitting on a Phase I calibration window) and returns the residuals.
#'
#' @param object A `karma2d_fit` object.
#' @param newdata Numeric matrix on the original scale. Defaults to the data
#'   used in the fit.
#' @param type Residual type, see [residuals.karma2d_fit()].
#' @param padding Border handling, see [predict.karma2d_fit()].
#' @return A matrix of residuals with the dimension of `newdata`.
#' @examples
#' y <- karma2d_simulate(40, 40, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 7)
#' fit <- karma2d_fit(y[21:40, 21:40], order = c(1, 0))
#' r <- karma2d_residuals(fit, newdata = y, padding = "duplicate")
#' dim(r)
#' @export
karma2d_residuals <- function(object, newdata = NULL,
                              type = c("quantile", "response", "link"),
                              padding = c("none", "duplicate")) {
  stopifnot(inherits(object, "karma2d_fit"))
  type <- match.arg(type)
  padding <- match.arg(padding)
  if (is.null(newdata)) newdata <- object$data
  mu <- stats::predict(object, newdata = newdata, padding = padding)
  sup <- object$support
  width <- sup[["upper"]] - sup[["lower"]]
  z <- to_unit_scale(newdata, sup[["lower"]], sup[["upper"]], "newdata")
  mu_unit <- (mu - sup[["lower"]]) / width
  switch(type,
    quantile = quantile_residuals(z, mu_unit,
                                  object$coefficients[["shape"]]),
    response = newdata - mu,
    link = {
      link_obj <- make_karma2d_link(object$link)
      link_obj$linkfun(z) - link_obj$linkfun(mu_unit)
    }
  )
}
