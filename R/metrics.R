#' Prediction accuracy metrics
#'
#' Computes the four image quality measures used in the article, comparing
#' observed values \eqn{y_i} with predictions \eqn{\hat y_i}, \eqn{i = 1,
#' \ldots, n}:
#' \describe{
#' \item{MAE}{\eqn{n^{-1}\sum |y_i - \hat y_i|};}
#' \item{MSE}{\eqn{n^{-1}\sum (y_i - \hat y_i)^2};}
#' \item{SMAPE}{\eqn{100\, n^{-1}\sum 2|y_i - \hat y_i| / (|y_i| +
#'   |\hat y_i|)}, in percent, with the convention 0/0 = 0;}
#' \item{RMSLE}{\eqn{\{n^{-1}\sum (\log(1 + y_i) - \log(1 + \hat
#'   y_i))^2\}^{1/2}}.}
#' }
#' Both inputs must be on the same scale; the article computes the metrics on
#' the rescaled unit interval. Pixels where either value is `NA` (e.g. the
#' unpredicted border) are dropped.
#'
#' @param observed,predicted Numeric vectors or matrices of equal size.
#' @return Named numeric vector `c(MAE, MSE, SMAPE, RMSLE)`, with attribute
#'   `n` (number of pixels used).
#' @examples
#' karma2d_metrics(c(0.2, 0.5, 0.9), c(0.25, 0.45, 0.8))
#' @export
karma2d_metrics <- function(observed, predicted) {
  if (!is.numeric(observed) || !is.numeric(predicted) ||
      length(observed) != length(predicted)) {
    stop("`observed` and `predicted` must be numeric with the same length.",
         call. = FALSE)
  }
  if (!is.null(dim(observed)) && !is.null(dim(predicted)) &&
      !identical(dim(observed), dim(predicted))) {
    stop("`observed` and `predicted` must have the same dimensions.",
         call. = FALSE)
  }
  ok <- !is.na(observed) & !is.na(predicted)
  y <- as.vector(observed)[ok]
  y_hat <- as.vector(predicted)[ok]
  if (length(y) == 0L) stop("No pixel has both values available.",
                            call. = FALSE)
  if (any(y <= -1 | y_hat <= -1)) {
    stop("RMSLE requires values greater than -1.", call. = FALSE)
  }
  abs_err <- abs(y - y_hat)
  denom <- abs(y) + abs(y_hat)
  smape_terms <- ifelse(denom == 0, 0, 2 * abs_err / denom)
  structure(
    c(MAE = mean(abs_err), MSE = mean((y - y_hat)^2),
      SMAPE = 100 * mean(smape_terms),
      RMSLE = sqrt(mean((log1p(y) - log1p(y_hat))^2))),
    n = length(y)
  )
}
