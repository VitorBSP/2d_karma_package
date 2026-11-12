# Pad an image as in Pereira, Palm and Bayer (2026, Sec. 3.1): the first 2w
# rows are duplicated above the image, then the first 2w columns of the
# result are duplicated to its left, giving an (N + 2w) x (M + 2w) image.
# The article calls this "mirror padding"; the rows and columns are copied
# without reversing their order, exactly as in the original scripts.
pad_duplicate <- function(z, w) {
  if (w == 0L) return(z)
  if (nrow(z) < 2L * w || ncol(z) < 2L * w) {
    stop("Duplicate padding needs at least 2 * max(p, q) rows and columns.",
         call. = FALSE)
  }
  z <- rbind(z[seq_len(2L * w), , drop = FALSE], z)
  cbind(z[, seq_len(2L * w), drop = FALSE], z)
}

# Linear predictor on the unit-scale image z, with border handling.
karma2d_filter_eta <- function(z, par, order, link_obj, padding) {
  w <- max(order)
  parts <- split_parameters(par, order)
  if (padding == "duplicate" && w > 0L) {
    zp <- pad_duplicate(z, w)
    rec <- karma2d_recursion(link_obj$linkfun(zp), parts$alpha, parts$ar,
                             parts$ma, order)
    keep_r <- (2L * w + 1L):nrow(zp)
    keep_c <- (2L * w + 1L):ncol(zp)
    return(rec$eta[keep_r, keep_c, drop = FALSE])
  }
  rec <- karma2d_recursion(link_obj$linkfun(z), parts$alpha, parts$ar,
                           parts$ma, order)
  eta <- rec$eta
  eta[!effective_region(nrow(z), ncol(z), w)] <- NA_real_
  eta
}

#' Predicted (filtered) image from a 2-D KARMA fit
#'
#' Computes the one-step predicted conditional medians
#' \deqn{\hat\mu[n,m] = g^{-1}\Big(\hat\alpha + \sum \hat\phi_{(i,j)}
#'   g(y[n-i,m-j]) + \sum \hat\theta_{(k,l)} \hat e[n-k,m-l]\Big)}
#' (Eq. 8 of the article), with \eqn{\hat e[n,m] = g(y[n,m]) - g(\hat\mu[n,m])}.
#'
#' @param object A `karma2d_fit` object.
#' @param newdata Optional matrix on the original scale (e.g. the full scene
#'   when the model was fitted on a calibration window). Defaults to the
#'   data used in the fit.
#' @param type `"median"` (original scale, default) or `"link"`
#'   (\eqn{\hat\eta}).
#' @param padding `"none"`: the first \eqn{w = \max(p, q)} rows and columns
#'   are not predicted (`NA`), as in the conditional likelihood.
#'   `"duplicate"`: the image is first extended by duplicating its first
#'   \eqn{2w} rows (on top) and then its first \eqn{2w} columns (on the left);
#'   the recursion runs on the extended image and the full \eqn{N \times M}
#'   prediction is returned. This is the border strategy of Section 3.1 of
#'   the article (called "mirror padding" there, although rows and columns
#'   are copied without reversal).
#' @param ... Unused.
#' @return A matrix with the dimension of `newdata`.
#' @examples
#' y <- karma2d_simulate(30, 30, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 42)
#' fit <- karma2d_fit(y, order = c(1, 0))
#' mu_hat <- predict(fit, padding = "duplicate")
#' karma2d_metrics(y, mu_hat)
#' @export
predict.karma2d_fit <- function(object, newdata = NULL,
                                type = c("median", "link"),
                                padding = c("none", "duplicate"), ...) {
  type <- match.arg(type)
  padding <- match.arg(padding)
  if (is.null(newdata)) newdata <- object$data
  sup <- object$support
  z <- to_unit_scale(newdata, sup[["lower"]], sup[["upper"]], "newdata")
  link_obj <- make_karma2d_link(object$link)
  eta <- karma2d_filter_eta(z, object$coefficients, object$order, link_obj,
                            padding)
  if (type == "link") return(eta)
  sup[["lower"]] + (sup[["upper"]] - sup[["lower"]]) * link_obj$linkinv(eta)
}

#' @rdname predict.karma2d_fit
#' @export
fitted.karma2d_fit <- function(object, ...) {
  object$fitted_median
}
