#' Methods for 2-D KARMA fits
#'
#' Standard accessors for objects returned by [karma2d_fit()].
#' `confint()` works through [stats::confint.default()] and gives Wald
#' intervals \eqn{\hat\gamma_i \pm z_{1-\alpha/2}\,se(\hat\gamma_i)}; `AIC()`
#' and `BIC()` work through `logLik()`.
#'
#' @param object,x A `karma2d_fit` (or `summary.karma2d_fit`) object.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` an object of
#'   class `logLik` with attributes `df` and `nobs`; `nobs()` the number of
#'   pixels in the likelihood; `summary()` an object of class
#'   `summary.karma2d_fit` holding the coefficient table.
#' @name karma2d_fit-methods
#' @examples
#' y <- karma2d_simulate(30, 30, alpha = 1.11, ma = c(0.35, 0.28, 0.18),
#'                       shape = 12.4, seed = 5)
#' fit <- karma2d_fit(y, order = c(0, 1))
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' logLik(fit)
#' AIC(fit)
#' summary(fit)
NULL

#' @rdname karma2d_fit-methods
#' @export
coef.karma2d_fit <- function(object, ...) object$coefficients

#' @rdname karma2d_fit-methods
#' @export
vcov.karma2d_fit <- function(object, ...) object$vcov

#' @rdname karma2d_fit-methods
#' @export
logLik.karma2d_fit <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients),
            nobs = object$nobs, class = "logLik")
}

#' @rdname karma2d_fit-methods
#' @export
nobs.karma2d_fit <- function(object, ...) object$nobs

model_label <- function(object) {
  sprintf("2-D KARMA(%d, %d)", object$order[["p"]], object$order[["q"]])
}

#' @rdname karma2d_fit-methods
#' @export
print.karma2d_fit <- function(x, digits = max(3L, getOption("digits") - 3L),
                              ...) {
  cat(model_label(x), "model with", x$link, "link\n")
  cat(sprintf("Image: %d x %d, effective region: %d pixels\n",
              nrow(x$data), ncol(x$data), x$nobs))
  cat("\nCoefficients:\n")
  print(format(x$coefficients, digits = digits), quote = FALSE)
  cat(sprintf("\nLog-likelihood: %s", format(x$loglik, digits = digits)))
  if (x$convergence != 0L) {
    cat(sprintf("\nWarning: BFGS did not converge (code %d).",
                x$convergence))
  }
  cat("\n")
  invisible(x)
}

#' @rdname karma2d_fit-methods
#' @export
summary.karma2d_fit <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  table <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(
    list(call = object$call, label = model_label(object), link = object$link,
         coefficients = table, loglik = object$loglik, nobs = object$nobs,
         aic = stats::AIC(object), bic = stats::BIC(object),
         convergence = object$convergence, inference = object$inference,
         residuals = as.vector(object$residuals)),
    class = "summary.karma2d_fit"
  )
}

#' @rdname karma2d_fit-methods
#' @export
print.summary.karma2d_fit <- function(x,
                                      digits = max(3L,
                                                   getOption("digits") - 3L),
                                      ...) {
  cat(x$label, "model with", x$link, "link\n\n")
  cat("Call:\n")
  print(x$call)
  r <- x$residuals[is.finite(x$residuals)]
  cat("\nQuantile residuals:\n")
  print(stats::setNames(stats::quantile(r), c("Min", "1Q", "Median", "3Q",
                                              "Max")), digits = digits)
  cat(sprintf("\nCoefficients (%s information):\n", x$inference$type))
  stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = TRUE,
                      P.values = TRUE)
  cat(sprintf("\nLog-likelihood: %s on %d parameters, %d pixels\n",
              format(x$loglik, digits = digits),
              nrow(x$coefficients), x$nobs))
  cat(sprintf("AIC: %s   BIC: %s\n", format(x$aic, digits = digits),
              format(x$bic, digits = digits)))
  if (x$convergence != 0L) {
    cat(sprintf("WARNING: BFGS did not converge (code %d).\n",
                x$convergence))
  }
  if (x$inference$status != "ok") {
    cat("WARNING:", x$inference$message,
        "Wald inference should not be trusted.\n")
  }
  invisible(x)
}

#' Log-likelihood, score vector and information matrix of the 2-D KARMA model
#'
#' Low-level functions evaluating the conditional log-likelihood (Eq. 5 of
#' the article), its analytic gradient (Appendix A), and the conditional
#' expected (Fisher) or observed information at arbitrary parameter values.
#' They are mainly useful for validation and for research extensions.
#'
#' The expected information is the conditional Fisher information
#' \eqn{K = \sum E[U_{nm} U_{nm}^\top \mid S[n,m]]} summed over the effective
#' region, where \eqn{U_{nm}} is the pixel score. Its \eqn{\gamma\varphi} and
#' \eqn{\varphi\varphi} entries differ from Appendix B of the article, which
#' contains sign errors; see `vignette("model-theory")`. The observed
#' information is minus the Hessian of the log-likelihood, computed by central
#' differences of the analytic score.
#'
#' @param y Numeric matrix with values in `(lower, upper)`.
#' @param parameters Parameter vector `c(alpha, phi, theta, shape)` (see
#'   [karma2d_lags()] for the ordering).
#' @param order Integer vector `c(p, q)`.
#' @param link,lower,upper As in [karma2d_fit()].
#' @param type `"expected"` or `"observed"`.
#' @return `karma2d_loglik()` a number (`-Inf` at invalid parameters),
#'   `karma2d_score()` a named vector, `karma2d_information()` a named
#'   symmetric matrix.
#' @examples
#' y <- karma2d_simulate(20, 20, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 1)
#' par <- c(0.66, 0.28, 0.19, 0.17, 19.6)
#' karma2d_loglik(y, par, order = c(1, 0))
#' karma2d_score(y, par, order = c(1, 0))
#' karma2d_information(y, par, order = c(1, 0))
#' @export
karma2d_loglik <- function(y, parameters, order, link = "logit", lower = 0,
                           upper = 1) {
  args <- prepare_evaluation(y, parameters, order, link, lower, upper)
  karma2d_loglik_internal(args$par, args$z, args$z_link, args$order,
                          args$link_obj)
}

#' @rdname karma2d_loglik
#' @export
karma2d_score <- function(y, parameters, order, link = "logit", lower = 0,
                          upper = 1) {
  args <- prepare_evaluation(y, parameters, order, link, lower, upper)
  stats::setNames(
    karma2d_score_internal(args$par, args$z, args$z_link, args$order,
                           args$link_obj),
    karma2d_coef_names(args$order)
  )
}

#' @rdname karma2d_loglik
#' @export
karma2d_information <- function(y, parameters, order,
                                type = c("expected", "observed"),
                                link = "logit", lower = 0, upper = 1) {
  type <- match.arg(type)
  args <- prepare_evaluation(y, parameters, order, link, lower, upper)
  fun <- if (type == "expected") karma2d_expected_information else
    karma2d_observed_information
  info <- fun(args$par, args$z, args$z_link, args$order, args$link_obj)
  nm <- karma2d_coef_names(args$order)
  dimnames(info) <- list(nm, nm)
  info
}

prepare_evaluation <- function(y, parameters, order, link, lower, upper) {
  order <- check_order(order)
  check_support(lower, upper)
  link_obj <- make_karma2d_link(link)
  z <- to_unit_scale(y, lower, upper)
  par <- as.numeric(parameters)
  split_parameters(par, order) # validates the length
  list(par = par, z = z, z_link = link_obj$linkfun(z), order = order,
       link_obj = link_obj)
}

#' Fisher information of a fitted 2-D KARMA model
#'
#' @param object A `karma2d_fit` object.
#' @param type `"expected"` (conditional Fisher information) or `"observed"`.
#' @return Named symmetric matrix evaluated at the estimates.
#' @examples
#' y <- karma2d_simulate(25, 25, alpha = 0.66, ar = c(0.28, 0.19, 0.17),
#'                       shape = 19.6, seed = 9)
#' fit <- karma2d_fit(y, order = c(1, 0))
#' karma2d_fisher(fit)
#' @export
karma2d_fisher <- function(object, type = c("expected", "observed")) {
  stopifnot(inherits(object, "karma2d_fit"))
  type <- match.arg(type)
  karma2d_information(object$data, object$coefficients, object$order,
                      type = type, link = object$link,
                      lower = object$support[["lower"]],
                      upper = object$support[["upper"]])
}
