#' Control parameters for fitting a 2-D KARMA model
#'
#' @param shape_transform How the shape parameter is handled during
#'   optimization. `"log"` (default) optimizes over \eqn{\log\varphi}, so that
#'   every trial value is valid; the score is transformed by the chain rule
#'   and the estimates are reported on the original scale. `"identity"`
#'   optimizes over \eqn{\varphi} directly, as in the original scripts of the
#'   article; invalid trial values then get an infinite objective, which BFGS
#'   handles by step reduction.
#' @param maxit Maximum number of BFGS iterations (passed to [stats::optim()]).
#' @param reltol Relative convergence tolerance (passed to [stats::optim()]).
#' @param information Information matrix used for the covariance matrix:
#'   `"expected"` (conditional Fisher information, default) or `"observed"`
#'   (minus the Hessian of the log-likelihood, by central differences of the
#'   analytic score).
#' @param trace Non-negative integer passed to [stats::optim()].
#' @return A list of class `karma2d_control`.
#' @examples
#' karma2d_control(maxit = 1000)
#' @export
karma2d_control <- function(shape_transform = c("log", "identity"),
                            maxit = 500L, reltol = sqrt(.Machine$double.eps),
                            information = c("expected", "observed"),
                            trace = 0L) {
  shape_transform <- match.arg(shape_transform)
  information <- match.arg(information)
  check_scalar_number(maxit, "maxit", positive = TRUE)
  check_scalar_number(reltol, "reltol", positive = TRUE)
  structure(
    list(shape_transform = shape_transform, maxit = as.integer(maxit),
         reltol = reltol, information = information,
         trace = as.integer(trace)),
    class = "karma2d_control"
  )
}

#' Starting values for the 2-D KARMA optimizer
#'
#' Implements the initialization of Pereira, Palm and Bayer (2026, Sec. 2.3):
#' \eqn{\alpha_0 = g(\bar y)}, all autoregressive and moving average
#' coefficients equal to zero, and
#' \eqn{\varphi_0 = 0.3\{\bar y (1 - \bar y)/S^2 - 1\}}, where \eqn{\bar y}
#' and \eqn{S^2} are the mean and the sample variance (denominator
#' \eqn{NM - 1}) of all pixels, on the unit scale.
#'
#' @param y Numeric matrix with values in `(lower, upper)`.
#' @param order Integer vector `c(p, q)`.
#' @param link Link function name.
#' @param lower,upper Bounds of the support.
#' @return Named numeric vector of starting values.
#' @examples
#' y <- karma2d_simulate(20, 20, alpha = 0.5, ar = c(0.3, 0.2, 0.1),
#'                       shape = 15, seed = 1)
#' karma2d_start(y, order = c(1, 0))
#' @export
karma2d_start <- function(y, order, link = "logit", lower = 0, upper = 1) {
  order <- check_order(order)
  check_support(lower, upper)
  link_obj <- make_karma2d_link(link)
  z <- to_unit_scale(y, lower, upper)
  y_bar <- mean(z)
  s2 <- stats::var(as.vector(z))
  shape0 <- 0.3 * (y_bar * (1 - y_bar) / s2 - 1)
  if (!is.finite(shape0) || shape0 <= 0) {
    warning("Moment-based starting value for `shape` is not positive; ",
            "using 1 instead.", call. = FALSE)
    shape0 <- 1
  }
  start <- c(link_obj$linkfun(y_bar),
             rep(0, n_lag_coef(order[1]) + n_lag_coef(order[2])),
             shape0)
  stats::setNames(start, karma2d_coef_names(order))
}

#' Fit a 2-D KARMA model by conditional maximum likelihood
#'
#' Fits the two-dimensional Kumaraswamy autoregressive moving average model
#' of Pereira, Palm and Bayer (2026) to a bounded image. Conditionally on the
#' strongly causal neighbourhood, each pixel follows a Kumaraswamy
#' distribution with median \eqn{\mu[n,m]} and shape \eqn{\varphi}, where
#' \deqn{g(\mu[n,m]) = \eta[n,m] = \alpha +
#'   \sum_{(i,j)} \phi_{(i,j)} g(y[n-i, m-j]) +
#'   \sum_{(k,l)} \theta_{(k,l)} e[n-k, m-l],}
#' \eqn{e[n,m] = g(y[n,m]) - \eta[n,m]}, and the sums run over
#' \eqn{\{0,\ldots,p\}^2 \setminus (0,0)} and
#' \eqn{\{0,\ldots,q\}^2 \setminus (0,0)} (see [karma2d_lags()] for the
#' ordering of the coefficients).
#'
#' The conditional log-likelihood sums over pixels with row and column index
#' greater than \eqn{w = \max(p, q)}; the first \eqn{w} rows and columns are
#' treated as fixed, with errors set to zero. Maximization uses
#' [stats::optim()] with method `"BFGS"` and the analytic score vector of
#' Appendix A of the article. Starting values follow [karma2d_start()].
#'
#' Standard errors are derived from the inverse of the conditional Fisher
#' information evaluated at the estimates (or of the observed information,
#' see [karma2d_control()]). The expected information implemented here was
#' re-derived and checked against numerical integration; it corrects sign
#' errors in Appendix B of the article (see `vignette("model-theory")`). If
#' the information matrix is not numerically positive definite, the
#' covariance matrix is set to `NA` and a warning is issued; no
#' regularization is ever applied silently.
#'
#' @param y Numeric matrix (image) with all values strictly inside
#'   `(lower, upper)`. Rows index the first spatial coordinate `n`.
#' @param order Integer vector `c(p, q)` of AR and MA orders. `c(0, 0)` gives
#'   independent Kumaraswamy pixels with a common median.
#' @param link Link function: `"logit"` (default), `"probit"`, `"cloglog"` or
#'   `"cauchit"`.
#' @param lower,upper Bounds of the support of `y`. The model is fitted to
#'   `(y - lower) / (upper - lower)`.
#' @param start Optional starting values on the original parameter scale, in
#'   the order of [karma2d_start()]. Defaults to [karma2d_start()].
#' @param control List from [karma2d_control()].
#' @return An object of class `karma2d_fit`, a list with components:
#'   \describe{
#'   \item{`coefficients`}{Named estimates `alpha`, `phi(i,j)`, `theta(k,l)`,
#'     `shape`.}
#'   \item{`vcov`}{Estimated covariance matrix, or a matrix of `NA` when the
#'     information matrix is not positive definite.}
#'   \item{`information`}{Information matrix at the estimates.}
#'   \item{`inference`}{List with `status` (`"ok"`, `"ill-conditioned"` or
#'     `"unavailable"`), `message`, `type` and `rcond`.}
#'   \item{`loglik`, `nobs`}{Maximized conditional log-likelihood and number
#'     of pixels in the effective region.}
#'   \item{`order`, `link`, `support`}{Model specification.}
#'   \item{`fitted_median`}{Matrix of fitted conditional medians on the
#'     original scale (`NA` on the border).}
#'   \item{`residuals`}{Matrix of quantile residuals (`NA` on the border).}
#'   \item{`data`}{The input image (needed by `predict()` and refits).}
#'   \item{`start`, `control`, `optim`, `convergence`}{Optimization details:
#'     starting values, control settings, the raw [stats::optim()] output on
#'     the optimization scale, and its convergence code.}
#'   \item{`effective_region`}{Rows and columns entering the likelihood.}
#'   \item{`call`}{The matched call.}
#'   }
#' @references Pereira, V. B. S., Palm, B. G. and Bayer, F. M. (2026). A
#'   Kumaraswamy-based spatial model for bounded image data. Manuscript.
#' @seealso [karma2d_simulate()], [predict.karma2d_fit()],
#'   [karma2d_residuals()], [karma2d_detect_anomalies()].
#' @examples
#' y <- karma2d_simulate(30, 30, alpha = 0.66,
#'                       ar = c(0.28, 0.19, 0.17), shape = 19.6, seed = 42)
#' fit <- karma2d_fit(y, order = c(1, 0))
#' fit
#' summary(fit)
#' confint(fit)
#' @export
karma2d_fit <- function(y, order = c(1L, 0L), link = "logit", lower = 0,
                        upper = 1, start = NULL,
                        control = karma2d_control()) {
  call <- match.call()
  order <- check_order(order)
  check_support(lower, upper)
  if (!inherits(control, "karma2d_control")) {
    control <- do.call(karma2d_control, as.list(control))
  }
  link_obj <- make_karma2d_link(link)
  z <- to_unit_scale(y, lower, upper)
  w <- max(order)
  if (nrow(z) <= w + 1L || ncol(z) <= w + 1L) {
    stop("The image is too small for the requested order.", call. = FALSE)
  }
  coef_names <- karma2d_coef_names(order)
  n_par <- length(coef_names)
  if (is.null(start)) {
    start <- karma2d_start(y, order, link, lower, upper)
  } else {
    if (!is.numeric(start) || length(start) != n_par || anyNA(start)) {
      stop(sprintf("`start` must be a numeric vector of length %d.", n_par),
           call. = FALSE)
    }
    if (start[n_par] <= 0) stop("The starting shape must be positive.",
                                call. = FALSE)
    start <- stats::setNames(as.numeric(start), coef_names)
  }
  z_link <- link_obj$linkfun(z)

  log_shape <- control$shape_transform == "log"
  to_par <- function(theta) {
    if (log_shape) theta[n_par] <- exp(theta[n_par])
    theta
  }
  objective <- function(theta) {
    value <- -karma2d_loglik_internal(to_par(theta), z, z_link, order,
                                      link_obj)
    if (is.finite(value)) value else Inf
  }
  gradient <- function(theta) {
    par <- to_par(theta)
    score <- karma2d_score_internal(par, z, z_link, order, link_obj)
    if (log_shape) score[n_par] <- score[n_par] * par[n_par]
    -score
  }
  theta0 <- unname(start)
  if (log_shape) theta0[n_par] <- log(theta0[n_par])
  if (!is.finite(objective(theta0))) {
    stop("The log-likelihood is not finite at the starting values.",
         call. = FALSE)
  }
  opt <- stats::optim(
    theta0, objective, gradient, method = "BFGS",
    control = list(maxit = control$maxit, reltol = control$reltol,
                   trace = control$trace)
  )
  estimate <- stats::setNames(to_par(opt$par), coef_names)
  if (opt$convergence != 0L) {
    warning(sprintf("BFGS did not converge (code %d%s).", opt$convergence,
                    if (is.null(opt$message)) "" else
                      paste0(": ", opt$message)),
            call. = FALSE)
  }

  inference <- karma2d_inference(estimate, z, z_link, order, link_obj,
                                 control$information)
  filtered <- karma2d_recursion(z_link, estimate[1],
                                split_parameters(estimate, order)$ar,
                                split_parameters(estimate, order)$ma, order)
  region <- effective_region(nrow(z), ncol(z), w)
  mu <- link_obj$linkinv(filtered$eta)
  mu[!region] <- NA_real_
  resid <- quantile_residuals(z, mu, estimate[n_par])

  structure(
    list(
      coefficients = estimate,
      vcov = inference$vcov,
      information = inference$information,
      inference = inference[c("status", "message", "type", "rcond")],
      loglik = -opt$value,
      nobs = sum(region),
      order = stats::setNames(order, c("p", "q")),
      link = link,
      support = c(lower = lower, upper = upper),
      fitted_median = lower + (upper - lower) * mu,
      residuals = resid,
      data = y,
      start = start,
      control = control,
      optim = opt,
      convergence = opt$convergence,
      effective_region = list(rows = (w + 1L):nrow(z),
                              cols = (w + 1L):ncol(z)),
      call = call
    ),
    class = "karma2d_fit"
  )
}

# Information matrix, covariance matrix and a diagnosis of their reliability.
karma2d_inference <- function(par, z, z_link, order, link_obj, type) {
  info <- if (type == "expected") {
    karma2d_expected_information(par, z, z_link, order, link_obj)
  } else {
    karma2d_observed_information(par, z, z_link, order, link_obj)
  }
  nm <- names(par)
  dimnames(info) <- list(nm, nm)
  vc <- matrix(NA_real_, length(par), length(par), dimnames = list(nm, nm))
  rc <- NA_real_
  status <- "unavailable"
  message <- NULL
  if (any(!is.finite(info))) {
    message <- "The information matrix has non-finite entries."
  } else {
    chol_info <- tryCatch(chol(info), error = function(e) NULL)
    if (is.null(chol_info)) {
      message <- "The information matrix is not positive definite."
    } else {
      rc <- rcond(info)
      vc[] <- chol2inv(chol_info)
      if (rc < 1e-12) {
        status <- "ill-conditioned"
        message <- sprintf(paste0(
          "The information matrix is ill-conditioned (reciprocal condition ",
          "number %.2g); standard errors may be numerically unreliable."
        ), rc)
      } else {
        status <- "ok"
      }
    }
  }
  if (!is.null(message)) {
    warning(message, if (status == "unavailable")
      " Standard errors are not available." else "", call. = FALSE)
  }
  list(information = info, vcov = vc, status = status,
       message = message, type = type, rcond = rc)
}
