#' Simulate an image from a 2-D KARMA model
#'
#' Generates an \eqn{N \times M} image from the 2-D KARMA model with the given
#' parameters, by the inversion method (each pixel is a Kumaraswamy quantile
#' of a uniform variate).
#'
#' **Generation scheme.** The simulation runs on an extended grid of
#' \eqn{(N + b) \times (M + b)} pixels, where \eqn{b} is `burn_in`. The first
#' \eqn{w = \max(p, q)} rows and columns of the grid are initialized with
#' \eqn{g(y) = \alpha} and \eqn{e = 0}; every other pixel is drawn
#' sequentially, row by row, from its conditional distribution; finally the
#' first \eqn{b} rows and columns are discarded. With the default
#' `burn_in = 2 * max(p, q)` this reproduces the original simulation code of
#' the article exactly (same random stream for the same seed). Larger values
#' reduce the influence of the initialization on the returned image.
#'
#' **Stability.** No condition for stationarity or stability of the 2-D
#' KARMA process is known (the article lists this as an open problem). For
#' parameters far from those used in the article the generated field may
#' drift towards the bounds of the support; pixels numerically equal to 0
#' or 1 make the simulation stop with an error.
#'
#' @param nrow,ncol Image dimensions \eqn{N} (rows) and \eqn{M} (columns).
#' @param alpha Intercept.
#' @param ar Autoregressive coefficients, of length \eqn{(p+1)^2 - 1}, in the
#'   order of [karma2d_lags()]. Use `numeric(0)` (default) for no AR part.
#' @param ma Moving average coefficients, of length \eqn{(q+1)^2 - 1}.
#' @param shape Shape parameter \eqn{\varphi > 0}.
#' @param link Link function name.
#' @param burn_in Number of leading rows and columns discarded; at least
#'   \eqn{w}. Defaults to \eqn{2w}.
#' @param seed Optional integer seed. When given, the simulation is
#'   reproducible and the global random number state is left unchanged.
#' @param lower,upper Support of the returned image; values are generated on
#'   the unit interval and mapped linearly to `(lower, upper)`.
#' @return A numeric matrix of dimension `nrow` x `ncol`, with attribute
#'   `"karma2d"` describing the generating model.
#' @examples
#' y <- karma2d_simulate(50, 50, alpha = 0.3, ar = c(0.3, -0.15, 0.23),
#'                       ma = c(0.15, -0.1, 0.25), shape = 17, seed = 2026)
#' range(y)
#' image(y, main = "2-D KARMA(1, 1)")
#' @export
karma2d_simulate <- function(nrow, ncol = nrow, alpha, ar = numeric(0),
                             ma = numeric(0), shape, link = "logit",
                             burn_in = NULL, seed = NULL, lower = 0,
                             upper = 1) {
  check_scalar_number(nrow, "nrow", positive = TRUE)
  check_scalar_number(ncol, "ncol", positive = TRUE)
  check_scalar_number(alpha, "alpha")
  check_scalar_number(shape, "shape", positive = TRUE)
  check_support(lower, upper)
  if (!is.numeric(ar) || any(!is.finite(ar))) stop("`ar` must be numeric.")
  if (!is.numeric(ma) || any(!is.finite(ma))) stop("`ma` must be numeric.")
  link_obj <- make_karma2d_link(link)
  order <- c(order_from_length(length(ar), "ar"),
             order_from_length(length(ma), "ma"))
  w <- max(order)
  if (is.null(burn_in)) burn_in <- 2L * w
  check_scalar_number(burn_in, "burn_in")
  if (burn_in < w) stop("`burn_in` must be at least max(p, q).",
                        call. = FALSE)
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  burn_in <- as.integer(burn_in)

  y_unit <- with_local_seed(seed, simulate_grid(
    nrow + burn_in, ncol + burn_in, alpha, ar, ma, shape, order, link_obj
  ))
  keep_r <- (burn_in + 1L):(nrow + burn_in)
  keep_c <- (burn_in + 1L):(ncol + burn_in)
  y <- lower + (upper - lower) * y_unit[keep_r, keep_c, drop = FALSE]
  attr(y, "karma2d") <- list(
    order = stats::setNames(order, c("p", "q")),
    parameters = stats::setNames(c(alpha, ar, ma, shape),
                                 karma2d_coef_names(order)),
    link = link, burn_in = burn_in, seed = seed,
    support = c(lower = lower, upper = upper)
  )
  y
}

simulate_grid <- function(n_rows, n_cols, alpha, ar, ma, shape, order,
                          link_obj) {
  w <- max(order)
  off_ar <- lag_offsets(order[1], n_rows)
  off_ma <- lag_offsets(order[2], n_rows)
  y_link <- rep(alpha, n_rows * n_cols)
  y <- numeric(n_rows * n_cols)
  error <- numeric(n_rows * n_cols)
  n_draw <- (n_rows - w) * (n_cols - w)
  u <- stats::runif(n_draw)
  draw <- 0L
  for (i in (w + 1L):n_rows) {
    for (j in (w + 1L):n_cols) {
      idx <- i + (j - 1L) * n_rows
      eta <- alpha + sum(ar * y_link[idx - off_ar]) +
        sum(ma * error[idx - off_ma])
      mu <- link_obj$linkinv(eta)
      draw <- draw + 1L
      y_ij <- qkuma(u[draw], mu, shape)
      if (!is.finite(y_ij) || y_ij <= 0 || y_ij >= 1) {
        stop(sprintf(paste0(
          "Simulated pixel [%d, %d] of the extended grid is numerically ",
          "equal to 0 or 1 (median %.3g); the parameters may define an ",
          "unstable process."
        ), i, j, mu), call. = FALSE)
      }
      y[idx] <- y_ij
      y_link[idx] <- link_obj$linkfun(y_ij)
      error[idx] <- y_link[idx] - eta
    }
  }
  matrix(y, n_rows, n_cols)
}
