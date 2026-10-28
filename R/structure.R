#' Spatial lags of the isotropic 2-D KARMA model
#'
#' Returns the spatial lags `(i, j)` associated with the autoregressive (or
#' moving average) coefficients of an isotropic model of order `order`, in the
#' order in which the coefficients are stored throughout the package.
#'
#' The strongly causal neighbourhood of pixel `[n, m]` (row `n`, column `m`)
#' contains the pixels `[n - i, m - j]` with `i, j = 0, ..., order` and
#' `(i, j) != (0, 0)`. Coefficient \eqn{\phi_{(i,j)}} multiplies
#' \eqn{g(y[n - i, m - j])}: `i` is the **row** lag and `j` the **column**
#' lag. Coefficients are ordered lexicographically with the row lag varying
#' slowest: `(0,1), (0,2), ..., (0,p), (1,0), (1,1), ..., (p,p)`, which is
#' the ordering \eqn{\phi = (\phi_{(0,1)}, \phi_{(0,2)}, \ldots,
#' \phi_{(p,p)})^\top} of Pereira, Palm and Bayer (2026). An order `p` model
#' therefore has \eqn{(p + 1)^2 - 1} coefficients.
#'
#' @param order Non-negative integer, the AR order `p` or MA order `q`.
#' @param prefix Prefix of the coefficient names.
#' @return A data frame with columns `row_lag`, `col_lag` and `name`, with
#'   \eqn{(order + 1)^2 - 1} rows.
#' @examples
#' karma2d_lags(1)
#' karma2d_lags(2, prefix = "theta")
#' @export
karma2d_lags <- function(order, prefix = "phi") {
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
        order < 0 || order != round(order)) {
    stop("`order` must be a single non-negative integer.", call. = FALSE)
  }
  grid <- expand.grid(col_lag = seq_len(order + 1) - 1L,
                      row_lag = seq_len(order + 1) - 1L)
  grid <- grid[!(grid$row_lag == 0L & grid$col_lag == 0L),
               c("row_lag", "col_lag")]
  rownames(grid) <- NULL
  grid$name <- sprintf("%s(%d,%d)", prefix, grid$row_lag, grid$col_lag)
  grid
}

n_lag_coef <- function(order) (order + 1L)^2 - 1L

# Order implied by the number of coefficients, or an error.
order_from_length <- function(n_coef, what) {
  order <- sqrt(n_coef + 1) - 1
  if (abs(order - round(order)) > 1e-8) {
    stop(sprintf(paste0(
      "`%s` has %d coefficients; an isotropic model of order p needs ",
      "(p + 1)^2 - 1 coefficients (3, 8, 15, ...)."
    ), what, n_coef), call. = FALSE)
  }
  as.integer(round(order))
}

karma2d_coef_names <- function(order) {
  c("alpha",
    karma2d_lags(order[1], "phi")$name,
    karma2d_lags(order[2], "theta")$name,
    "shape")
}

# Split a parameter vector c(alpha, phi, theta, shape).
split_parameters <- function(par, order) {
  p1 <- n_lag_coef(order[1])
  q1 <- n_lag_coef(order[2])
  if (length(par) != p1 + q1 + 2L) {
    stop(sprintf(
      "A 2-D KARMA(%d, %d) model has %d parameters; got %d.",
      order[1], order[2], p1 + q1 + 2L, length(par)
    ), call. = FALSE)
  }
  par <- unname(par)
  list(
    alpha = par[1],
    ar = par[seq_len(p1) + 1L],
    ma = par[seq_len(q1) + 1L + p1],
    shape = par[p1 + q1 + 2L]
  )
}

# Offsets, in column-major linear indexing of an n-row matrix, of the lagged
# pixels [i - row_lag, j - col_lag] relative to [i, j].
lag_offsets <- function(order, n_rows) {
  lags <- karma2d_lags(order)
  as.integer(lags$row_lag + lags$col_lag * n_rows)
}

# Logical matrix flagging the effective region used by the conditional
# likelihood: rows and columns w + 1, ..., where w = max(p, q).
effective_region <- function(n_rows, n_cols, w) {
  (row(matrix(0, n_rows, n_cols)) > w) & (col(matrix(0, n_rows, n_cols)) > w)
}

#' Recursive linear predictor (reference R implementation)
#'
#' Computes \eqn{\eta[n,m]}, the errors \eqn{e[n,m] = g(y[n,m]) - \eta[n,m]}
#' and, optionally, the derivatives \eqn{\partial\eta/\partial\gamma} for
#' \eqn{\gamma = (\alpha, \phi^\top, \theta^\top)^\top} (Appendix A of the
#' article). Pixels are visited row by row; pixels in the first `w` rows or
#' columns are not modelled: their predictor, error and derivatives are
#' zero, i.e. the recursion is conditional on the border.
#'
#' @param y_link Matrix \eqn{g(y)}.
#' @param derivatives Logical; compute derivatives?
#' @return List with `eta`, `error` (matrices) and `derivatives`
#'   (`(1 + p1 + q1) x (n * k)` matrix, columns in column-major pixel order,
#'   or `NULL`).
#' @noRd
karma2d_recursion_r <- function(y_link, alpha, ar, ma, order,
                                derivatives = FALSE) {
  n_rows <- nrow(y_link)
  n_cols <- ncol(y_link)
  p <- order[1]
  q <- order[2]
  w <- max(p, q)
  off_ar <- lag_offsets(p, n_rows)
  off_ma <- lag_offsets(q, n_rows)
  y_vec <- as.vector(y_link)
  eta <- numeric(n_rows * n_cols)
  error <- numeric(n_rows * n_cols)
  n_par <- 1L + length(off_ar) + length(off_ma)
  deriv <- if (derivatives) matrix(0, n_par, n_rows * n_cols) else NULL

  if (n_rows > w && n_cols > w) {
    for (i in (w + 1L):n_rows) {
      for (j in (w + 1L):n_cols) {
        idx <- i + (j - 1L) * n_rows
        ar_lag <- y_vec[idx - off_ar]
        ma_lag <- error[idx - off_ma]
        eta_ij <- alpha + sum(ar * ar_lag) + sum(ma * ma_lag)
        eta[idx] <- eta_ij
        error[idx] <- y_vec[idx] - eta_ij
        if (derivatives) {
          base <- c(1, ar_lag, ma_lag)
          if (q > 0L) {
            base <- base - drop(deriv[, idx - off_ma, drop = FALSE] %*% ma)
          }
          deriv[, idx] <- base
        }
      }
    }
  }
  list(
    eta = matrix(eta, n_rows, n_cols),
    error = matrix(error, n_rows, n_cols),
    derivatives = deriv
  )
}
