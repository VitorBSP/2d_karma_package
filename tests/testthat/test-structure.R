test_that("an order-p model has (p + 1)^2 - 1 lags", {
  expect_equal(vapply(0:4, function(p) nrow(karma2d_lags(p)), integer(1)),
               (0:4 + 1L)^2 - 1L)
})

test_that("lag (0, 0) is excluded", {
  lags <- karma2d_lags(3)
  expect_false(any(lags$row_lag == 0 & lags$col_lag == 0))
})

test_that("first-order lags follow the article ordering", {
  expect_equal(karma2d_lags(1)$name, c("phi(0,1)", "phi(1,0)", "phi(1,1)"))
})

test_that("second-order lags vary the column lag fastest", {
  lags <- karma2d_lags(2)
  expect_equal(paste(lags$row_lag, lags$col_lag),
               c("0 1", "0 2", "1 0", "1 1", "1 2", "2 0", "2 1", "2 2"))
})

test_that("coefficient names cover alpha, AR, MA and shape", {
  expect_equal(karma2d:::karma2d_coef_names(c(1L, 1L)),
               c("alpha", "phi(0,1)", "phi(1,0)", "phi(1,1)",
                 "theta(0,1)", "theta(1,0)", "theta(1,1)", "shape"))
})

test_that("a wrong number of AR coefficients is rejected", {
  expect_error(karma2d_simulate(10, 10, alpha = 0, ar = c(0.1, 0.2),
                                shape = 5), "\\(p \\+ 1\\)\\^2 - 1")
})

test_that("phi(i, j) multiplies g(y[n - i, m - j])", {
  # With a single non-zero coefficient the predictor at [n, m] is
  # alpha + phi * g(y[n - i, m - j]).
  set.seed(10)
  y <- matrix(stats::runif(36, 0.2, 0.8), 6, 6)
  g <- stats::qlogis(y)
  ar <- c(0, 0.5, 0) # phi(1,0): previous row, same column
  rec <- karma2d:::karma2d_recursion_r(g, 0.1, ar, numeric(0), c(1L, 0L))
  expect_equal(rec$eta[4, 5], 0.1 + 0.5 * g[3, 5], tolerance = 1e-14)
})

test_that("theta(k, l) multiplies e[n - k, m - l]", {
  set.seed(11)
  y <- matrix(stats::runif(36, 0.2, 0.8), 6, 6)
  g <- stats::qlogis(y)
  ma <- c(0.4, 0, 0) # theta(0,1): same row, previous column
  rec <- karma2d:::karma2d_recursion_r(g, 0.2, numeric(0), ma, c(0L, 1L))
  expect_equal(rec$eta[3, 4], 0.2 + 0.4 * rec$error[3, 3],
               tolerance = 1e-14)
})

test_that("errors are g(y) minus the predictor", {
  set.seed(12)
  y <- matrix(stats::runif(25, 0.2, 0.8), 5, 5)
  g <- stats::qlogis(y)
  rec <- karma2d:::karma2d_recursion_r(g, 0.1, c(0.2, 0.1, 0.05),
                                       c(0.1, 0.1, 0.1), c(1L, 1L))
  expect_equal(rec$error[2:5, 2:5], g[2:5, 2:5] - rec$eta[2:5, 2:5],
               tolerance = 1e-14)
})

test_that("border pixels are not modelled", {
  set.seed(13)
  g <- matrix(stats::rnorm(49), 7, 7)
  rec <- karma2d:::karma2d_recursion_r(g, 0.3, numeric(0),
                                       stats::runif(8, -0.1, 0.1), c(0L, 2L))
  expect_true(all(rec$eta[1:2, ] == 0) && all(rec$eta[, 1:2] == 0))
})

test_that("pure AR predictor is linear in the lagged values", {
  set.seed(14)
  y <- matrix(stats::runif(64, 0.1, 0.9), 8, 8)
  g <- stats::qlogis(y)
  ar <- c(0.2, -0.1, 0.3)
  rec <- karma2d:::karma2d_recursion_r(g, 0.5, ar, numeric(0), c(1L, 0L))
  manual <- 0.5 + 0.2 * g[2:8, 1:7] - 0.1 * g[1:7, 2:8] + 0.3 * g[1:7, 1:7]
  expect_equal(rec$eta[2:8, 2:8], manual, tolerance = 1e-14)
})
