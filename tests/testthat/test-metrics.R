test_that("MAE and MSE follow their definitions", {
  m <- karma2d_metrics(c(0.2, 0.5, 0.9), c(0.3, 0.5, 0.7))
  expect_equal(unname(m[c("MAE", "MSE")]), c(0.1, (0.01 + 0.04) / 3))
})

test_that("SMAPE follows its definition", {
  m <- karma2d_metrics(c(0.2, 0.5), c(0.3, 0.5))
  expect_equal(m[["SMAPE"]], 100 * mean(c(2 * 0.1 / 0.5, 0)))
})

test_that("SMAPE treats 0 / 0 as zero", {
  expect_equal(karma2d_metrics(c(0, 0.5), c(0, 0.5))[["SMAPE"]], 0)
})

test_that("RMSLE follows its definition", {
  m <- karma2d_metrics(c(0.2, 0.5), c(0.3, 0.4))
  expect_equal(m[["RMSLE"]],
               sqrt(mean((log1p(c(0.2, 0.5)) - log1p(c(0.3, 0.4)))^2)))
})

test_that("metrics are zero for a perfect prediction", {
  expect_equal(unname(karma2d_metrics(c(0.1, 0.9), c(0.1, 0.9))),
               c(0, 0, 0, 0), ignore_attr = TRUE)
})

test_that("pixels with NA are dropped", {
  expect_equal(attr(karma2d_metrics(c(NA, 0.5, 0.6), c(0.1, NA, 0.6)), "n"),
               1L)
})

test_that("mismatched dimensions are rejected", {
  expect_error(karma2d_metrics(matrix(0.5, 2, 3), matrix(0.5, 3, 2)),
               "dimensions")
})
