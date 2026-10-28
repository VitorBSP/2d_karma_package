y_full <- simulate_scenario(paper_scenarios$ar, 40, seed = 11)
fit <- karma2d_fit(y_full[21:40, 21:40], order = c(1, 0))

test_that("prediction on the fitted data equals the fitted medians", {
  expect_equal(predict(fit), fitted(fit))
})

test_that("prediction without padding leaves the border NA", {
  p <- predict(fit, newdata = y_full)
  expect_true(all(is.na(p[1, ])) && all(is.na(p[, 1])) && !anyNA(p[-1, -1]))
})

test_that("duplicate padding predicts the whole image", {
  expect_false(anyNA(predict(fit, newdata = y_full, padding = "duplicate")))
})

test_that("padding does not change interior predictions of a pure AR model", {
  # Without MA terms the predictor depends only on observed neighbours.
  a <- predict(fit, newdata = y_full)
  b <- predict(fit, newdata = y_full, padding = "duplicate")
  expect_equal(a[-1, -1], b[-1, -1], tolerance = 1e-14)
})

test_that("duplicate padding copies the first 2w rows and columns", {
  z <- matrix(1:16, 4, 4)
  p <- karma2d:::pad_duplicate(z, 1L)
  expect_equal(p, cbind(rbind(z[1:2, ], z)[, 1:2], rbind(z[1:2, ], z)))
})

test_that("predicted medians stay inside the support", {
  p <- predict(fit, newdata = y_full, padding = "duplicate")
  expect_true(all(p > 0 & p < 1))
})

test_that("link predictions are the link of median predictions", {
  expect_equal(predict(fit, type = "link"), stats::qlogis(predict(fit)),
               tolerance = 1e-12)
})

test_that("quantile residuals are Phi^-1 of the conditional cdf", {
  r <- residuals(fit)
  mu <- fitted(fit)
  manual <- stats::qnorm(pkuma(fit$data, mu, coef(fit)[["shape"]]))
  expect_equal(r[-1, -1], manual[-1, -1], tolerance = 1e-10)
})

test_that("quantile residuals stay finite for extreme observations", {
  r <- karma2d:::quantile_residuals(matrix(c(1e-12, 1 - 1e-12), 1),
                                    matrix(0.5, 1, 2), 5)
  expect_true(all(is.finite(r)))
})

test_that("quantile residuals are roughly standard normal under the model", {
  r <- as.vector(residuals(fit))
  r <- r[is.finite(r)]
  expect_equal(c(mean(r), stats::sd(r)), c(0, 1), tolerance = 0.15)
})

test_that("response residuals are y minus the fitted median", {
  expect_equal(residuals(fit, type = "response"), fit$data - fitted(fit))
})

test_that("link residuals equal the recursion errors", {
  rec <- karma2d:::karma2d_recursion_r(stats::qlogis(fit$data),
                                       coef(fit)[[1]], coef(fit)[2:4],
                                       numeric(0), c(1L, 0L))
  expect_equal(residuals(fit, type = "link")[-1, -1], rec$error[-1, -1],
               tolerance = 1e-10)
})

test_that("karma2d_residuals on the fitted data equals residuals()", {
  expect_equal(karma2d_residuals(fit), residuals(fit))
})

test_that("karma2d_residuals works on a larger image", {
  expect_equal(dim(karma2d_residuals(fit, y_full, padding = "duplicate")),
               dim(y_full))
})

test_that("diagnostics report residual moments", {
  d <- karma2d_diagnostics(fit)
  expect_named(d$residual_moments,
               c("mean", "sd", "skewness", "kurtosis", "prop_outside_3"))
})

test_that("diagnostics report three first-order spatial autocorrelations", {
  expect_equal(nrow(karma2d_diagnostics(fit)$spatial_acf), 3L)
})

test_that("diagnostics print", {
  expect_output(print(karma2d_diagnostics(fit)), "autocorrelation")
})

test_that("karma2d_fisher matches karma2d_information at the estimate", {
  expect_equal(karma2d_fisher(fit),
               karma2d_information(fit$data, coef(fit), c(1, 0)))
})
