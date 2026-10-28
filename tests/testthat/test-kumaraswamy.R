# Tolerances: quantities computed in closed form are compared at 1e-12
# (a few ulps after several transcendental operations); numerical
# integrals at 1e-8 (requested rel.tol of integrate() is 1e-10).

test_that("density integrates to one over a range of shapes and medians", {
  # Substituting y = 1 - t^k with k = max(1, 1 / delta) removes the
  # integrable singularity (1 - y)^(delta - 1) at y = 1 when delta < 1.
  # Cases with delta < 0.5 are excluded: there a non-negligible part of the
  # mass lies within 1e-15 of y = 1, beyond double precision resolution.
  grid <- expand.grid(median = c(0.05, 0.3, 0.5, 0.8, 0.97),
                      shape = c(0.5, 1, 3, 20))
  grid <- grid[log(0.5) / log1p(-grid$median^grid$shape) >= 0.5, ]
  total <- mapply(function(m, s) {
    k <- max(1, log1p(-m^s) / log(0.5))
    stats::integrate(function(t) dkuma(1 - t^k, m, s) * k * t^(k - 1),
                     0, 1, rel.tol = 1e-10)$value
  }, grid$median, grid$shape)
  expect_equal(total, rep(1, nrow(grid)), tolerance = 1e-8)
})

test_that("density integrates to one on a general support", {
  total <- stats::integrate(dkuma, -1, 1, median = 0.4, shape = 4,
                            lower = -1, upper = 1)$value
  expect_equal(total, 1, tolerance = 1e-8)
})

test_that("cdf equals the integral of the density", {
  x <- c(0.1, 0.4, 0.75)
  numeric_cdf <- vapply(x, function(b) {
    stats::integrate(dkuma, 0, b, median = 0.6, shape = 5,
                     rel.tol = 1e-11)$value
  }, numeric(1))
  expect_equal(pkuma(x, 0.6, 5), numeric_cdf, tolerance = 1e-8)
})

test_that("the median parameter is the median", {
  m <- c(0.01, 0.2, 0.5, 0.9, 0.999)
  expect_equal(pkuma(m, m, shape = c(0.3, 19.6, 3, 1, 50)), rep(0.5, 5),
               tolerance = 1e-12)
})

test_that("qkuma(0.5) returns the median even when median^shape is tiny", {
  expect_equal(qkuma(0.5, median = 0.2, shape = 19.6), 0.2,
               tolerance = 1e-12)
})

test_that("cdf is non-decreasing", {
  x <- seq(0, 1, length.out = 501)
  expect_true(all(diff(pkuma(x, 0.7, 2.5)) >= 0))
})

test_that("quantile function inverts the cdf", {
  p <- c(1e-10, 0.01, 0.3, 0.5, 0.9, 1 - 1e-10)
  expect_equal(pkuma(qkuma(p, 0.35, 7), 0.35, 7), p, tolerance = 1e-10)
})

test_that("cdf inverts the quantile function", {
  x <- c(0.05, 0.3, 0.6, 0.95)
  expect_equal(qkuma(pkuma(x, 0.8, 12), 0.8, 12), x, tolerance = 1e-10)
})

test_that("upper tail and log scale are consistent", {
  x <- c(0.2, 0.5, 0.9)
  expect_equal(pkuma(x, 0.5, 3, lower.tail = FALSE, log.p = TRUE),
               log1p(-pkuma(x, 0.5, 3)), tolerance = 1e-12)
})

test_that("extreme upper tail probabilities keep relative accuracy", {
  # 1 - F(x) = (1 - x^shape)^delta in closed form
  x <- 0.9999
  delta <- log(0.5) / log1p(-0.5^3)
  expect_equal(pkuma(x, 0.5, 3, lower.tail = FALSE),
               (-expm1(3 * log(x)))^delta, tolerance = 1e-10)
})

test_that("quantiles accept log probabilities and upper tails", {
  expect_equal(qkuma(log(0.25), 0.5, 3, lower.tail = FALSE, log.p = TRUE),
               qkuma(0.75, 0.5, 3), tolerance = 1e-12)
})

test_that("log density equals the log of the density", {
  x <- c(0.1, 0.5, 0.9)
  expect_equal(dkuma(x, 0.4, 2, log = TRUE), log(dkuma(x, 0.4, 2)),
               tolerance = 1e-12)
})

test_that("general support is a location-scale transform of the unit case", {
  x <- c(-0.8, 0, 0.7)
  expect_equal(dkuma(x, 0.2, 3, lower = -1, upper = 1),
               dkuma((x + 1) / 2, 0.6, 3) / 2, tolerance = 1e-12)
})

test_that("quantiles on a general support are on the original scale", {
  expect_equal(qkuma(0.5, median = 40, shape = 3, lower = 0, upper = 100),
               40, tolerance = 1e-12)
})

test_that("density is zero at and outside the bounds", {
  expect_equal(dkuma(c(-0.1, 0, 1, 1.5), 0.5, 2), c(0, 0, 0, 0))
})

test_that("cdf is 0 below and 1 above the support", {
  expect_equal(pkuma(c(-1, 0, 1, 2), 0.5, 2), c(0, 0, 1, 1))
})

test_that("quantiles at p = 0 and p = 1 are the bounds", {
  expect_equal(qkuma(c(0, 1), 3, 2, lower = 2, upper = 5), c(2, 5))
})

test_that("arguments are recycled", {
  expect_equal(dkuma(0.5, c(0.2, 0.5, 0.8), 3),
               c(dkuma(0.5, 0.2, 3), dkuma(0.5, 0.5, 3), dkuma(0.5, 0.8, 3)))
})

test_that("zero-length input gives zero-length output", {
  expect_length(pkuma(numeric(0), 0.5, 2), 0L)
})

test_that("missing values propagate", {
  expect_equal(is.na(dkuma(c(0.5, NA), c(NA, 0.5), 2)), c(TRUE, TRUE))
})

test_that("invalid shape gives NaN with a warning", {
  expect_warning(out <- dkuma(0.5, 0.5, -1), "NaN")
  expect_true(is.nan(out))
})

test_that("median outside the support gives NaN with a warning", {
  expect_warning(out <- pkuma(0.5, 1.2, 2), "NaN")
  expect_true(is.nan(out))
})

test_that("probabilities outside [0, 1] give NaN with a warning", {
  expect_warning(out <- qkuma(1.5, 0.5, 2), "NaN")
  expect_true(is.nan(out))
})

test_that("random deviates lie strictly inside the support", {
  set.seed(1)
  x <- rkuma(1e4, median = 0.95, shape = 30)
  expect_true(all(x > 0 & x < 1))
})

test_that("random deviates have the requested median", {
  set.seed(2)
  x <- rkuma(2e4, median = 0.3, shape = 4)
  # standard error of the sample median ~ 1 / (2 f(m) sqrt(n)) ~ 0.0013
  expect_equal(stats::median(x), 0.3, tolerance = 0.01)
})

test_that("random deviates follow the distribution (KS test)", {
  set.seed(3)
  x <- rkuma(5000, median = 0.6, shape = 2.5)
  p <- suppressWarnings(stats::ks.test(x, pkuma, median = 0.6,
                                       shape = 2.5)$p.value)
  expect_gt(p, 0.001)
})

test_that("rkuma with length(n) > 1 uses length(n)", {
  set.seed(4)
  expect_length(rkuma(c(1, 2, 3), 0.5, 2), 3L)
})
