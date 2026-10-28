block_mask <- function() {
  m <- matrix(FALSE, 10, 10)
  m[3:6, 4:8] <- TRUE
  m
}

test_that("erosion by a 3 x 3 square shrinks a block by one pixel", {
  e <- morph_erode(block_mask(), karma2d_structuring_element(3))
  expected <- matrix(FALSE, 10, 10)
  expected[4:5, 5:7] <- TRUE
  expect_equal(e, expected)
})

test_that("dilation by a 3 x 3 square grows a block by one pixel", {
  d <- morph_dilate(block_mask(), karma2d_structuring_element(3))
  expected <- matrix(FALSE, 10, 10)
  expected[2:7, 3:9] <- TRUE
  expect_equal(d, expected)
})

test_that("opening removes isolated pixels", {
  m <- block_mask()
  m[9, 2] <- TRUE
  expect_false(morph_open(m, karma2d_structuring_element(3))[9, 2])
})

test_that("opening preserves a block larger than the element", {
  expect_equal(morph_open(block_mask(), karma2d_structuring_element(3)),
               block_mask())
})

test_that("closing fills a one-pixel hole", {
  m <- block_mask()
  m[4, 6] <- FALSE
  expect_true(morph_close(m, karma2d_structuring_element(3))[4, 6])
})

test_that("a 3 x 5 element erodes rows by 1 and columns by 2", {
  e <- morph_erode(block_mask(), karma2d_structuring_element(3, 5))
  expected <- matrix(FALSE, 10, 10)
  expected[4:5, 6] <- TRUE
  expect_equal(e, expected)
})

test_that("objects touching the border are not eroded by the border", {
  m <- matrix(FALSE, 6, 6)
  m[1:3, 1:3] <- TRUE
  expect_true(morph_erode(m, karma2d_structuring_element(3))[1, 1])
})

test_that("an empty structuring element is rejected", {
  expect_error(morph_erode(block_mask(), matrix(FALSE, 3, 3)), "empty")
})

test_that("thresholding flags residuals beyond the limits", {
  r <- matrix(c(0, 3.5, -3.2, 2.9), 2)
  d <- karma2d_detect_anomalies(r, limit = 3, morphology = "none")
  expect_equal(d$threshold_mask, matrix(c(FALSE, TRUE, TRUE, FALSE), 2))
})

test_that("NA residuals are never flagged", {
  r <- matrix(c(NA, 5, 0, 0), 2)
  d <- karma2d_detect_anomalies(r, morphology = "none")
  expect_false(d$threshold_mask[1, 1])
})

test_that("asymmetric limits are supported", {
  r <- matrix(c(-2.5, 2.5), 1)
  d <- karma2d_detect_anomalies(r, limit = c(-2, 3), morphology = "none")
  expect_equal(as.vector(d$threshold_mask), c(TRUE, FALSE))
})

test_that("an injected anomaly is detected after opening", {
  y <- simulate_scenario(paper_scenarios$ar, 40, seed = 12)
  fit <- karma2d_fit(y[21:40, 21:40], order = c(1, 0))
  y[5:10, 5:12] <- 0.2
  r <- karma2d_residuals(fit, y, padding = "duplicate")
  det <- karma2d_detect_anomalies(r, limit = 3)
  # most of the anomaly is recovered and nothing far from it is flagged
  inside <- det$detected_mask[5:10, 5:12]
  outside <- det$detected_mask
  outside[3:12, 3:14] <- FALSE
  expect_true(mean(inside) > 0.5 && !any(outside))
})

test_that("anomaly detection prints a summary", {
  d <- karma2d_detect_anomalies(matrix(c(0, 4), 1), morphology = "none")
  expect_output(print(d), "Pixels outside the limits: 1")
})
