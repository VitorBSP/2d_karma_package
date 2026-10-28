test_that("simulation returns an image of the requested size", {
  y <- simulate_scenario(paper_scenarios$arma, size = 12)
  expect_equal(dim(y), c(12L, 12L))
})

test_that("non-square images are supported", {
  y <- karma2d_simulate(8, 13, alpha = 0.5, ar = c(0.2, 0.1, 0.1),
                        shape = 10, seed = 1)
  expect_equal(dim(y), c(8L, 13L))
})

test_that("simulated values lie strictly inside the unit interval", {
  y <- simulate_scenario(paper_scenarios$ma, size = 30)
  expect_true(all(y > 0 & y < 1))
})

test_that("simulated values respect a general support", {
  y <- karma2d_simulate(10, 10, alpha = 0.5, ar = c(0.2, 0.1, 0.1),
                        shape = 10, lower = -1, upper = 1, seed = 2)
  expect_true(all(y > -1 & y < 1))
})

test_that("the same seed gives the same image", {
  s <- paper_scenarios$arma
  expect_identical(simulate_scenario(s, seed = 5),
                   simulate_scenario(s, seed = 5))
})

test_that("different seeds give different images", {
  s <- paper_scenarios$arma
  expect_false(identical(simulate_scenario(s, seed = 5),
                         simulate_scenario(s, seed = 6)))
})

test_that("a seed does not change the global random stream", {
  set.seed(99)
  before <- .Random.seed
  simulate_scenario(paper_scenarios$ar, seed = 1)
  expect_identical(.Random.seed, before)
})

test_that("without a seed the global stream is used", {
  s <- paper_scenarios$ar
  set.seed(7)
  a <- simulate_scenario(s, seed = NULL)
  set.seed(7)
  b <- simulate_scenario(s, seed = NULL)
  expect_identical(a, b)
})

test_that("a larger burn-in still returns the requested size", {
  y <- karma2d_simulate(10, 10, alpha = 0.3, ma = c(0.1, 0.1, 0.1),
                        shape = 8, burn_in = 30, seed = 3)
  expect_equal(dim(y), c(10L, 10L))
})

test_that("burn-in smaller than the order is rejected", {
  expect_error(karma2d_simulate(10, 10, alpha = 0.3,
                                ma = rep(0.05, 8), shape = 8, burn_in = 1),
               "burn_in")
})

test_that("order (0, 0) gives iid pixels with the requested median", {
  y <- karma2d_simulate(80, 80, alpha = stats::qlogis(0.7), shape = 6,
                        seed = 4)
  expect_equal(stats::median(y), 0.7, tolerance = 0.01)
})

test_that("an explosive process stops with an informative error", {
  expect_error(karma2d_simulate(30, 30, alpha = 2, ar = c(0.9, 0.9, 0.9),
                                shape = 5, seed = 1), "unstable")
})

test_that("the generating model is stored as an attribute", {
  y <- simulate_scenario(paper_scenarios$ar)
  expect_equal(attr(y, "karma2d")$order, c(p = 1L, q = 0L))
})

test_that("simulation is independent of the working directory", {
  s <- paper_scenarios$ar
  a <- simulate_scenario(s, seed = 8)
  old <- setwd(tempdir())
  on.exit(setwd(old))
  expect_identical(simulate_scenario(s, seed = 8), a)
})
