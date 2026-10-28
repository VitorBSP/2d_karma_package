# Tolerance for analytic vs numerical derivatives: numDeriv's Richardson
# extrapolation is accurate to about 1e-7 relative here; we use 1e-5
# relative to the gradient norm to stay robust across platforms.

score_check <- function(s, par, size = 15, seed = 1) {
  y <- simulate_scenario(s, size, seed)
  analytic <- karma2d_score(y, par, s$order)
  numeric <- numDeriv::grad(function(p) karma2d_loglik(y, p, s$order), par)
  max(abs(analytic - numeric)) / max(1, max(abs(numeric)))
}

perturbed <- function(s, delta) {
  par <- scenario_parameters(s)
  par + delta * c(rep(0.05, length(par) - 1L), 2)
}

test_that("analytic score matches numerical gradient: AR(1)", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$ar
  expect_lt(score_check(s, scenario_parameters(s)), 1e-5)
})

test_that("analytic score matches numerical gradient: MA(1)", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$ma
  expect_lt(score_check(s, scenario_parameters(s)), 1e-5)
})

test_that("analytic score matches numerical gradient: ARMA(1,1)", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$arma
  expect_lt(score_check(s, scenario_parameters(s)), 1e-5)
})

test_that("analytic score matches numerically away from the truth", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$arma
  expect_lt(max(score_check(s, perturbed(s, 1)),
                score_check(s, perturbed(s, -1), seed = 2)), 1e-5)
})

test_that("analytic score matches numerically for order 2", {
  skip_if_not_installed("numDeriv")
  s <- list(order = c(2L, 1L), alpha = 0.4,
            ar = c(0.2, 0.05, 0.15, 0.1, 0.02, 0.05, 0.02, 0.01),
            ma = c(0.1, 0.1, 0.05), shape = 10)
  expect_lt(score_check(s, scenario_parameters(s), size = 14), 1e-5)
})

test_that("analytic score matches numerically for other links", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$arma
  y <- simulate_scenario(s, 14, 3)
  par <- scenario_parameters(s)
  err <- vapply(c("probit", "cloglog", "cauchit"), function(link) {
    a <- karma2d_score(y, par, s$order, link = link)
    n <- numDeriv::grad(function(p) karma2d_loglik(y, p, s$order, link), par)
    max(abs(a - n)) / max(1, max(abs(n)))
  }, numeric(1))
  expect_lt(max(err), 1e-5)
})

test_that("pixel information K_mu_mu matches numerical integration", {
  cfg <- list(c(0.5, 3), c(0.8, 17), c(0.66, 19.6), c(0.9, 2), c(0.95, 1.2))
  ratio <- vapply(cfg, function(v) {
    package_pixel_information(v[1], v[2])[["mu_mu"]] /
      numeric_pixel_information(v[1], v[2])[["mu_mu"]]
  }, numeric(1))
  expect_equal(ratio, rep(1, length(cfg)), tolerance = 1e-7)
})

test_that("pixel information K_mu_shape matches numerical integration", {
  cfg <- list(c(0.5, 3), c(0.8, 17), c(0.66, 19.6), c(0.9, 2), c(0.95, 1.2))
  ratio <- vapply(cfg, function(v) {
    package_pixel_information(v[1], v[2])[["mu_shape"]] /
      numeric_pixel_information(v[1], v[2])[["mu_shape"]]
  }, numeric(1))
  expect_equal(ratio, rep(1, length(cfg)), tolerance = 1e-7)
})

test_that("pixel information K_shape_shape matches numerical integration", {
  cfg <- list(c(0.5, 3), c(0.8, 17), c(0.66, 19.6), c(0.9, 2), c(0.95, 1.2))
  ratio <- vapply(cfg, function(v) {
    package_pixel_information(v[1], v[2])[["shape_shape"]] /
      numeric_pixel_information(v[1], v[2])[["shape_shape"]]
  }, numeric(1))
  expect_equal(ratio, rep(1, length(cfg)), tolerance = 1e-7)
})

test_that("information is accurate at the removable singularities", {
  # delta = 1 and delta = 2 make individual terms 0/0
  shape <- 2
  mu <- (1 - 0.5^(1 / c(1 + 1e-7, 2 - 1e-7)))^(1 / shape)
  pk <- vapply(mu, function(m) package_pixel_information(m, shape)[3],
               numeric(1))
  nk <- vapply(mu, function(m) numeric_pixel_information(m, shape)[3],
               numeric(1))
  expect_equal(pk, nk, tolerance = 1e-7)
})

test_that("expected information is symmetric", {
  s <- paper_scenarios$arma
  k <- karma2d_information(simulate_scenario(s), scenario_parameters(s),
                           s$order)
  expect_equal(k, t(k), tolerance = 1e-14)
})

test_that("expected information is positive definite at the truth", {
  s <- paper_scenarios$arma
  k <- karma2d_information(simulate_scenario(s, 25), scenario_parameters(s),
                           s$order)
  expect_gt(min(eigen(k, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("observed information equals minus the numerical Hessian", {
  skip_if_not_installed("numDeriv")
  s <- paper_scenarios$ma
  y <- simulate_scenario(s)
  par <- scenario_parameters(s)
  obs <- karma2d_information(y, par, s$order, type = "observed")
  hess <- numDeriv::hessian(function(p) karma2d_loglik(y, p, s$order), par)
  expect_equal(unname(obs), -hess, tolerance = 1e-6)
})

test_that("expected and observed information agree on average", {
  # E[observed information] = expected information; averaging the ratio of
  # traces over replicated images removes most of the sampling noise.
  s <- paper_scenarios$ar
  par <- scenario_parameters(s)
  ratio <- vapply(1:20, function(r) {
    y <- simulate_scenario(s, 20, seed = 100 + r)
    o <- karma2d_information(y, par, s$order, type = "observed")
    e <- karma2d_information(y, par, s$order)
    o[5, 5] / e[5, 5]
  }, numeric(1))
  # sampling sd of the mean ratio is about 0.01 here
  expect_equal(mean(ratio), 1, tolerance = 0.05)
})

test_that("the log-likelihood is -Inf for a non-positive shape", {
  s <- paper_scenarios$ar
  par <- scenario_parameters(s)
  par[5] <- -1
  expect_equal(karma2d_loglik(simulate_scenario(s), par, s$order), -Inf)
})

test_that("a wrong parameter length is an error", {
  s <- paper_scenarios$ar
  expect_error(karma2d_loglik(simulate_scenario(s), c(1, 2, 3), s$order),
               "parameters")
})

test_that("log-likelihood equals the sum of log densities", {
  s <- paper_scenarios$ma
  y <- simulate_scenario(s)
  par <- scenario_parameters(s)
  link <- stats::make.link("logit")
  rec <- karma2d:::karma2d_recursion_r(link$linkfun(y), par[1],
                                       numeric(0), par[2:4], s$order)
  mu <- link$linkinv(rec$eta[2:15, 2:15])
  expect_equal(karma2d_loglik(y, par, s$order),
               sum(dkuma(y[2:15, 2:15], mu, par[5], log = TRUE)),
               tolerance = 1e-12)
})
