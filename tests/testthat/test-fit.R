# Fitted objects are shared across tests to keep the suite fast.
fit_ar <- karma2d_fit(simulate_scenario(paper_scenarios$ar, 30, seed = 1),
                      order = c(1, 0))
fit_ma <- karma2d_fit(simulate_scenario(paper_scenarios$ma, 30, seed = 2),
                      order = c(0, 1))
fit_arma <- karma2d_fit(simulate_scenario(paper_scenarios$arma, 40,
                                          seed = 3), order = c(1, 1))

test_that("AR fit converges", {
  expect_equal(fit_ar$convergence, 0L)
})

test_that("MA fit converges", {
  expect_equal(fit_ma$convergence, 0L)
})

test_that("ARMA fit converges", {
  expect_equal(fit_arma$convergence, 0L)
})

test_that("the score vanishes at the estimate", {
  s <- karma2d_score(fit_arma$data, coef(fit_arma), c(1, 1))
  # relative to the pixel count; reltol of optim leaves O(1e-3) gradients
  expect_lt(max(abs(s)) / nobs(fit_arma), 1e-3)
})

test_that("estimates are close to the truth on a 40 x 40 image", {
  truth <- scenario_parameters(paper_scenarios$arma)
  se <- sqrt(diag(vcov(fit_arma)))
  expect_true(all(abs(coef(fit_arma) - truth) < 4 * se))
})

test_that("coefficients are named in the documented order", {
  expect_equal(names(coef(fit_ma)),
               c("alpha", "theta(0,1)", "theta(1,0)", "theta(1,1)", "shape"))
})

test_that("log and identity shape transforms reach the same maximum", {
  y <- fit_ar$data
  alt <- karma2d_fit(y, c(1, 0),
                     control = karma2d_control(shape_transform = "identity"))
  expect_equal(alt$loglik, fit_ar$loglik, tolerance = 1e-6)
})

test_that("the shape estimate is positive", {
  expect_gt(coef(fit_ma)[["shape"]], 0)
})

test_that("vcov is the inverse of the expected information", {
  expect_equal(vcov(fit_arma) %*% fit_arma$information,
               diag(8), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("observed information can be requested", {
  f <- karma2d_fit(fit_ar$data, c(1, 0),
                   control = karma2d_control(information = "observed"))
  expect_equal(f$inference$type, "observed")
})

test_that("summary reports Wald statistics", {
  tab <- summary(fit_ar)$coefficients
  expect_equal(tab[, "z value"], tab[, "Estimate"] / tab[, "Std. Error"])
})

test_that("confint gives Wald intervals", {
  ci <- confint(fit_ar)
  se <- sqrt(diag(vcov(fit_ar)))
  expect_equal(unname(ci[, 2] - coef(fit_ar)), unname(stats::qnorm(0.975) *
                                                       se))
})

test_that("logLik carries df and nobs", {
  ll <- logLik(fit_arma)
  expect_equal(c(attr(ll, "df"), attr(ll, "nobs")), c(8, 39^2))
})

test_that("AIC is available through logLik", {
  expect_equal(AIC(fit_ar), -2 * fit_ar$loglik + 2 * 5)
})

test_that("print shows the model order", {
  expect_output(print(fit_arma), "2-D KARMA\\(1, 1\\)")
})

test_that("summary prints a coefficient table", {
  expect_output(print(summary(fit_ma)), "theta\\(1,1\\)")
})

test_that("fitted medians are NA on the border only", {
  f <- fitted(fit_ma)
  expect_true(all(is.na(f[1, ])) && !anyNA(f[-1, -1]))
})

test_that("the paper's starting values are used by default", {
  y <- fit_ar$data
  st <- karma2d_start(y, c(1, 0))
  expected <- c(stats::qlogis(mean(y)), 0, 0, 0,
                0.3 * (mean(y) * (1 - mean(y)) / stats::var(as.vector(y)) -
                         1))
  expect_equal(unname(st), expected)
})

test_that("a general support gives the same coefficients as rescaling", {
  y <- fit_ar$data
  f <- karma2d_fit(2 * y - 1, c(1, 0), lower = -1, upper = 1)
  expect_equal(coef(f), coef(fit_ar), tolerance = 1e-8)
})

test_that("data outside the support are rejected", {
  y <- fit_ar$data
  y[1, 1] <- 1
  expect_error(karma2d_fit(y, c(1, 0)), "strictly inside")
})

test_that("missing values are rejected", {
  y <- fit_ar$data
  y[2, 2] <- NA
  expect_error(karma2d_fit(y, c(1, 0)), "missing")
})

test_that("non-convergence is reported with a warning", {
  expect_warning(karma2d_fit(fit_ar$data, c(1, 0),
                             control = karma2d_control(maxit = 2)),
                 "did not converge")
})

test_that("a non positive definite information gives NA covariance", {
  info <- matrix(c(1, 2, 2, 1), 2) # eigenvalues 3 and -1
  expect_warning(out <- karma2d:::covariance_from_information(info),
                 "not positive definite")
  expect_true(all(is.na(out$vcov)))
})

test_that("an ill-conditioned information is flagged", {
  info <- diag(c(1, 1e-14))
  expect_warning(out <- karma2d:::covariance_from_information(info),
                 "ill-conditioned")
  expect_equal(out$status, "ill-conditioned")
})

test_that("no regularization is applied to a valid information", {
  info <- matrix(c(2, 0.5, 0.5, 1), 2)
  out <- karma2d:::covariance_from_information(info)
  expect_equal(out$vcov, solve(info), tolerance = 1e-14)
})
