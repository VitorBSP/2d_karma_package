# The fixture was produced by data-raw/legacy-baseline.R, which runs the
# original (unmodified) scripts in analysis/legacy/. Tolerances reflect
# floating point reordering only, except where noted.

baseline <- readRDS(test_path("fixtures", "legacy-baseline.rds"))

for (case_name in names(baseline$cases)) {
  cs <- baseline$cases[[case_name]]
  s <- cs$scenario
  np <- length(cs$parameters)

  test_that(paste("simulation reproduces the legacy image:", case_name), {
    y <- karma2d_simulate(nrow(cs$y), ncol(cs$y), alpha = s$alpha,
                          ar = s$ar, ma = s$ma, shape = s$shape,
                          seed = cs$seed)
    expect_equal(unclass(y)[, ], cs$y, tolerance = 1e-10,
                 ignore_attr = TRUE)
  })

  test_that(paste("log-likelihood reproduces legacy:", case_name), {
    expect_equal(karma2d_loglik(cs$y, cs$parameters, s$order), cs$loglik,
                 tolerance = 1e-12)
  })

  test_that(paste("score reproduces legacy:", case_name), {
    expect_equal(unname(karma2d_score(cs$y, cs$parameters, s$order)),
                 cs$score_base, tolerance = 1e-10)
  })

  label <- paste("gamma blocks of the information reproduce legacy:",
                 case_name)
  test_that(label, {
    k <- karma2d_information(cs$y, cs$parameters, s$order)
    expect_equal(unname(k[, -np]), cs$fisher_legacy[, -np],
                 tolerance = 1e-10)
  })

  label <- paste("legacy shape-shape entry was the observed one:", case_name)
  test_that(label, {
    # Documents the audit finding: the legacy "Fisher" (shape, shape) entry
    # equals the observed information, not the expected one.
    o <- karma2d_information(cs$y, cs$parameters, s$order,
                             type = "observed")
    expect_equal(o[np, np], cs$fisher_legacy[np, np], tolerance = 1e-6)
  })
}

test_that("density reproduces legacy where legacy is accurate", {
  g <- baseline$kumaraswamy$grid
  accurate <- g$mu^g$shape > 1e-8 # legacy used log(1 - u), not log1p(-u)
  expect_equal(with(g[accurate, ], dkuma(x, mu, shape)),
               baseline$kumaraswamy$density[accurate], tolerance = 1e-10)
})

test_that("quantiles reproduce legacy where legacy is accurate", {
  g <- baseline$kumaraswamy$grid
  accurate <- g$mu^g$shape > 1e-8
  expect_equal(with(g[accurate, ], qkuma(x, mu, shape)),
               baseline$kumaraswamy$quantile[accurate], tolerance = 1e-10)
})
