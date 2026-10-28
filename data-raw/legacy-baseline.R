# Reference values produced by the ORIGINAL (legacy) implementation.
#
# This script sources the legacy code in analysis/legacy/ without modifying
# it and records its numerical output. The result is stored in
# tests/testthat/fixtures/legacy-baseline.rds and is used by the test suite
# to check that the refactored implementation reproduces the legacy behaviour
# wherever the legacy behaviour was audited as correct (see
# docs/technical-audit.md).
#
# Run from the repository root:
#   Rscript data-raw/legacy-baseline.R
#
# Requires numDeriv. The legacy code is only sourced here; it is never part of
# the installed package.

legacy_dir <- file.path("analysis", "legacy")

# Legacy files call source() with paths relative to an old layout. Evaluate
# them with those source() calls removed; the needed files are sourced
# explicitly below.
source_without_source_calls <- function(path, envir) {
  code <- readLines(path, warn = FALSE, encoding = "UTF-8")
  code <- code[!grepl("^\\s*source\\(", code)]
  eval(parse(text = code), envir = envir)
}

legacy_env <- function(files) {
  env <- new.env(parent = globalenv())
  for (f in files) source_without_source_calls(file.path(legacy_dir, f), env)
  env
}

# Kumaraswamy functions -----------------------------------------------------
dist_env <- legacy_env("src/functions.R")
x_grid <- c(0.01, 0.1, 0.35, 0.5, 0.8, 0.99)
mu_grid <- c(0.2, 0.5, 0.75)
shape_grid <- c(0.7, 3, 19.6)
grid <- expand.grid(x = x_grid, mu = mu_grid, shape = shape_grid)
kuma <- list(
  grid = grid,
  density = with(grid, dist_env$dkuma(x, mu, shape)),
  cdf = with(grid, dist_env$pkuma(x, mu, shape)),
  quantile = with(grid, dist_env$qkuma(x, mu, shape))
)

# Model quantities ----------------------------------------------------------
model_files <- c(
  "src/functions.R",
  "functions/2dkarma_base.R",
  "functions/log_vero_2dkarma.R",
  "functions/new_escore_2dkarma.R",
  "functions/fisher_2dkarma.R",
  "functions/fisher_2dkarma_vectorized.R",
  "functions/simu_2dkarma.R"
)
env <- legacy_env(model_files)
reformed_env <- legacy_env(c("src/functions.R",
                             "functions/escore_reformed_2dkarma.R"))

scenarios <- list(
  ar1 = list(order = c(1, 0), alpha = 0.6609,
             ar = c(0.2819, 0.1897, 0.1653), ma = numeric(0),
             shape = 19.6065),
  ma1 = list(order = c(0, 1), alpha = 1.1114, ar = numeric(0),
             ma = c(0.3467, 0.2819, 0.1848), shape = 12.3857),
  arma11 = list(order = c(1, 1), alpha = 0.3,
                ar = c(0.3, -0.15, 0.23), ma = c(0.15, -0.1, 0.25),
                shape = 17)
)

image_size <- 12L
cases <- list()
for (name in names(scenarios)) {
  s <- scenarios[[name]]
  set.seed(20260101)
  invisible(utils::capture.output(
    y <- env$simulacao(
      alpha = s$alpha, phi = s$ar, theta = s$ma,
      n = image_size, k = image_size,
      AR = s$order[1], MA = s$order[2], prec = s$shape
    )
  ))
  p1 <- (s$order[1] + 1)^2 - 1
  q1 <- (s$order[2] + 1)^2 - 1
  z <- c(s$alpha, s$ar, s$ma, s$shape)
  args <- list(z = z, ynew = y, ar = s$order[1], ma = s$order[2],
               n = image_size, k = image_size, p1 = p1, q1 = q1)
  loglik_fun <- function(par) {
    a <- args
    a$z <- par
    do.call(env$loglik, a)
  }
  cases[[name]] <- list(
    scenario = s,
    seed = 20260101,
    y = y,
    parameters = z,
    loglik = do.call(env$loglik, args),
    score_base = do.call(env$escore_reformed, args),
    score_reformed = do.call(reformed_env$escore, args),
    score_numeric = numDeriv::grad(loglik_fun, z),
    hessian_numeric = numDeriv::hessian(loglik_fun, z),
    fisher_legacy = do.call(env$fisher_information_2dkarma, args),
    fisher_legacy_vectorized =
      do.call(env$fisher_information_2dkarma_vectorized, args)
  )
}

baseline <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  numderiv_version = as.character(utils::packageVersion("numDeriv")),
  kumaraswamy = kuma,
  cases = cases
)

out <- file.path("tests", "testthat", "fixtures", "legacy-baseline.rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(baseline, out, version = 2)
message("Wrote ", out)
