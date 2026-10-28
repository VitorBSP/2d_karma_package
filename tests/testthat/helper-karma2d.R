# Parameter sets of Tables 1-3 of the article.
paper_scenarios <- list(
  ar = list(order = c(1L, 0L), alpha = 0.6609,
            ar = c(0.2819, 0.1897, 0.1653), ma = numeric(0),
            shape = 19.6065),
  ma = list(order = c(0L, 1L), alpha = 1.1114, ar = numeric(0),
            ma = c(0.3467, 0.2819, 0.1848), shape = 12.3857),
  arma = list(order = c(1L, 1L), alpha = 0.3, ar = c(0.3, -0.15, 0.23),
              ma = c(0.15, -0.1, 0.25), shape = 17)
)

scenario_parameters <- function(s) c(s$alpha, s$ar, s$ma, s$shape)

simulate_scenario <- function(s, size = 15, seed = 1) {
  karma2d_simulate(size, size, alpha = s$alpha, ar = s$ar, ma = s$ma,
                   shape = s$shape, seed = seed)
}

# Pixel-level expected information of (mu, shape) by numerical integration.
# Written in terms of E = -delta log(1 - Y^shape) ~ Exp(1), which avoids
# evaluating the density near y = 1; independent of the package code.
numeric_pixel_information <- function(mu, shape) {
  log1mu <- log1p(-mu^shape)
  delta <- log(0.5) / log1mu
  scores <- function(t) {
    r <- exp(-t / delta)
    x <- -expm1(-t / delta) # y^shape
    log_x <- log1p(-r)
    log_y <- log_x / shape
    ratio <- ifelse(r > 1e-12, log_x / r, -1 - r / 2) # log(x) / (1 - x)
    c_term <- mu^(shape - 1) / ((1 - mu^shape) * log1mu) * (1 - t)
    cbind(shape * c_term,
          1 / shape + log_y + c_term * mu * log(mu) -
            (delta - 1) * x * ratio / shape)
  }
  entry <- function(a, b) {
    stats::integrate(function(t) {
      s <- scores(t)
      s[, a] * s[, b] * exp(-t)
    }, 0, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value
  }
  c(mu_mu = entry(1, 1), mu_shape = entry(1, 2), shape_shape = entry(2, 2))
}

# Pixel-level expected information from the package formulas.
package_pixel_information <- function(mu, shape) {
  y <- matrix(0.5, 2, 2) # one modelled pixel for a (0, 0) model
  link <- stats::make.link("logit")
  par <- c(link$linkfun(mu), shape)
  k <- karma2d:::karma2d_expected_information(
    par, y, link$linkfun(y), c(0L, 0L), link
  )
  # K is in (alpha, shape); alpha = eta and d mu / d eta = mu (1 - mu)
  d <- link$mu.eta(par[1])
  c(mu_mu = k[1, 1] / d^2 / 4, mu_shape = k[1, 2] / d / 4,
    shape_shape = k[2, 2] / 4)
}
