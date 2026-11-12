# Conditional log-likelihood, score vector and information matrices of the
# 2-D KARMA model. All functions here work on the unit scale (0, 1) and on
# the original parametrization gamma = (alpha, phi, theta, shape).
#
# Notation (per pixel of the effective region): u is mu^shape, delta is
# log(0.5) / log(1 - u), A is shape mu^(shape - 1) / ((1 - u) log(1 - u)),
# so that the derivative of the log-density with respect to mu is
# A (1 + delta log(1 - y^shape)), i.e. shape times c[n, m] of Appendix A,
# and B is the derivative of log(delta) with respect to the shape,
# u log(mu) / ((1 - u) log(1 - u)).
# Under the model, E := -delta * log(1 - Y^shape) ~ Exp(1) and
# X := Y^shape ~ Beta(1, delta); the expected information below follows from
# these two facts. See docs/technical-audit.md, Section "Fisher information".

euler_gamma <- -digamma(1)

# Evaluate f(delta) avoiding the removable singularity at delta = at by
# averaging the values at at +/- h (error O(h^2), f is smooth).
avoid_removable <- function(f, delta, at, h = 1e-4) {
  out <- f(delta)
  near <- !is.na(delta) & abs(delta - at) < h
  if (any(near)) {
    out[near] <- 0.5 * (f(delta[near] - h) + f(delta[near] + h))
  }
  out
}

# h(delta) = E[X log X / (1 - X)], X ~ Beta(1, delta)
#          = (psi(2) - psi(delta + 1)) / (delta - 1)
kuma_h <- function(delta) {
  avoid_removable(function(d) (digamma(2) - digamma(d + 1)) / (d - 1),
                  delta, at = 1)
}

# t(delta) = (delta - 1) / delta * E[X log(X)^2 / (1 - X)^2], X ~ Beta(1, delta)
#          = [(psi(delta) - psi(2))^2 + psi'(2) - psi'(delta)] / (delta - 2)
kuma_t <- function(delta) {
  avoid_removable(function(d) {
    ((digamma(d) - digamma(2))^2 + trigamma(2) - trigamma(d)) / (d - 2)
  }, delta, at = 2)
}

# Pixel-level quantities shared by the likelihood, score and information.
# Returns NULL when the parameters give an invalid model (non-positive shape,
# medians numerically equal to 0 or 1, or non-finite values).
karma2d_quantities <- function(par, y, y_link, order, link_obj,
                               derivatives = FALSE) {
  parts <- split_parameters(par, order)
  if (!is.finite(parts$shape) || parts$shape <= 0 ||
        any(!is.finite(c(parts$alpha, parts$ar, parts$ma)))) {
    return(NULL)
  }
  rec <- karma2d_recursion(y_link, parts$alpha, parts$ar, parts$ma, order,
                           derivatives = derivatives)
  w <- max(order)
  eff <- which(effective_region(nrow(y), ncol(y), w))
  eta <- rec$eta[eff]
  mu <- link_obj$linkinv(eta)
  shape <- parts$shape
  if (any(!is.finite(eta)) || any(!is.finite(mu)) ||
        any(mu <= 0 | mu >= 1)) {
    return(NULL)
  }
  u <- mu^shape
  log1mu <- log1p(-u)
  delta <- log(0.5) / log1mu
  if (any(!is.finite(delta)) || any(delta <= 0)) {
    return(NULL)
  }
  y_eff <- y[eff]
  log_y <- log(y_eff)
  shape_log_y <- shape * log_y
  list(
    parts = parts,
    eff = eff,
    rec = rec,
    eta = eta,
    mu = mu,
    mu_eta = link_obj$mu.eta(eta),
    y = y_eff,
    shape = shape,
    u = u,
    delta = delta,
    log_y = log_y,
    # log of one minus y^shape, and y^shape over one minus y^shape
    log1my = log1mexp(shape_log_y),
    yphi_ratio = exp(shape_log_y) / -expm1(shape_log_y),
    ratio_u = ratio_u_log1mu(u)
  )
}

karma2d_loglik_internal <- function(par, y, y_link, order, link_obj) {
  qt <- karma2d_quantities(par, y, y_link, order, link_obj)
  if (is.null(qt)) {
    return(-Inf)
  }
  ll <- log(qt$shape) + log(qt$delta) + (qt$shape - 1) * qt$log_y +
    (qt$delta - 1) * qt$log1my
  out <- sum(ll)
  if (is.na(out)) -Inf else out
}

score_from_quantities <- function(qt) {
  shape <- qt$shape
  mu <- qt$mu
  a_term <- shape * qt$ratio_u / (mu * (1 - qt$u))
  b_term <- qt$ratio_u * log(mu) / (1 - qt$u)
  centred <- 1 + qt$delta * qt$log1my # 1 - E, mean zero under the model
  score_mu <- a_term * centred
  deriv <- qt$rec$derivatives[, qt$eff, drop = FALSE]
  score_gamma <- drop(deriv %*% (score_mu * qt$mu_eta))
  score_shape <- sum(1 / shape + qt$log_y + b_term * centred -
                       (qt$delta - 1) * qt$yphi_ratio * qt$log_y)
  c(score_gamma, score_shape)
}

karma2d_score_internal <- function(par, y, y_link, order, link_obj) {
  qt <- karma2d_quantities(par, y, y_link, order, link_obj,
                           derivatives = TRUE)
  if (is.null(qt)) {
    return(rep(NA_real_, length(par)))
  }
  score_from_quantities(qt)
}

# Expected (conditional) Fisher information, summed over the effective region.
karma2d_expected_information <- function(par, y, y_link, order, link_obj) {
  qt <- karma2d_quantities(par, y, y_link, order, link_obj,
                           derivatives = TRUE)
  n_par <- length(par)
  if (is.null(qt)) {
    return(matrix(NA_real_, n_par, n_par))
  }
  shape <- qt$shape
  mu <- qt$mu
  delta <- qt$delta
  a_term <- shape * qt$ratio_u / (mu * (1 - qt$u))
  b_term <- qt$ratio_u * log(mu) / (1 - qt$u)
  h <- kuma_h(delta)
  t_term <- kuma_t(delta)

  k_mu_mu <- a_term^2
  k_mu_shape <- a_term * (b_term + delta * h / shape)
  k_shape_shape <- 1 / shape^2 + b_term^2 + 2 * delta * b_term * h / shape +
    delta * t_term / shape^2

  deriv <- qt$rec$derivatives[, qt$eff, drop = FALSE]
  k_gamma_gamma <- deriv %*% (t(deriv) * (k_mu_mu * qt$mu_eta^2))
  k_gamma_shape <- drop(deriv %*% (k_mu_shape * qt$mu_eta))
  info <- rbind(cbind(k_gamma_gamma, k_gamma_shape),
                c(k_gamma_shape, sum(k_shape_shape)))
  dimnames(info) <- NULL
  0.5 * (info + t(info))
}

# Observed information: minus the Hessian of the log-likelihood, obtained by
# central differences of the analytic score.
karma2d_observed_information <- function(par, y, y_link, order, link_obj,
                                         rel_step = 1e-5) {
  n_par <- length(par)
  hess <- matrix(NA_real_, n_par, n_par)
  for (i in seq_len(n_par)) {
    h <- rel_step * max(abs(par[i]), 1)
    if (i == n_par) h <- min(h, par[i] / 4) # keep shape positive
    up <- down <- par
    up[i] <- par[i] + h
    down[i] <- par[i] - h
    hess[, i] <- (karma2d_score_internal(up, y, y_link, order, link_obj) -
                    karma2d_score_internal(down, y, y_link, order, link_obj)) /
      (2 * h)
  }
  -0.5 * (hess + t(hess))
}
