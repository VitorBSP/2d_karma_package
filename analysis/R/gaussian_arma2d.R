# Gaussian 2-D ARMA benchmark of the article (Section 4.2).
#
# y[n,m] = alpha + sum phi(i,j) y[n-i,m-j] + sum theta(k,l) e[n-k,m-l] + e[n,m],
# e ~ N(0, sigma^2), with the same strongly causal neighbourhood, lag
# ordering and conditioning on the first w rows/columns as 2-D KARMA, fitted
# to the rescaled (0, 1) data on the identity scale. Fitted by conditional
# maximum likelihood with BFGS over (alpha, phi, theta, log sigma) using the
# analytic gradient. This replaces analysis/legacy/functions/escore_arma.R,
# whose gradient used a different neighbourhood from its log-likelihood.
#
# Uses the internal recursion of karma2d (karma2d:::karma2d_recursion_r).

gaussian_arma2d_fit <- function(y, order) {
  p1 <- (order[1] + 1)^2 - 1
  q1 <- (order[2] + 1)^2 - 1
  w <- max(order)
  eff <- row(y) > w & col(y) > w
  n_par <- p1 + q1 + 2
  unpack <- function(par) list(alpha = par[1], ar = par[seq_len(p1) + 1],
                               ma = par[seq_len(q1) + 1 + p1],
                               sigma = exp(par[n_par]))
  run <- function(par, derivatives = FALSE) {
    u <- unpack(par)
    karma2d:::karma2d_recursion_r(y, u$alpha, u$ar, u$ma, order,
                                  derivatives = derivatives)
  }
  negloglik <- function(par) {
    e <- run(par)$error[eff]
    -sum(stats::dnorm(e, sd = exp(par[n_par]), log = TRUE))
  }
  gradient <- function(par) {
    rec <- run(par, derivatives = TRUE)
    e <- rec$error[eff]
    s2 <- exp(2 * par[n_par])
    d <- rec$derivatives[, which(eff), drop = FALSE]
    # d e / d gamma = - d eta / d gamma
    -c(drop(d %*% e) / s2, sum(e^2 / s2 - 1))
  }
  start <- c(mean(y), rep(0, p1 + q1), log(stats::sd(as.vector(y))))
  opt <- stats::optim(start, negloglik, gradient, method = "BFGS",
                      control = list(maxit = 500), hessian = TRUE)
  est <- opt$par
  # delta method for sigma = exp(log sigma)
  vc <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  jac <- diag(c(rep(1, n_par - 1), exp(est[n_par])))
  v <- if (is.null(vc)) rep(NA_real_, n_par) else diag(jac %*% vc %*% jac)
  se <- ifelse(is.finite(v) & v > 0, sqrt(pmax(v, 0)), NA_real_)
  coef <- c(est[-n_par], exp(est[n_par]))
  names(coef) <- c("alpha", karma2d_lags(order[1], "phi")$name,
                   karma2d_lags(order[2], "theta")$name, "sigma")
  list(coefficients = coef, se = stats::setNames(se, names(coef)),
       loglik = -opt$value, convergence = opt$convergence, order = order)
}

gaussian_arma2d_predict <- function(fit, newdata, padding = "duplicate") {
  order <- fit$order
  w <- max(order)
  cf <- fit$coefficients
  p1 <- (order[1] + 1)^2 - 1
  q1 <- (order[2] + 1)^2 - 1
  z <- if (padding == "duplicate" && w > 0)
    karma2d:::pad_duplicate(newdata, w) else newdata
  rec <- karma2d:::karma2d_recursion_r(z, cf[1], cf[seq_len(p1) + 1],
                                       cf[seq_len(q1) + 1 + p1], order)
  if (padding == "duplicate" && w > 0) {
    keep_r <- (2 * w + 1):nrow(z)
    keep_c <- (2 * w + 1):ncol(z)
    return(rec$eta[keep_r, keep_c])
  }
  eta <- rec$eta
  eta[!(row(eta) > w & col(eta) > w)] <- NA
  eta
}
