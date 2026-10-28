source("./2DKARMA/functions/2dkarma_base.R")

fisher_information_2dkarma <- function(z, ynew, ar, ma, n, k, p1, q1) {
  quantities <- build_2dkarma_quantities(
    z = z,
    ynew = ynew,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  mu <- quantities$mu
  y1 <- quantities$y1
  mu_eta <- quantities$mu_eta
  delta <- quantities$delta
  prec <- quantities$prec
  derivative_matrix <- quantities$derivative_matrix
  
  lambda1 <- mu^(prec - 2) /
    ((1 - mu^prec) * log(1 - mu^prec))
  
  lambda2 <- mu^(2 * prec - 2) /
    ((1 - mu^prec)^2 * (log(1 - mu^prec))^2)
  
  euler_constant <- 0.5772156649
  
  information_weight <- -prec^2 * lambda2 * mu_eta^2
  
  d_term <- -prec * mu * log(mu) * lambda2 -
    prec * mu * delta * lambda1 *
    ((1 - digamma(delta + 1) - euler_constant) / ((delta - 1) * prec))
  
  cross_weight <- d_term * mu_eta
  
  c_term <- quantities$vc / prec
  
  l1 <- 1 / (prec^2)
  l2 <- delta * (mu^2) * lambda2 * (log(mu)^2) * log(1 - y1^prec)
  l3 <- 2 * delta * (mu^2) * log(mu) * lambda1 *
    ((y1^prec * log(y1)) / (1 - y1^prec))
  l4 <- (delta - 1) * (y1^prec * (log(y1)^2)) / ((1 - y1^prec)^2)
  l5 <- c_term * mu * (log(mu)^2) *
    (lambda1 * (mu^2) + 1 / (1 - mu^prec))
  
  precision_information <- -(l1 - l2 + l3 + l4 - l5)
  
  gamma_gamma_block <- -crossprod(
    derivative_matrix,
    information_weight * derivative_matrix
  )
  
  gamma_precision_block <- -matrix(
    colSums(derivative_matrix * cross_weight),
    ncol = 1
  )
  
  precision_precision_block <- -matrix(sum(precision_information), nrow = 1)
  
  rbind(
    cbind(gamma_gamma_block, gamma_precision_block),
    cbind(t(gamma_precision_block), precision_precision_block)
  )
}