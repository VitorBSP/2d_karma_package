fisher_information_2dkarma_vectorized <- function(
    z,
    ynew,
    ar,
    ma,
    n,
    k,
    p1,
    q1,
    log_prec_scale = FALSE
) {
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
  prec <- quantities$prec
  delta <- quantities$delta
  mu_eta <- quantities$mu_eta
  
  normalize_block <- function(block, n_rows) {
    if (is.null(block)) {
      return(matrix(0, nrow = n_rows, ncol = 0))
    }
    
    if (is.vector(block)) {
      return(matrix(block, ncol = 1))
    }
    
    block
  }
  
  n_obs <- length(quantities$mu)
  
  alpha_derivative_vector <- quantities$alpha_derivative_vector
  phi_derivative_matrix <- normalize_block(
    quantities$phi_derivative_matrix,
    n_rows = n_obs
  )
  theta_derivative_matrix <- normalize_block(
    quantities$theta_derivative_matrix,
    n_rows = n_obs
  )
  
  y1 <- quantities$y1
  
  euler_constant <- 0.5772156649
  
  lambda_1 <- mu^(prec - 2) /
    ((1 - mu^prec) * log(1 - mu^prec))
  
  lambda_2 <- mu^(2 * prec - 2) /
    ((1 - mu^prec)^2 * (log(1 - mu^prec))^2)
  
  weight_gamma_gamma <- -(prec^2) * lambda_2 * (mu_eta^2)
  
  weight_gamma_prec <- -prec * mu * (
    delta * lambda_1 *
      ((1 - digamma(delta + 1) - euler_constant) / ((delta - 1) * prec)) +
      log(mu) * lambda_2
  ) * mu_eta
  
  c_term <- quantities$vc / prec
  
  l1 <- 1 / (prec^2)
  l2 <- delta * (mu^2) * lambda_2 * (log(mu)^2) * log(1 - y1^prec)
  l3 <- 2 * delta * (mu^2) * log(mu) * lambda_1 *
    ((y1^prec * log(y1)) / (1 - y1^prec))
  l4 <- (delta - 1) * (y1^prec * (log(y1)^2)) / ((1 - y1^prec)^2)
  l5 <- c_term * mu * (log(mu)^2) *
    (lambda_1 * (mu^2) + 1 / (1 - mu^prec))
  
  weight_prec_prec <- l1 - l2 + l3 + l4 - l5
  
  weighted_alpha <- weight_gamma_gamma * alpha_derivative_vector
  
  Kaa <- matrix(
    -sum(alpha_derivative_vector * weighted_alpha),
    nrow = 1,
    ncol = 1
  )
  
  if (ncol(phi_derivative_matrix) > 0) {
    weighted_phi <- weight_gamma_gamma * phi_derivative_matrix
    Kap <- matrix(
      -colSums(alpha_derivative_vector * weighted_phi),
      nrow = 1
    )
    Kpp <- -crossprod(phi_derivative_matrix, weighted_phi)
  } else {
    Kap <- matrix(0, nrow = 1, ncol = 0)
    Kpp <- matrix(0, nrow = 0, ncol = 0)
  }
  
  if (ncol(theta_derivative_matrix) > 0) {
    weighted_theta <- weight_gamma_gamma * theta_derivative_matrix
    Kat <- matrix(
      -colSums(alpha_derivative_vector * weighted_theta),
      nrow = 1
    )
    Ktt <- -crossprod(theta_derivative_matrix, weighted_theta)
  } else {
    Kat <- matrix(0, nrow = 1, ncol = 0)
    Ktt <- matrix(0, nrow = 0, ncol = 0)
  }
  
  if (ncol(phi_derivative_matrix) > 0 && ncol(theta_derivative_matrix) > 0) {
    weighted_theta <- weight_gamma_gamma * theta_derivative_matrix
    Kpt <- -crossprod(phi_derivative_matrix, weighted_theta)
  } else {
    Kpt <- matrix(0, nrow = ncol(phi_derivative_matrix), ncol = ncol(theta_derivative_matrix))
  }
  
  Kaprec <- matrix(
    -sum(alpha_derivative_vector * weight_gamma_prec),
    nrow = 1,
    ncol = 1
  )
  
  if (ncol(phi_derivative_matrix) > 0) {
    Kpprec <- matrix(
      -colSums(phi_derivative_matrix * weight_gamma_prec),
      ncol = 1
    )
  } else {
    Kpprec <- matrix(0, nrow = 0, ncol = 1)
  }
  
  if (ncol(theta_derivative_matrix) > 0) {
    Ktprec <- matrix(
      -colSums(theta_derivative_matrix * weight_gamma_prec),
      ncol = 1
    )
  } else {
    Ktprec <- matrix(0, nrow = 0, ncol = 1)
  }
  
  Kprecprec <- matrix(
    sum(weight_prec_prec),
    nrow = 1,
    ncol = 1
  )
  
  upper_block <- cbind(Kaa, Kap, Kat, Kaprec)
  middle_phi_block <- cbind(t(Kap), Kpp, Kpt, Kpprec)
  middle_theta_block <- cbind(t(Kat), t(Kpt), Ktt, Ktprec)
  lower_block <- cbind(t(Kaprec), t(Kpprec), t(Ktprec), Kprecprec)
  
  fisher_matrix <- rbind(
    upper_block,
    middle_phi_block,
    middle_theta_block,
    lower_block
  )
  
  fisher_matrix <- 0.5 * (fisher_matrix + t(fisher_matrix))
  
  if (log_prec_scale) {
    jacobian_matrix <- diag(
      c(rep(1, nrow(fisher_matrix) - 1), prec)
    )
    fisher_matrix <- jacobian_matrix %*% fisher_matrix %*% jacobian_matrix
  }
  
  fisher_matrix
}