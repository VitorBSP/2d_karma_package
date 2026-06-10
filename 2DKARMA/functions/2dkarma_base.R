build_2dkarma_quantities <- function(z, ynew, ar, ma, n, k, p1, q1, link = 'logit') {
  p <- max(ar)
  q <- max(ma)
  w <- max(p, q)
  
  if (p > 0 && q > 0) {
    alpha <- z[1]
    phi <- z[2:(p1 + 1)]
    theta <- z[(p1 + 2):(p1 + q1 + 1)]
    prec <- z[p1 + q1 + 2]
  } else if (p > 0) {
    alpha <- z[1]
    phi <- z[2:(p1 + 1)]
    theta <- numeric(0)
    prec <- z[p1 + 2]
  } else {
    alpha <- z[1]
    phi <- numeric(0)
    theta <- z[2:(q1 + 1)]
    prec <- z[q1 + 2]
  }
  
  link_name <- link
  link_stats <- make.link(link_name)
  linkfun <- link_stats$linkfun
  linkinv <- link_stats$linkinv
  mu_eta_function <- link_stats$mu.eta
  
  y <- ynew
  y_link <- linkfun(ynew)
  
  eta <- matrix(0, nrow = n, ncol = k)
  error <- matrix(0, nrow = n, ncol = k)
  
  alpha_derivative <- matrix(0, nrow = n, ncol = k)
  phi_derivative <- array(0, c(n, k, p1))
  theta_derivative <- array(0, c(n, k, q1))
  
  for (i in (w + 1):n) {
    for (j in (w + 1):k) {
      if (p > 0) {
        phi_window <- as.vector(t(y_link[(i - p):i, (j - p):j]))
        phi_lag <- matrix(phi_window[(length(phi_window) - 1):1], ncol = 1)
        ar_component <- drop(phi %*% phi_lag)
      } else {
        ar_component <- 0
      }
      
      if (q > 0) {
        theta_window <- as.vector(t(error[(i - q):i, (j - q):j]))
        theta_lag <- matrix(theta_window[(length(theta_window) - 1):1], ncol = 1)
        ma_component <- drop(theta %*% theta_lag)
      } else {
        ma_component <- 0
      }
      
      eta[i, j] <- alpha + ar_component + ma_component
      error[i, j] <- y_link[i, j] - eta[i, j]
      
      if (q > 0) {
        alpha_window <- as.vector(t(alpha_derivative[(i - q):i, (j - q):j]))
        alpha_lag <- matrix(alpha_window[(length(alpha_window) - 1):1], ncol = 1)
        alpha_derivative[i, j] <- 1 - drop(theta %*% alpha_lag)
      } else {
        alpha_derivative[i, j] <- 1
      }
      
      if (p > 0) {
        phi_index <- 1
        
        for (ii in 0:p) {
          for (jj in 0:p) {
            if (ii == 0 && jj == 0) {
              next
            }
            
            if (q > 0) {
              phi_derivative_window <- as.vector(
                t(phi_derivative[(i - q):i, (j - q):j, phi_index])
              )
              phi_derivative_lag <- matrix(
                phi_derivative_window[(length(phi_derivative_window) - 1):1],
                ncol = 1
              )
              correction_term <- drop(theta %*% phi_derivative_lag)
            } else {
              correction_term <- 0
            }
            
            phi_derivative[i, j, phi_index] <- y_link[i - ii, j - jj] - correction_term
            phi_index <- phi_index + 1
          }
        }
      }
      
      if (q > 0) {
        theta_index <- 1
        
        for (ii in 0:q) {
          for (jj in 0:q) {
            if (ii == 0 && jj == 0) {
              next
            }
            
            theta_derivative_window <- as.vector(
              t(theta_derivative[(i - q):i, (j - q):j, theta_index])
            )
            theta_derivative_lag <- matrix(
              theta_derivative_window[(length(theta_derivative_window) - 1):1],
              ncol = 1
            )
            
            theta_derivative[i, j, theta_index] <- error[i - ii, j - jj] -
              drop(theta %*% theta_derivative_lag)
            
            theta_index <- theta_index + 1
          }
        }
      }
    }
  }
  
  eta_effective <- eta[(w + 1):n, (w + 1):k, drop = FALSE]
  mu <- as.vector(t(linkinv(eta_effective)))
  y1 <- as.vector(t(y[(w + 1):n, (w + 1):k, drop = FALSE]))
  mu_eta <- as.vector(t(mu_eta_function(eta_effective)))
  
  delta <- log(0.5) / log(1 - mu^prec)
  c1 <- mu^(prec - 1) / ((1 - mu^prec) * log(1 - mu^prec))
  c2 <- delta * log(1 - y1^prec) + 1
  vc <- prec * c1 * c2
  
  alpha_derivative_vector <- as.vector(
    t(alpha_derivative[(w + 1):n, (w + 1):k])
  )
  
  if (p > 0) {
    phi_derivative_matrix <- do.call(
      cbind,
      lapply(
        seq_len(p1),
        function(index) {
          as.vector(t(phi_derivative[(w + 1):n, (w + 1):k, index]))
        }
      )
    )
  } else {
    phi_derivative_matrix <- matrix(
      0,
      nrow = length(alpha_derivative_vector),
      ncol = 0
    )
  }
  
  if (q > 0) {
    theta_derivative_matrix <- do.call(
      cbind,
      lapply(
        seq_len(q1),
        function(index) {
          as.vector(t(theta_derivative[(w + 1):n, (w + 1):k, index]))
        }
      )
    )
  } else {
    theta_derivative_matrix <- matrix(
      0,
      nrow = length(alpha_derivative_vector),
      ncol = 0
    )
  }
  
  derivative_matrix <- cbind(
    matrix(alpha_derivative_vector, ncol = 1),
    phi_derivative_matrix,
    theta_derivative_matrix
  )
  
  precision_score_term <- {
    a1 <- 1 / prec
    a2 <- log(y1)
    a3 <- vc * mu * log(mu) / prec
    a4 <- delta - 1
    a5 <- (y1^prec * log(y1)) / (1 - y1^prec)
    as.vector(a1 + a2 + a3 - (a4 * a5))
  }
  
  list(
    p = p,
    q = q,
    w = w,
    alpha = alpha,
    phi = phi,
    theta = theta,
    prec = prec,
    eta = eta,
    error = error,
    mu = mu,
    y1 = y1,
    mu_eta = mu_eta,
    delta = delta,
    vc = vc,
    alpha_derivative_vector = alpha_derivative_vector,
    phi_derivative_matrix = phi_derivative_matrix,
    theta_derivative_matrix = theta_derivative_matrix,
    derivative_matrix = derivative_matrix,
    precision_score_term = precision_score_term
  )
}