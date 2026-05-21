loglik_gaussiana <- function(z, ynew, AR, MA, n, k, p1, q1, a = 0, b = 1) {
  ar <- max(AR)
  ma <- max(MA)
  w <- max(ar, ma)
  
  if (ar > 0 && ma > 0) {
    alpha <- z[1]
    phi <- z[2:(p1 + 1)]
    theta <- z[(p1 + 2):(p1 + q1 + 1)]
    log_sd <- z[p1 + q1 + 2]
  } else if (ar > 0) {
    alpha <- z[1]
    phi <- z[2:(p1 + 1)]
    theta <- numeric(0)
    log_sd <- z[p1 + 2]
  } else {
    alpha <- z[1]
    phi <- numeric(0)
    theta <- z[2:(q1 + 1)]
    log_sd <- z[q1 + 2]
  }
  
  sd_value <- exp(log_sd)
  
  eta <- matrix(0, nrow = n, ncol = k)
  error <- matrix(0, nrow = n, ncol = k)
  
  for (i in (w + 1):n) {
    for (j in (w + 1):k) {
      ar_component <- 0
      ma_component <- 0
      
      if (ar > 0) {
        ar_window <- as.vector(t(ynew[(i - ar):i, (j - ar):j, drop = FALSE]))
        y_lag <- matrix(ar_window[(length(ar_window) - 1):1], ncol = 1)
        ar_component <- drop(phi %*% y_lag)
      }
      
      if (ma > 0) {
        ma_window <- as.vector(t(error[(i - ma):i, (j - ma):j, drop = FALSE]))
        error_lag <- matrix(ma_window[(length(ma_window) - 1):1], ncol = 1)
        ma_component <- drop(theta %*% error_lag)
      }
      
      eta[i, j] <- alpha + ar_component + ma_component
      error[i, j] <- ynew[i, j] - eta[i, j]
    }
  }
  
  mu <- as.vector(t(eta[(w + 1):n, (w + 1):k, drop = FALSE]))
  y_obs <- as.vector(t(ynew[(w + 1):n, (w + 1):k, drop = FALSE]))
  
  sum(dnorm(y_obs, mean = mu, sd = sd_value, log = TRUE))
}