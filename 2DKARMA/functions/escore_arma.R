escore_gaussiana <- function(z, ynew, AR, MA, n, k, p1, q1, a = 0, b = 1){
  ar = max(AR)
  ma = max(MA)
  
  if(ar > 0 && ma > 0){
    alpha <- z[1]
    phi = z[2:(p1+1)] 
    theta = z[(p1+2):(p1+q1+1)]
    sd <- z[p1+q1+2] 
  } else if(ar > 0){
    alpha <- z[1]
    phi = z[2:(p1+1)] 
    theta = 0
    sd <- z[p1+2] 
  } else {
    alpha <- z[1]
    phi = 0
    theta = z[2:(q1+1)]
    sd <- z[q1+2] 
  }
  
  m = max(ar, ma)
  eta = error = matrix(0, ncol = k, nrow = n)
  
  # Inicializar gradientes com zero
  num_params = length(z)
  grad = rep(0, num_params)
  
  for (i in (m+1):n){
    for (j in (m+1):k){
      if(ar > 0){
        phi_ = phi
        xx = as.vector(t(ynew[(i-ar):(i-1),(j-ar):(j-1)]))
        y_new1 = as.matrix(xx[length(xx):1])
      } else {
        phi_ = 0
        y_new1 = 0
      }
      
      if(ma > 0){
        theta_ = theta
        yy = as.vector(t(error[(i-ma):(i-1),(j-ma):(j-1)]))
        error_new = as.matrix(yy[length(yy):1])
      } else {
        theta_ = 0
        error_new = 0
      }
      
      eta[i,j]  = alpha + phi_%*%y_new1 + theta_%*%error_new
      error[i,j] = ynew[i,j] - eta[i,j]
      
      # Cálculo dos gradientes
      idx_alpha = 1
      d_eta_alpha = 1
      grad[idx_alpha] = grad[idx_alpha] - (1 / sd^2) * error[i,j] * d_eta_alpha
      
      # Parâmetros AR (phi)
      if(ar > 0){
        idx_phi_start = 2
        for(k_phi in 1:length(phi)){
          d_eta_phi = y_new1[k_phi]
          grad[idx_phi_start + k_phi - 1] = grad[idx_phi_start + k_phi - 1] - (1 / sd^2) * error[i,j] * d_eta_phi
        }
      }
      
      # Parâmetros MA (theta)
      if(ma > 0){
        idx_theta_start = if(ar > 0) (p1 + 2) else 2
        for(k_theta in 1:length(theta)){
          d_eta_theta = error_new[k_theta]
          grad[idx_theta_start + k_theta - 1] = grad[idx_theta_start + k_theta - 1] - (1 / sd^2) * error[i,j] * d_eta_theta
        }
      }
    }
  }
  
  # Gradiente em relação a sd
  idx_sd = length(z)
  N = (n - m) * (k - m)
  sum_error2 = sum(error[(m+1):n, (m+1):k]^2)
  grad[idx_sd] = -N * (1 / sd) + (1 / sd^3) * sum_error2
  
  return(grad)
}