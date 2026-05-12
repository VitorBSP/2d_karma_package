loglik <- function(z, ynew, ar, ma, n, k, p1, q1, a = 0, b = 1){
  
  p = max(ar) # AR order
  q = max(ma) # MA order
  
  if(p > 0 && q > 0){
    alpha <- z[1]
    phi = z[2:(p1+1)] 
    theta = z[(p1+2):(p1+q1+1)]
    prec <- z[p1+q1+2] 
  }else{
    if(p > 0){
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = 0
      prec <- z[p1+2] 
    }else{
      alpha <- z[1]
      phi = 0
      theta = z[2:(q1+1)]
      prec <- z[q1+2] 
    }
  }
  
  eta = error = matrix(0, ncol=k,nrow=n)
  
  for (i in (m+1):(n)){
    for (j in (m+1):(k)){
      phi_ = phi
      theta_ = theta
      if(ar > 0){
        xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
      }else{
        phi_ = 0
        y_new1 = 0
      }
      
      if(ma > 0){
        yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
        error_new = as.matrix(yy[(length(yy)-1):1])
      }else{
        theta_ = 0
        error_new = 0
      }
      
      eta[i,j]  = alpha + phi_%*%y_new1 + theta_%*%error_new
      error[i,j] = ynew[i,j]-eta[i,j]
      
    }
    
  }
  
  mu = as.vector(t(linkinv(eta[(m+1):n, (m+1):k])))
  y1 = as.vector(t(ynew[(m+1):n, (m+1):k]))
  ll = suppressWarnings(log(dkuma(y1,mu,prec, a=a, b=b)))
  return(sum(ll))
} 
