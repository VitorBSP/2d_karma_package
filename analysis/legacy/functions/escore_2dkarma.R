escore <- function(z, ynew, AR, MA, n, k, p1, q1, a = 0, b = 1){
  ar = max(AR)
  ma = max(MA)
  if(ar > 0 && ma > 0){
    alpha <- z[1]
    phi = z[2:(p1+1)] 
    theta = z[(p1+2):(p1+q1+1)]
    prec <- z[p1+q1+2] 
  }else{
    if(ar > 0){
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
  
  XX = c()
  
  for (i in (m+1):n)
  {
    for (j in (m+1):k)
    {
      xx1 = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
      XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
    }
  }
  
  P = XX[,2:dim(XX)[2]]
  
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
  y1 = as.vector(t(y[(m+1):n, (m+1):k]))
  
  #############################
  
  # dalpha
  deta_dalpha = matrix(0, nrow = n, ncol=k)
  
  for(i in (m+1):(n-m)){
    for(j in (m+1):(k-m))
    {
      xx = c(deta_dalpha[(i-ma):(i),(j-ma):(j)])
      tam = length(xx)
      #deta_dalpha_anterior = as.matrix(xx[(tam-1):1])
      deta_dalpha_anterior = as.matrix(xx[1:(tam-1)])
      deta_dalpha[i,j] = 1 - theta%*%deta_dalpha_anterior
    }
    
  }
  
  v = as.vector(deta.dalpha[(m+1):(n-m),(m+1):(k-m)])
  
  # dphi 
  P1 = rbind(rep(0,p1),P)
  deta.dphi = matrix(0, ncol=p1,nrow=dim(P1)[1])
  dsum.phi = matrix(0, ncol=p1,nrow=q1)
  
  for(i in m:dim(P1)[1])
  {
    if( i == m)
    {
      deta.dphi[i,]<- P1[i,] 
    }else{
      
      for(j in 1:q1)
      {
        dsum.phi[j,] = theta[j] * deta.dphi[i-ma,]
      }
      deta.dphi[i,]<- P1[i,] - apply(dsum.phi,2,mean)
    }
  }
  
  rP <- deta.dphi[-m,]
  
  #dtheta
  XX = c()
  
  for (i in (m+1):n)
  {
    for (j in (m+1):k)
    {
      xx1 = as.vector(t(error[(i-ma):i,(j-ma):j]))
      XX = rbind(XX,t(as.matrix(xx1[(length(xx1)):1])))
    }
  }
  
  R = rbind(rep(0,q1),XX[,2:dim(XX)[2]])
  deta.dtheta = matrix(0, ncol=q1,nrow=dim(R)[1])
  dsum.theta = matrix(0, ncol=q1,nrow=q1)
  
  for(i in m:dim(R)[1])
  {
    if( i == m)
    {
      deta.dtheta[i,]<- R[i,] 
    }else{
      for(j in 1:q1)
      {
        dsum.theta[j,] = theta[j] * deta.dtheta[i-ma,]
      }
      deta.dtheta[i,]<- R[i,] - apply(dsum.theta,2,mean)
    }
  }
  rR <- deta.dtheta[-m,]
  
  ####
  
  dmu = as.vector(((pi*(y1^2))/(2*(mu^3))-(2)/(mu)))
  
  mT = diag(mu)
  
  
  delta <- log(0.5)/log(1-mu^prec)
  c1 <- mu^(prec-1)/((1-mu^prec)*log(1-mu^prec))
  c2 <- delta*log(1-y1^prec)+1
  vc <- prec*c1*c2
  
  a1 <- 1/prec
  a2 <- log(y1)
  a3 <- vc*mu*log(mu)/prec
  a4 <- delta-1
  a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
  a <- as.vector(a1+a2+a3-(a4*a5))
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  mT <- diag(mu.eta(eta[(m+1):n]))
  
  alpha <- t(v) %*% mT %*% vc
  Uphi <-   t(rP) %*% mT %*% vc
  Utheta <- t(rR) %*% mT %*% vc
  Uprec <- sum(a) 
  
  rval <- c(Ualpha,Uphi,Utheta,Uprec)
}