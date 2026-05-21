escore <- function(z, ynew, ar, ma, n, k, p1, q1, a = 0, b = 1){
  
  p = max(ar) # AR order
  q = max(ma) # MA order
  w <- max(p, q)
  
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
  
  
  link = 'logit'
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog"))){
    stats <- make.link(linktemp)
  } else {
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\"")) 
  }
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv, 
                         mu.eta = stats$mu.eta, 
                         diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  y = ynew
  ynew = linkfun(ynew)
  eta = error = matrix(0, ncol=k,nrow=n)
  dimarray = 1
  
  deta.dalpha = matrix(0, nrow = n, ncol=k)
  deta.dphi = array(0, c(n,k,p1))
  deta.dtheta = array(0, c(n,k,q1))
  
  for (i in (w + 1):n) {
    for (j in (w + 1):k) {
      xx = as.vector(t(ynew[(i - p):i, (j - p):j]))
      yy = as.vector(t(error[(i - q):i, (j - q):j]))
      y_new1 = as.matrix(xx[(length(xx) - 1):1])
      error_new = as.matrix(yy[(length(yy) - 1):1])
      
      eta[i, j] = alpha + phi %*% y_new1 + theta %*% error_new
      error[i, j] = ynew[i, j] - eta[i, j]
      
      # deta.dalpha
      de.da = as.vector(t(deta.dalpha[(i - q):i, (j - q):j]))
      deta.dalpha_anterior = as.matrix(de.da[(length(de.da) - 1):1])
      deta.dalpha[i, j] = 1 - theta %*% deta.dalpha_anterior
      
      # deta.dphi
      for (ii in 0:p) {
        for (jj in 0:p) {
          if (ii == 0 & jj == 0) {
            dimarray = 1
            next
          }
          de.dp = as.vector(t(deta.dphi[(i - q):i, (j - q):j, dimarray]))
          de.dp_anterior = as.matrix(de.dp[(length(de.dp) - 1):1])
          deta.dphi[i, j, dimarray] = ynew[i - ii, j - jj] - theta %*% de.dp_anterior
          dimarray = dimarray + 1
        }
      }
      
      # deta.dtheta
      for (ii in 0:q) {
        for (jj in 0:q) {
          if (ii == 0 & jj == 0) {
            dimarray = 1
            next
          }
          de.dt = as.vector(t(deta.dtheta[(i - q):i, (j - q):j, dimarray]))
          de.dt_anterior = as.matrix(de.dt[(length(de.dt) - 1):1])
          deta.dtheta[i, j, dimarray] = error[i - ii, j - jj] - theta %*% de.dt_anterior
          dimarray = dimarray + 1
        }
      }
    }
  }
  
  
  # VERIFICAR
  
  
  # 1/g'(mu)
  #mu = as.vector(exp(t(eta[(m+1):n,(m+1):k])))
  mu_eta = as.vector(t(mu.eta(eta[(w + 1):n, (w + 1):k])))
  
  # Aqui será o mu.eta
  
  # dell.dmu
  #y1 = as.vector(t(y[(m+1):n,(m+1):k]))
  #dmu = as.vector(((pi * (y1^2)) / (2 * (mu^3)) - (2) / (mu)))
  
  mu <- as.vector(t(linkinv(eta[(w + 1):n, (w + 1):k])))
  y1 <- as.vector(t(y[(w + 1):n, (w + 1):k]))
  
  delta <- log(0.5)/log(1-mu^prec)
  c1 <- mu^(prec-1)/((1-mu^prec)*log(1-mu^prec))
  c2 <- delta*log(1-y1^prec)+1
  vc <- prec*c1*c2
  
  a1 <- 1/prec
  a2 <- log(y1)
  a3 <- vc*mu*log(mu)/prec #PQ TA DIVIDIDO PELO PARAMETRO DE PRECISAO
  a4 <- delta-1
  a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
  prec_estim <- as.vector(a1+a2+a3-(a4*a5)) 
  
  # Aqui será o vc
  
  # deta.dgama 
  #- deta.dalpha
  #- deta.dphi
  #- deta.theta
  
  a = as.vector(deta.dalpha[(w + 1):n, (w + 1):k])
  
  if(p > 0 && q > 0){
    rP <- do.call(
      cbind,
      lapply(1:p1, function(i) as.vector(t(deta.dphi[(w + 1):n, (w + 1):k, i])))
    )
    rT <- do.call(
      cbind,
      lapply(1:q1, function(i) as.vector(t(deta.dtheta[(w + 1):n, (w + 1):k, i])))
    )
    deta.dgamma = cbind(a, rP, rT)
  }else{
    if(p > 0){
      rP = do.call(cbind,  lapply(1:p1,function(i){ as.vector(t(deta.dphi[(w+1):n,(w+1):k,i])) }))
      deta.dgamma = cbind(a, rP)
    }else{
      rT = do.call(cbind,  lapply(1:q1,function(i){ as.vector(t(deta.dtheta[(w+1):n,(w+1):k,i])) }))
      deta.dgamma = cbind(a, rT)
    }
  }

  
  
  
  gamma_par = colSums(vc * mu_eta * deta.dgamma)
  
  # vc * mu.eta * deta.dgamma
  
  return(c(gamma_par, sum(prec_estim)))
}
