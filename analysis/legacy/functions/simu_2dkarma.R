simulacao = function (alpha = NA, phi = NA, theta = NA, n = NA, k = NA, AR = NA, 
                      MA = NA, link = 'logit', prec = 10, ne = 1){
  
  # n is the number of rows
  # k is the number of columns 
  # ne is related to the employed neighborhood 
  # ne = 1: strongly causal neighborhood (used in Palm, Bayer, and Cintra ?)
  # ne = 2: non causal region
  # ne = 3: semi-causal
  
  # Apenas realizando com ne = 1
  
  source("./2DKARMA/src/functions.R")
  #ar = p
  #ma = q
  
  ar = max(AR)
  ma = max(MA)
  p1 = (ar+1)^2-1
  q1 = (ma+1)^2-1
 
  if(p1 != length(phi)){
    num_phi = sum(AR*2 + 1)
    if(num_phi != length(phi)){
      stop("For an autoregressive model ", ar, ' order needs ', num_phi, 
           ' parameters for phi.')  
    }
    
  } 
  
  if(q1 != length(theta)){
    num_theta = sum(MA*2 + 1)
    if(num_theta != length(theta)){
      stop("For an moving average model ", ma, ' order needs ', num_theta, 
           ' parameters for theta.')  
    }
    
  } 
  
  
  
  
  # q <- max(ma)
  # n <- length(y)
  # 
  # m <- max(p,q,na.rm=T)
  # 
  # p1 <- length(ar)
  # q1 <- length(ma)
  
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
  
  # if(any(is.na(phi)==F))
  # {
  #   ar = 1:p
  # }
  # 
  # if(any(is.na(theta)==F))
  # {
  #   ma = 1:q
  # }
  
  ###### ARMA model
  
    m = max(ar,ma)
    
    ynew = matrix(1,ncol=(k+2*m),nrow=(n+2*m)) * alpha 
    
    ynew1 = c()
    
    y = matrix(0, ncol=(k+2*m),nrow=(n+2*m))
    eta = error = mu = y 
    
    # p1 = max(ar)
    # p2 = max(ma)
    for (i in (m+1):(n+2*m))
    {
      for (j in (m+1):(k+2*m))
      {
        theta_ = theta
        phi_ = phi
        if(length(phi) == 0){
          phi_ = 0
          y_new1 = 0
        }else{
          
          if(length(phi) == p1){ 
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
          } else {
            
            xx = as.vector(t(ynew[(i-ar):i,(j-ar):j]))
            y_new1 = as.matrix(xx[(length(xx)-1):1])
            
            full_ar = 1:ar
            ar_insignificant = full_ar[!full_ar %in% AR]
            ar_index = matrix(TRUE, nrow(y_new1))
            
            for(not_ar in ar_insignificant){
              first_discarted_index = ((not_ar-1)+1)^2
              last_discarted_index = ((not_ar)+1)^2 -1
              ar_index[(first_discarted_index): last_discarted_index] = FALSE
            }
            
            y_new1 = as.matrix(y_new1[ar_index])
          }
        }
        
        
        if(length(theta) == 0){
          theta_ = 0
          error_new = 0
        }else{
          if(length(theta) == q1){
            
            yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
          } else{
            yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
            error_new = as.matrix(yy[(length(yy)-1):1])
            
            full_ma = 1:ma
            ma_insignificant = full_ma[!full_ma %in% MA]
            ma_index = matrix(TRUE, nrow(error_new))
            
            
            for(not_ma in ma_insignificant){
              first_discarted_index = ((not_ma-1)+1)^2
              last_discarted_index = ((not_ma)+1)^2 - 1
              ma_index[(first_discarted_index): last_discarted_index] = FALSE
            }
            
            error_new = as.matrix(error_new[ma_index])
          }
        }
        
        
        eta[i,j]  = alpha + phi_%*%y_new1 + theta_%*%error_new
        
        mu[i,j] = linkinv(eta[i,j])
        y[i,j] = rkuma(mu[i,j], prec)
        ynew[i,j] = linkfun(y[i,j])
        
        error[i,j] = ynew[i,j]-linkfun(mu[i,j])
        
        print(i)
        print(j)
        print(y_new1)
      }
      
    }
  
  y = y[(2*m+1):(n+2*m), (2*m+1):(k+2*m)]
  return(y)
}


# simulacao(alpha = 0.5, phi = c(0.3, 0.5, 0.23), theta = c(0.3, 0.7, 0.25), n = 10, k = 10, p = c(1), q = 1)