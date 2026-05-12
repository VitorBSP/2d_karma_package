con <- file("resultado_2dkarma10.txt", open = "wt")
sink(con)


library(numDeriv) 
library(tidyverse)
library(moments)
library(fastrep)

#Valores:

replicas = 10000 #replicas
n_ = c(10, 25, 50) # tamanho das imagens

# Parametros
alpha =  0.67 # intercept (beta)
phi = c(0.31, 0.21, -0.15)
theta = c() # MA parameters
ar = c(1)
ma = c(0)
prec = 23


# Funções básicas

dkuma <- function(x, mu, phi, a = 0, b = 1, quantil = 0.5){
  
  
  if(all((x > 0) & (x < 1))){
    # range = (1/(b-a))
    # not_kernel = (phi*log(quantil))
    # log_not_kernel = (log(1-mu^phi))
    # kernel_expressao1 = x^(phi-1)
    # kernel_expressao2 = (1-x^phi)^((log(quantil)/log_not_kernel) - 1)
    # d = range*(not_kernel/log_not_kernel)*kernel_expressao1*kernel_expressao2
    d<-(1/(b-a))*((phi*log(0.5))/(log(1-mu^phi)))*x^(phi-1)*(1-x^phi)^( ((log(0.5))/(log(1-mu^phi)))-1 )
    
    return(d)
    
    # not_kernel = phi ^ delta / (b - a)
    # kernel = (x - a) / (b - a)
    #
    # not_kernel * kernel ^ (phi - 1) * (1 - kernel ^ phi) ^ (delta - 1)
    
  }else{
    stop("x must be in unit range.")
  }
  
}

pkuma <- function(x, mu, phi, a = 0, b = 1, quantil = 0.5){
  
  
  
  if(all((x > 0) & (x < 1))){
    kernel = (1 - x^phi)
    potencia = (log(quantil)/(log(1-mu^phi)))
    p = 1 - kernel^(potencia)
    return(p)
    # kernel = (x - a) / (b - a)
    # 1 - (1 - kernel ^ phi) ^ delta
    
  }else{
    stop("x must be in unit range.")
  }
  
  
}

qkuma <- function(p, mu, phi, a = 0, b = 1, quantil = 0.5){
  
  # kernel = (1-p)
  # potencia = (log(1-mu^phi))/(log(quantil))
  #
  # q = a + (b-a)*(1 - kernel^potencia)^(1/phi)
  
  q<- a + (b-a)*( 1-(1-p)^( (log(1-mu^phi))/(log(0.5)) ) )^(1/phi)
  return(q)
  
  
  # kernel = (1 - p) ^ (1/mu)
  #
  # a + (b - a)*(1 - kernel) ^ (1/phi)
  
}


rkuma <- function(mu, phi, a = 0, b = 1, quantil = 0.5){
  
  n = length(mu)
  
  u = runif(n)
  qkuma(u, mu, phi, a, b, quantil)
}

# Simulação

simulacao = function (alpha = NA, phi = NA, theta = NA, n = NA, k = NA, AR = NA, 
                      MA = NA, link = 'logit', prec = 10, ne = 1){
  
  # n is the number of rows
  # k is the number of columns 
  # ne is related to the employed neighborhood 
  # ne = 1: strongly causal neighborhood (used in Palm, Bayer, and Cintra ?)
  # ne = 2: non causal region
  # ne = 3: semi-causal
  
  # Apenas realizando com ne = 1
  
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
    }
    
  }
  
  y = y[(2*m+1):(n+2*m), (2*m+1):(k+2*m)]
  return(y)
}



# Verossimilhança

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
  y1 = as.vector(t(y[(m+1):n, (m+1):k]))
  ll = suppressWarnings(log(dkuma(y1,mu,prec, a=a, b=b)))
  return(sum(ll))
} 


# Escore

escore <- function(z, ynew, ar, ma, n, k, p1, q1, a = 0, b = 1){
  
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
  
  eta = error = matrix(0, ncol=k,nrow=n)
  deta.dalpha = matrix(0, nrow = n, ncol=k)
  deta.dphi = deta.dtheta = array(0, c(n,k,p1))
  dimarray = 1
  
  deta.dalpha = matrix(0, nrow = n, ncol=k)
  deta.dphi = array(0, c(n,k,p1))
  deta.dtheta = array(0, c(n,k,q1))
  
  for (i in (m + 1):n) {
    for (j in (m + 1):k) {
      xx = as.vector(t(ynew[(i - ar):i, (j - ar):j]))
      y_new1 = as.matrix(xx[(length(xx) - 1):1])
      
      yy = as.vector(t(error[(i - ma):i, (j - ma):j]))
      error_new = as.matrix(yy[(length(yy) - 1):1])
      
      eta[i, j] = alpha + phi %*% y_new1 + theta %*% error_new
      error[i, j] = ynew[i, j] - eta[i, j]
      
      # deta.dalpha
      de.da = as.vector(deta.dalpha[(i - ma):i, (j - ma):j])
      deta.dalpha_anterior = as.matrix(de.da[(length(de.da) - 1):1])
      deta.dalpha[i, j] = 1 - theta %*% deta.dalpha_anterior
      
      # deta.dphi
      for (ii in 0:ar) {
        for (jj in 0:ar) {
          if (ii == 0 & jj == 0) {
            dimarray = 1
            next
          }
          de.dp = as.vector(t(deta.dphi[(i - ma):i, (j - ma):j, dimarray]))
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
          de.dt = as.vector(t(deta.dtheta[(i - ma):i, (j - ma):j, dimarray]))
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
  mu_eta = as.vector(t(mu.eta(eta[(m+1):n, (m+1):k])))
  
  # Aqui será o mu.eta
  
  # dell.dmu
  #y1 = as.vector(t(y[(m+1):n,(m+1):k]))
  #dmu = as.vector(((pi * (y1^2)) / (2 * (mu^3)) - (2) / (mu)))
  
  y1 = as.vector(t(y[(m+1):n,(m+1):k]))
  mu = as.vector(linkinv(t(eta[(m+1):n,(m+1):k])))
  
  delta <- log(0.5)/log(1-mu^prec)
  c1 <- mu^(prec-1)/((1-mu^prec)*log(1-mu^prec))
  c2 <- delta*log(1-y1^prec)+1
  vc <- prec*c1*c2
  
  a1 <- 1/prec
  a2 <- log(y1)
  a3 <- vc*mu*log(mu)/prec
  a4 <- delta-1
  a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
  prec_estim <- as.vector(a1+a2+a3-(a4*a5)) 
  
  # Aqui será o vc
  
  # deta.dgama 
  #- deta.dalpha
  #- deta.dphi
  #- deta.theta
  
  a = as.vector(deta.dalpha[(m+1):n,(m+1):k])
  
  if(p > 0 && q > 0){
    rP = do.call(cbind,  lapply(1:p1,function(i){ as.vector(t(deta.dphi[(m+1):n,(m+1):k,i])) }))
    rT = do.call(cbind,  lapply(1:q1,function(i){ as.vector(t(deta.dtheta[(m+1):n,(m+1):k,i])) }))
    deta.dgamma = cbind(a, rP, rT)
  }else{
    if(p > 0){
      rP = do.call(cbind,  lapply(1:p1,function(i){ as.vector(t(deta.dphi[(m+1):n,(m+1):k,i])) }))
      deta.dgamma = cbind(a, rP)
    }else{
      rT = do.call(cbind,  lapply(1:q1,function(i){ as.vector(t(deta.dtheta[(m+1):n,(m+1):k,i])) }))
      deta.dgamma = cbind(a, rT)
    }
  }
  
  
  
  
  gamma_par = colSums(vc * mu_eta * deta.dgamma)
  
  # vc * mu.eta * deta.dgamma
  
  return(c(gamma_par, sum(prec_estim)))
}


# Monte carlo

parameter_estimation <- function(reg){
  opt = optim(reg, loglik, 
              escore,
              method = "BFGS", 
              hessian = TRUE,
              control = list(fnscale = -1), ynew = ynew, ar = ar, ma = ma, n = n,
              k = k, p1 = p1, q1=q1)
  
  print(opt$par)
  return(opt$par)
}


comparator_deriv <- function(reg){
  opt = optim(reg, loglik, 
              method = "BFGS", 
              hessian = FALSE,
              control = list(fnscale = -1), ynew = ynew, ar = ar, ma = ma, n = n,
              k = k, p1 = p1, q1=q1)
  
  print(opt$par)
  
  escore_analitic_numeric = cbind(analitic = escore(opt$par, ynew = ynew, ar = ar, ma = ma, n = n,
                                                    k = k, p1 = p1, q1=q1),
                                  numeric = grad(loglik, opt$par, ynew = ynew, ar = ar, ma = ma, n = n,
                                                 k = k, p1 = p1, q1=q1))
  
  print(escore_analitic_numeric)
  # print(opt_escore$par)
  return(opt$par)
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



estimate_params_10 = matrix(nrow=replicas,ncol=length(c(alpha, phi, theta, prec)))
estimate_params_25 = matrix(nrow=replicas,ncol=length(c(alpha, phi, theta, prec)))
estimate_params_50 = matrix(nrow=replicas,ncol=length(c(alpha, phi, theta, prec)))

teste1 = NULL

for(n in n_){
  for(i in 1:replicas){
    y = simulacao(alpha = alpha, phi = phi, theta = theta, 
                  n = n, k = n, AR = ar, MA = ma, prec = prec)
    p = max(ar) # AR order
    q = max(ma) # MA order
    n = dim(y)[1] # n is the number of rows
    k = dim(y)[2] # k is the number of columns 
    m = max(p,q,na.rm=T)
    p1 = (p+1)^2-1
    q1 = (q+1)^2-1
    ynew = linkfun(y)
    var_y = sd(y)^2
    mean_y = mean(y)
    chute_prec = ((mean_y*(1-mean_y))/var_y - 1)*0.3
    teste1 = append(teste1, chute_prec) 
    reg = c(linkfun(mean(y)), rep(0, p1), rep(0, q1), chute_prec)
    true_param = c(alpha, phi, theta, prec)
    estim = try(parameter_estimation(reg))
    if((n == n_[2]) & (i == 1)){
      try(comparator_deriv(reg))
    }
    if(class(estim) != 'try-error'){
      if(n == 10){
        print("foi")
        estimate_params_10[i,] = estim
        
      }else{
        if(n == 25){
          estimate_params_25[i,] = estim
        }else{
          estimate_params_50[i,] = estim
        }
      }
    }
  }
  
  
}

if(n == n_[1]){
  print(estimate_params_10 %>% data.frame() %>%
          colMeans(na.rm=T) %>%
          cbind(apply(estimate_params_10 %>% data.frame(), 2, sd, na.rm=T)) %>%
          as.data.frame() %>% 
          rename("Média" = ".", "Desvio Padrão" = "V2") %>%
          mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
          mutate(Viés = Média - Parâmetros) %>%
          mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
                 Assimetria = moments::skewness(estimate_params_10 %>% data.frame(), 
                                                na.rm = T),
                 Curtose = moments::kurtosis(estimate_params_10 %>% data.frame(),
                                             na.rm=T)) %>%
          mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
          select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`, 
                 EQM, Assimetria, Curtose)  %>% 
          fastrep::tbl(format='latex'))
  
}else{
  if(n == n_[2]){
    print(estimate_params_25 %>% data.frame() %>%
            colMeans(na.rm=T) %>%
            cbind(apply(estimate_params_25 %>% data.frame(), 2, sd, na.rm=T)) %>%
            as.data.frame() %>% 
            rename("Média" = ".", "Desvio Padrão" = "V2") %>%
            mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
            mutate(Viés = Média - Parâmetros) %>%
            mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
                   Assimetria = moments::skewness(estimate_params_25 %>% 
                                                    data.frame(), na.rm = T),
                   Curtose = moments::kurtosis(estimate_params_25 %>% 
                                                 data.frame(), na.rm=T)) %>%
            mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
            select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`,
                   EQM, Assimetria, Curtose) %>% 
            fastrep::tbl(format='latex'))
  }else{
    print(estimate_params_50 %>% data.frame() %>%
            colMeans(na.rm=T) %>%
            cbind(apply(estimate_params_50 %>% data.frame(), 2, sd, na.rm=T)) %>%
            as.data.frame() %>% 
            rename("Média" = ".", "Desvio Padrão" = "V2") %>%
            mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
            mutate(Viés = Média - Parâmetros) %>%
            mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
                   Assimetria = moments::skewness(estimate_params_50 %>% data.frame(), 
                                                  na.rm = T),
                   Curtose = moments::kurtosis(estimate_params_50 %>% data.frame(),
                                               na.rm=T)) %>%
            mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
            select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`, EQM, 
                   Assimetria, Curtose) %>%
            fastrep::tbl(format='latex'))
  }
}



c(alpha, phi, theta, prec)


estimate_params_50 %>% data.frame() %>%
  colMeans(na.rm=T) %>%
  cbind(apply(estimate_params_50 %>% data.frame(), 2, sd, na.rm=T)) %>%
  as.data.frame() %>% 
  rename("Média" = ".", "Desvio Padrão" = "V2") %>%
  mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
  mutate(Viés = Média - Parâmetros) %>%
  mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
         Assimetria = moments::skewness(estimate_params_50 %>% data.frame(), 
                                        na.rm = T),
         Curtose = moments::kurtosis(estimate_params_50 %>% data.frame(),
                                     na.rm=T)) %>%
  mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
  select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`, EQM, 
         Assimetria, Curtose)

estimate_params_25 %>% data.frame() %>%
  colMeans(na.rm=T) %>%
  cbind(apply(estimate_params_25 %>% data.frame(), 2, sd, na.rm=T)) %>%
  as.data.frame() %>% 
  rename("Média" = ".", "Desvio Padrão" = "V2") %>%
  mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
  mutate(Viés = Média - Parâmetros) %>%
  mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
         Assimetria = moments::skewness(estimate_params_25 %>% data.frame(), 
                                        na.rm = T),
         Curtose = moments::kurtosis(estimate_params_25 %>% data.frame(),
                                     na.rm=T)) %>%
  mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
  select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`, EQM, 
         Assimetria, Curtose)

estimate_params_10 %>% data.frame() %>%
  colMeans(na.rm=T) %>%
  cbind(apply(estimate_params_10 %>% data.frame(), 2, sd, na.rm=T)) %>%
  as.data.frame() %>% 
  rename("Média" = ".", "Desvio Padrão" = "V2") %>%
  mutate(Parâmetros = c(alpha, phi, theta, prec)) %>%
  mutate(Viés = Média - Parâmetros) %>%
  mutate(EQM =  `Viés`^2 + `Desvio Padrão`^2,
         Assimetria = moments::skewness(estimate_params_10 %>% data.frame(), 
                                        na.rm = T),
         Curtose = moments::kurtosis(estimate_params_10 %>% data.frame(),
                                     na.rm=T)) %>%
  mutate(`Viés Relativo` =  Viés/Parâmetros*100) %>%
  select(Parâmetros, Média, `Desvio Padrão`, `Viés`, `Viés Relativo`, EQM, 
         Assimetria, Curtose)

sink(type = "message")
close(con)
