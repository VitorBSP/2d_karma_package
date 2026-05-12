dkuma <- function(x, mu, phi, a = 0, b = 1, quantil = 0.5){
  
  
  if(all((x > a) & (x < b))){
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

