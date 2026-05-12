############################# IMPORTS ############################# 
library(abind)
library(raster)
library(rgdal)
library(EBImage)
library(tidyverse)
library(reshape2)
source("./2DKARMA/functions/log_vero_arma.R")

############################# FUNÇÕES AUXILIARES ############################# 

preencher_com_media <- function(matriz) {
  media <- mean(matriz, na.rm = TRUE)
  matriz_preenchida <- matriz
  
  matriz_preenchida[is.na(matriz)] <- media
  
  return(matriz_preenchida)
}

parameter_estimation_arma <- function(reg, ynew, AR, MA, n, k, p1, q1, a = 0, b = 1){
  opt = optim(reg, loglik_gaussiana, 
              method = "BFGS", 
              hessian = FALSE,
              control = list(fnscale = -1), 
              ynew = ynew, AR = AR, MA = MA, n = n, k = k, 
              p1 = p1, q1 = q1, a = a, b = b)
  return(opt$par)
}


mse <- function(verdadeiro, predito){
  metrica = mean((verdadeiro - predito)^2)
  return(metrica)
}

mape <- function(verdadeiro, predito){
  metrica = mean(abs(verdadeiro - predito / verdadeiro)) * 100
  return(metrica)
}

mae <- function(verdadeiro, predito){
  metrica = sum(abs(verdadeiro - predito))/length(as.vector(predito))
  return(metrica)
}


r_2 <- function(verdadeiro, predito){
  metrica = cor(as.vector(verdadeiro),as.vector(predito))^2
  return(metrica)
}

rmsle <- function(verdadeiro, predito){
  metrica = sqrt(mean((log1p(verdadeiro) - log1p(predito))^2))
  return(metrica)
}

arma_mle_inference <- function(reg, ynew, AR, MA, n, k, p1, q1, a = 0, b = 1, 
                               chol = FALSE, simetrizacao = FALSE ) {
  opt <- optim(
    par     = reg,
    fn      = loglik_gaussiana,
    method  = "BFGS",
    hessian = TRUE,                      # <- pede Hessiana numérica no ponto ótimo
    control = list(fnscale = -1),        # <- maximiza loglik (minimizando -loglik)
    ynew = ynew, AR = AR, MA = MA, n = n, k = k, p1 = p1, q1 = q1, a = a, b = b
  )
  
  theta_hat <- opt$par
  if (is.null(names(theta_hat))) names(theta_hat) <- paste0("par", seq_along(theta_hat))
  
  H <- opt$hessian
  # Informação observada: -ℓ''(θ̂)
  
  if (simetrizacao == FALSE){
    HESSIAN = -H
  }
  else {
    HESSIAN <- -0.5 * (H + t(H)) 
  }
  
  if ((chol == FALSE)){
    cov_theta <- tryCatch(
      solve(HESSIAN),
      error = function(e) {
        warning("Hessian mal-condicionada.")
        matrix(NA_real_, length(theta_hat), length(theta_hat))
      }
    )
  } else {
    cov_theta <- tryCatch(
      chol2inv(chol(HESSIAN)),     # Var(θ̂) ≈ A^{-1}
      error = function(e) {
        lam <- 1e-8 * max(abs(diag(A)), na.rm = TRUE)
        A_r <- A + diag(lam, nrow(A))
        tryCatch(
          chol2inv(chol(A_r)),
          error = function(e2) {
            warning("Hessian ainda mal-condicionada.")
            if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(A_r) else matrix(NA_real_, nrow(A), ncol(A))
          }
        )
      }
    )
  } 
  
  
  
  se  <- sqrt(diag(cov_theta))
  z   <- theta_hat / se
  p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  infer_table <- data.frame(
    estimate = theta_hat,
    std_error = se,
    z_value = z,
    p_value = p,
    row.names = names(theta_hat),
    check.names = FALSE
  )
  
  return(list(
    coef        = theta_hat,
    vcov        = cov_theta,
    table       = infer_table,
    convergence = opt$convergence,
    message     = opt$message,
    counts      = opt$counts,
    value       = opt$value  # valor do objetivo no ótimo (aqui, -loglik)
  ))
}

############################# GRAFICOS ############################# 

def_graph <- function(df, ar, ma){
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = "Times", size = 14, 
                                     margin = margin(0, 0, 0, 0.2, unit = "cm")),
          axis.text.y = element_text(family = "Times", size = 14, 
                                     margin = margin(0.2, 0, 0, 0, unit = "cm"))) +
    scale_fill_gradient2(low = "blue", mid = "lightgreen", high = "red", midpoint = 0) +
    labs(fill = "") + 
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "cm"))
}

def_threshold <- function(df, ar, ma){
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = "Times", size = 14, 
                                     margin = margin(0, 0, 0, 0.2, unit = "cm")),
          axis.text.y = element_text(family = "Times", size = 14, 
                                     margin = margin(0.2, 0, 0, 0, unit = "cm"))) +
    scale_fill_gradient2(low = "white", mid = "black", high = "white", midpoint = 0) +
    labs(fill = "") + 
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "cm"))
}

def_graph_ndvi <- function(df, ar, ma) {
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = "Times", size = 14, 
                                     margin = margin(0, 0, 0, 0.2, unit = "cm")),
          axis.text.y = element_text(family = "Times", size = 14, 
                                     margin = margin(0.2, 0, 0, 0, unit = "cm"))) +
    scale_fill_gradient2(low = "blue", mid = "red", high = "darkgreen", midpoint = 0) +
    labs(fill = "") + 
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "cm"))
}


def_hist <- function(dados, binwidth = NULL){
  ggplot(dados, aes(x = value)) + 
    geom_histogram(aes(y = after_stat(count)),
                   color = "black", fill = "transparent",
                   binwidth = binwidth) +  # defina para controlar a largura das barras
    scale_x_continuous(breaks = c(-10, -5, -3 , 0, 3, 5, 10)) +
    ylab("Count") +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(family = "Times", size = 12),
      axis.text.y  = element_text(family = "Times", size = 12),
      axis.title.y = element_text(family = "Times", size = 16)
    )
}



############################# LINK ############################# 

link = 'identity'
linktemp <- substitute(link)
if (!is.character(linktemp)){
  linktemp <- deparse(linktemp)
  if (linktemp == "link")
    linktemp <- eval(link)
}
if (any(linktemp == c("logit", "probit", "cloglog", 'identity'))){
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

############################# DADOS ############################# 

All<- list.files("./2DKARMA/ndvi/ndvi", full.names = TRUE, 
                 # full.names = TRUE, 
                 pattern = "*.tif$"
                 # pattern = glob2rx("NDVI_*.tif$")
)
All
stack_all<-stack(All)
elem = as.matrix(stack_all[[32]])

y_true_selected = (elem - (-1))/(1 - (-1))

rownames(y_true_selected) = 1:dim(y_true_selected)[1]
colnames(y_true_selected) = 1:dim(y_true_selected)[2]

y_train = y_true_selected[51:100, 51:100]
y_anomaly = y_true_selected
y_anomaly[5:20, 15:30] = 0.999999
n = dim(y_train)[1]
k = dim(y_train)[2]
#ynew = linkfun(y)
sd_y = sd(y_train)
mean_y = mean(y_train)


graph_anomaly <- melt(y_anomaly)
def_graph(graph_anomaly, 0, 0)
def_graph_ndvi(graph_anomaly, 0, 0)
graph_full <- melt(y_true_selected)
def_graph(graph_full) 
graph_anomaly_yellow_lines <- melt(aperm(y_anomaly[nrow(y_anomaly):1,], c(2, 1)))
def_graph_ndvi(graph_anomaly_yellow_lines, 0, 0) +
  geom_segment(aes(x = 1, y = 25, xend = 25, yend = 25), linetype = "dotted", color = "yellow", size = 1.5) +
  geom_segment(aes(x = 25, y = 25, xend = 25, yend = 1), linetype = "dotted", color = "yellow", size = 1.5)


p=1
q=1

if(p == 0) {
  AR <- c(0)
} else {
  AR <- seq(p)
}
if(q == 0) {
  MA <- c(0)
} else {
  MA <- seq(q)
}

if(p == 0 && q == 0){
  next
}
ar = max(AR)
ma = max(MA)
m = max(ar,ma,na.rm=T)
p1 = (ar+1)^2-1
q1 = (ma+1)^2-1
reg = c(mean_y, rep(0, p1), rep(0, q1), sd_y)
set.seed(15)
estim = try(parameter_estimation_arma(reg, y_train, AR, MA, 
                                      n, k, p1, q1))
estim
estim_chol <- arma_mle_inference(reg, y_train, AR, MA, n, k, p1, q1, chol = TRUE)
estim_chol$table


if(ar > 0 && ma > 0){
  alpha <- estim[1]
  phi = estim[2:(p1+1)]
  theta = estim[(p1+2):(p1+q1+1)]
  prec <- estim[p1+q1+2]
}else{
  if(ar > 0){
    alpha <- estim[1]
    phi = estim[2:(p1+1)]
    theta = 0
    prec <- estim[p1+2]
  }else{
    alpha <- estim[1]
    phi = 0
    theta = estim[2:(q1+1)]
    prec <- estim[q1+2]
  }
}


m = max(ar,ma)

############################# PREDIÇÃO ############################# 

n_prev = 50
k_prev = 50

y_col = y_anomaly[1:(2*m), ]
y_row = rbind(y_col, y_anomaly)[, 1:(2*m)]
y_new = rbind(y_col, y_anomaly)
y_new = cbind(y_row, y_new)
y_new = preencher_com_media(y_new)
y_withma = matrix(NA, ncol=dim(y_train)[2] + k_prev, nrow=dim(y_train)[1] + n_prev)
ynew = rbind(y_col, y_withma)
ynew = cbind(y_row, ynew)
ynew = preencher_com_media(ynew)
ynew1 = c()

eta = error = mu = matrix(0, ncol=(dim(y_train)[2] + k_prev+2*m),
                          nrow=(dim(y_train)[1] + n_prev+2*m))

# p1 = max(ar)
# p2 = max(ma)
for (i in (m+1):(n+2*m+n_prev))
{
  for (j in (m+1):(k+2*m + k_prev))
  {
    theta_ = theta
    phi_ = phi
    if(p1 == 0){
      phi_ = 0
      y_new1 = 0
    }else{
      if(length(phi) == p1){
        xx = as.vector(t(y_new[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
      } else {
        
        xx = as.vector(t(y_new[(i-ar):i,(j-ar):j]))
        y_new1 = as.matrix(xx[(length(xx)-1):1])
        
        full_ar = 1:ar
        ar_insignificant = full_ar[!full_ar %in% p]
        ar_index = matrix(TRUE, nrow(y_new1))
        
        for(not_ar in ar_insignificant){
          first_discarted_index = sum(1:(not_ar-1))*2 + (not_ar-1) + 1
          last_discarted_index = sum(1:not_ar)*2 + (not_ar)
          ar_index[(first_discarted_index): last_discarted_index] = FALSE
        }
        
        y_new1 = as.matrix(y_new1[ar_index])
      }
    }
    
    
    if(q1 == 0){
      theta_ = 0
      error_new = 0
    }else{
      if(length(theta) == q1){
        
        yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
        error_new = as.matrix(yy[(length(yy)-1):1])
        
        if(length(theta) == 0){
          theta_ = 0
        }
      } else{
        yy = as.vector(t(error[(i-ma):i,(j-ma):j]))
        error_new = as.matrix(yy[(length(yy)-1):1])
        
        full_ma = 1:ma
        ma_insignificant = full_ma[!full_ma %in% 1]
        ma_index = matrix(TRUE, nrow(error_new))
        
        
        for(not_ma in ma_insignificant){
          first_discarted_index = sum(1:(not_ma-1))*2 + (not_ma-1) + 1
          last_discarted_index = sum(1:not_ma)*2 + (not_ma)
          ma_index[(first_discarted_index): last_discarted_index] = FALSE
        }
        
        error_new = as.matrix(error_new[ma_index])
      }
    }
    
    eta[i,j]  = alpha + phi_%*%y_new1 + theta_%*%error_new
    error[i,j] = y_new[i,j]-eta[i,j]
    ynew[i,j] = linkinv(eta[i,j])
  }
  
}

############################# RESULTADOS ############################# 

y_prev_result = ynew[(2*m+1):(n+n_prev+2*m), (2*m+1):(k+k_prev+2*m)]

y_anomaly_removing_na = preencher_com_media(y_anomaly)

residuo = y_anomaly_removing_na - y_prev_result

soma_residuos_quadrado <- sum(residuo^2)

raiz_quadrada <- sqrt(soma_residuos_quadrado)

residuo_padronizado = residuo/sd(residuo)

estimated_ndvi = ynew[(2*m+1):(n+n_prev+2*m), (2*m+1):(k+k_prev+2*m)]*(1 - (-1)) - 1

residuo_quantilico = qnorm(pnorm(y_anomaly, 
                                 y_prev_result, 
                                 prec))

residuo_quantilico[is.infinite(residuo_quantilico)] <- max(as.vector(residuo_quantilico)[is.finite(as.vector(residuo_quantilico))])

w1<-5
h1<-5

graph_residuo_padronizado <- melt(aperm(residuo_padronizado, c(2, 1)))
def_graph(graph_residuo_padronizado, 0, 0) +
  theme(legend.position="none")


hist_residuo_padronizado <- melt(residuo_padronizado[51:100, 51:100])
grafico_hist = try(def_hist(hist_residuo_padronizado))
grafico_hist

hist_residuo_quant <- melt(residuo_quantilico[51:100, 51:100])
grafico_hist_quant = try(def_hist(hist_residuo_quant))
grafico_hist_quant

############################# MÉTRICAS ############################# 
mae(y_anomaly, estimated_ndvi)
mse(y_anomaly, estimated_ndvi)
mape(y_anomaly, estimated_ndvi)
r_2(y_anomaly, estimated_ndvi)
rmsle(y_anomaly, estimated_ndvi)

mae(y_train, estimated_ndvi[51:100, 51:100])
mse(y_train, estimated_ndvi[51:100, 51:100])
mape(y_train, estimated_ndvi[51:100, 51:100])
r_2(y_train, estimated_ndvi[51:100, 51:100])
rmsle(y_train, estimated_ndvi[51:100, 51:100])


############################# RESÍDUOS ############################# 
residuo_threshold<- ifelse(residuo_padronizado < -3 | residuo_padronizado > 3, 1, 0)
melt_threshold <- melt(aperm(residuo_threshold[nrow(residuo_threshold):1,], c(2, 1)))
def_threshold(melt_threshold, ar, ma)

residuo_erodido <- erode(residuo_threshold)
melt_erodido <- melt(residuo_erodido)
def_threshold(melt_erodido, ar, ma)

residuo_dilatado <- dilate(residuo_erodido)
melt_dilatado <- melt(aperm(residuo_dilatado[nrow(residuo_dilatado):1,], c(2, 1)))
def_threshold(melt_dilatado, ar, ma)


original_ndvi = y_anomaly*(1 - (-1)) - 1
melt_ndvi <- melt(aperm(original_ndvi[nrow(original_ndvi):1,], c(2, 1)))
def_graph_ndvi(melt_ndvi, ar, ma)


melt_estimated <- melt(aperm(estimated_ndvi[nrow(estimated_ndvi):1,], c(2, 1)))
def_graph_ndvi(melt_estimated, ar, ma) 

calibracao = melt(aperm(estimated_ndvi[51:100, 51:100][nrow(estimated_ndvi[51:100, 51:100]):1,], c(2, 1)))
def_graph_ndvi(calibracao)

