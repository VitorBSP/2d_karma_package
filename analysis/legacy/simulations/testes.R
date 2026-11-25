source("./2DKARMA/functions/2dkarma_base.R")
source("./2DKARMA/functions/fisher_2dkarma.R")
source("./2DKARMA/functions/new_escore_2dkarma.R")
source("./2DKARMA/functions/simu_2dkarma.R")
source("./2DKARMA/functions/log_vero_2dkarma.R")
source("./2DKARMA/functions/fisher_2dkarma_vectorized.R")

fit_only_loglik <- function(start, y_sim, ar, ma, n, k, p1, q1) {
  optim(
    par = start,
    fn = loglik,
    method = "BFGS",
    control = list(fnscale = -1),
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
}

fit_with_score <- function(start, y_sim, ar, ma, n, k, p1, q1) {
  optim(
    par = start,
    fn = loglik,
    gr = escore_reformed,
    method = "BFGS",
    control = list(fnscale = -1),
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
}

regularize_fisher <- function(information_matrix, ridge = 1e-4) {
  symmetric_information <- 0.5 * (information_matrix + t(information_matrix))
  symmetric_information + diag(ridge, nrow(symmetric_information))
}

fit_with_fisher_polishing <- function(
    start,
    y_sim,
    ar,
    ma,
    n,
    k,
    p1,
    q1,
    tolerance = 1e-6,
    max_iterations = 50,
    min_step_factor = 1e-4,
    ridge = 1e-4
) {
  current_parameter <- start
  current_loglik <- loglik(
    z = current_parameter,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  for (iteration in seq_len(max_iterations)) {
    current_score <- escore_reformed(
      z = current_parameter,
      ynew = y_sim,
      ar = ar,
      ma = ma,
      n = n,
      k = k,
      p1 = p1,
      q1 = q1
    )
    
    current_information <- fisher_information_2dkarma(
      z = current_parameter,
      ynew = y_sim,
      ar = ar,
      ma = ma,
      n = n,
      k = k,
      p1 = p1,
      q1 = q1
    )
    
    current_information <- regularize_fisher(
      information_matrix = current_information,
      ridge = ridge
    )
    
    proposed_step <- tryCatch(
      solve(current_information, current_score),
      error = function(e) qr.solve(current_information, current_score)
    )
    
    step_factor <- 0.25
    accepted <- FALSE
    
    repeat {
      candidate_parameter <- current_parameter + step_factor * as.vector(proposed_step)
      
      if (is_valid_2dkarma(candidate_parameter, y_sim, ar, ma, n, k, p1, q1)) {
        candidate_loglik <- loglik(
          z = candidate_parameter,
          ynew = y_sim,
          ar = ar,
          ma = ma,
          n = n,
          k = k,
          p1 = p1,
          q1 = q1
        )
        
        if (is.finite(candidate_loglik) && candidate_loglik >= current_loglik) {
          accepted <- TRUE
          break
        }
      }
      
      step_factor <- step_factor / 2
      
      if (step_factor < min_step_factor) {
        break
      }
    }
    
    if (!accepted) {
      return(
        list(
          par = current_parameter,
          iterations = iteration - 1,
          converged = FALSE,
          value = current_loglik,
          message = "Fisher polishing interrompido."
        )
      )
    }
    
    parameter_change <- max(abs(candidate_parameter - current_parameter))
    current_parameter <- candidate_parameter
    current_loglik <- candidate_loglik
    
    if (parameter_change < tolerance) {
      return(
        list(
          par = current_parameter,
          iterations = iteration,
          converged = TRUE,
          value = current_loglik,
          message = "Convergência atingida."
        )
      )
    }
  }
  
  list(
    par = current_parameter,
    iterations = max_iterations,
    converged = FALSE,
    value = current_loglik,
    message = "Número máximo de iterações atingido."
  )
}

is_valid_2dkarma <- function(par, y_sim, ar, ma, n, k, p1, q1, eps = 1e-8) {
  quantities <- tryCatch(
    build_2dkarma_quantities(
      z = par,
      ynew = y_sim,
      ar = ar,
      ma = ma,
      n = n,
      k = k,
      p1 = p1,
      q1 = q1
    ),
    error = function(e) NULL
  )
  
  if (is.null(quantities)) {
    return(FALSE)
  }
  
  mu_ok <- all(is.finite(quantities$mu)) &&
    all(quantities$mu > eps) &&
    all(quantities$mu < 1 - eps)
  
  prec_ok <- is.finite(quantities$prec) && quantities$prec > eps
  delta_ok <- all(is.finite(quantities$delta)) && all(quantities$delta > eps)
  
  mu_ok && prec_ok && delta_ok
}


fit_with_fisher_puro_core <- function(
    start,
    y_sim,
    ar,
    ma,
    n,
    k,
    p1,
    q1,
    fisher_function,
    score_function,
    tolerance = 1e-8,
    max_iterations = 100,
    ridge = 0
) {
  current_parameter <- start
  
  for (iteration in seq_len(max_iterations)) {
    current_score <- tryCatch(
      score_function(
        z = current_parameter,
        ynew = y_sim,
        ar = ar,
        ma = ma,
        n = n,
        k = k,
        p1 = p1,
        q1 = q1
      ),
      error = function(e) {
        attr(e, "where") <- "score"
        e
      }
    )
    
    if (inherits(current_score, "error")) {
      return(list(
        par = current_parameter,
        iterations = iteration - 1,
        converged = FALSE,
        value = NA_real_,
        message = paste("Erro no score:", current_score$message)
      ))
    }
    
    current_information <- tryCatch(
      fisher_function(
        z = current_parameter,
        ynew = y_sim,
        ar = ar,
        ma = ma,
        n = n,
        k = k,
        p1 = p1,
        q1 = q1
      ),
      error = function(e) {
        attr(e, "where") <- "fisher"
        e
      }
    )
    
    if (inherits(current_information, "error")) {
      return(list(
        par = current_parameter,
        iterations = iteration - 1,
        converged = FALSE,
        value = NA_real_,
        message = paste("Erro na Fisher:", current_information$message)
      ))
    }
    
    current_information <- 0.5 * (current_information + t(current_information))
    
    if (ridge > 0) {
      current_information <- current_information +
        diag(ridge, nrow(current_information))
    }
    
    parameter_step <- tryCatch(
      solve(current_information, current_score),
      error = function(e) e
    )
    
    if (inherits(parameter_step, "error")) {
      return(list(
        par = current_parameter,
        iterations = iteration - 1,
        converged = FALSE,
        value = NA_real_,
        message = paste("Fisher singular:", parameter_step$message)
      ))
    }
    
    next_parameter <- current_parameter + as.vector(parameter_step)
    
    if (max(abs(next_parameter - current_parameter)) < tolerance) {
      return(list(
        par = next_parameter,
        iterations = iteration,
        converged = TRUE,
        value = loglik(
          z = next_parameter,
          ynew = y_sim,
          ar = ar,
          ma = ma,
          n = n,
          k = k,
          p1 = p1,
          q1 = q1
        ),
        message = "Convergência atingida."
      ))
    }
    
    current_parameter <- next_parameter
  }
  
  list(
    par = current_parameter,
    iterations = max_iterations,
    converged = FALSE,
    value = loglik(
      z = current_parameter,
      ynew = y_sim,
      ar = ar,
      ma = ma,
      n = n,
      k = k,
      p1 = p1,
      q1 = q1
    ),
    message = "Número máximo de iterações atingido."
  )
}


set.seed(123)

ar <- c(1)
ma <- c(1)
p1 <- (ar + 1)^2 - 1
q1 <- (ma + 1)^2 - 1
AR = max(ar)
MA = max(ma)


alpha = 0.30
phi = c(0.30, -0.15, 0.23)
theta = c(0.15, -0.10, 0.25)
prec = 17

true = c(alpha, phi, theta, prec)

y_sim <- simulacao(
  n = 50,
  k = 50,
  alpha = alpha,
  phi = phi,
  theta = theta,
  prec = prec,
  AR = max(ar),
  MA = max(ma)
)

linktemp <- substitute('logit')
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

var_y = sd(y_sim)^2
mean_y = mean(y_sim)
chute_prec = ((mean_y*(1-mean_y))/var_y - 1)*0.3
start <- c(linkfun(mean_y), rep(0, p1), rep(0, q1), chute_prec)

system.time(fit_a <- fit_only_loglik(start, y_sim, AR,MA, 50, 50, p1, q1))
fit_a
true
system.time(fit_b <- fit_with_score(start, y_sim, AR, MA, 50, 50, p1, q1))
fit_b
true
system.time(fit_c <- fit_with_fisher_polishing(start, y_sim, AR, MA, 50, 50, p1, q1))
fit_c
true
system.time(fit_c <- fit_with_fisher_polishing(fit_b$par, y_sim, AR, MA, 50, 50, p1, q1))
fit_c
true


system.time(
  fit_d <- fit_with_fisher_puro_core(
    start = start,
    y_sim = y_sim,
    ar = ar,
    ma = ma,
    n = 50,
    k = 50,
    p1 = p1,
    q1 = q1,
    fisher_function = fisher_information_2dkarma_vectorized,
    score_function = escore_reformed,
    tolerance = 1e-8,
    max_iterations = 100,
    ridge = 1e-6
  )
)

fit_d
true

system.time(
  fit_d <- fit_with_fisher_puro_core(
    start = fit_b$par,
    y_sim = y_sim,
    ar = ar,
    ma = ma,
    n = 50,
    k = 50,
    p1 = p1,
    q1 = q1,
    fisher_function = fisher_information_2dkarma_vectorized,
    score_function = escore_reformed,
    tolerance = 1e-8,
    max_iterations = 100,
    ridge = 1e-6
  )
)

fit_d
true


system.time(
  fit_e <- fit_with_fisher_puro_core(
    start = start,
    y_sim = y_sim,
    ar = ar,
    ma = ma,
    n = 50,
    k = 50,
    p1 = p1,
    q1 = q1,
    fisher_function = fisher_information_2dkarma,
    score_function = escore_reformed,
    tolerance = 1e-8,
    max_iterations = 100,
    ridge = 1e-6
  )
)

fit_e
true

system.time(
  fit_e <- fit_with_fisher_puro_core(
    start = fit_b$par,
    y_sim = y_sim,
    ar = ar,
    ma = ma,
    n = 50,
    k = 50,
    p1 = p1,
    q1 = q1,
    fisher_function = fisher_information_2dkarma,
    score_function = escore_reformed,
    tolerance = 1e-8,
    max_iterations = 100,
    ridge = 1e-6
  )
)

fit_e
true

compare_fishers <- function(theta0, y_sim, ar, ma, n, k, p1, q1) {
  score_value <- escore_reformed(
    z = theta0,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  fisher_compact <- fisher_information_2dkarma(
    z = theta0,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  fisher_vectorized <- fisher_information_2dkarma_vectorized(
    z = theta0,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  fisher_compact <- 0.5 * (fisher_compact + t(fisher_compact))
  fisher_vectorized <- 0.5 * (fisher_vectorized + t(fisher_vectorized))
  
  difference_matrix <- fisher_vectorized - fisher_compact
  
  list(
    max_abs_score = max(abs(score_value)),
    max_abs_diff = max(abs(difference_matrix)),
    mean_abs_diff = mean(abs(difference_matrix)),
    eig_compact = eigen(fisher_compact, symmetric = TRUE, only.values = TRUE)$values,
    eig_vectorized = eigen(fisher_vectorized, symmetric = TRUE, only.values = TRUE)$values
  )
}

comparison_fisher <- compare_fishers(
  theta0 = fit_b$par,
  y_sim = y_sim,
  ar = ar,
  ma = ma,
  n = 50,
  k = 50,
  p1 = p1,
  q1 = q1
)

comparison_fisher

one_fisher_step <- function(
    theta0,
    y_sim,
    ar,
    ma,
    n,
    k,
    p1,
    q1,
    fisher_function,
    ridge = 1e-6,
    step_factor = 1e-3
) {
  score_value <- escore_reformed(
    z = theta0,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  fisher_value <- fisher_function(
    z = theta0,
    ynew = y_sim,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  fisher_value <- 0.5 * (fisher_value + t(fisher_value))
  fisher_value <- fisher_value + diag(ridge, nrow(fisher_value))
  
  step_direction <- solve(fisher_value, score_value)
  theta1 <- theta0 + step_factor * as.vector(step_direction)
  
  list(
    theta1 = theta1,
    step_direction = step_direction,
    valid = is_valid_2dkarma(theta1, y_sim, ar, ma, n, k, p1, q1),
    loglik0 = loglik(theta0, y_sim, ar, ma, n, k, p1, q1),
    loglik1 = loglik(theta1, y_sim, ar, ma, n, k, p1, q1)
  )
}

step_compact <- one_fisher_step(
  theta0 = fit_b$par,
  y_sim = y_sim,
  ar = ar,
  ma = ma,
  n = 50,
  k = 50,
  p1 = p1,
  q1 = q1,
  fisher_function = fisher_information_2dkarma
)

step_vectorized <- one_fisher_step(
  theta0 = fit_b$par,
  y_sim = y_sim,
  ar = ar,
  ma = ma,
  n = 50,
  k = 50,
  p1 = p1,
  q1 = q1,
  fisher_function = fisher_information_2dkarma_vectorized
)

step_compact
step_vectorized

