source("./2DKARMA/functions/2dkarma_base.R")

escore_reformed <- function(z, ynew, ar, ma, n, k, p1, q1, a = 0, b = 1) {
  quantities <- build_2dkarma_quantities(
    z = z,
    ynew = ynew,
    ar = ar,
    ma = ma,
    n = n,
    k = k,
    p1 = p1,
    q1 = q1
  )
  
  gamma_score <- colSums(
    quantities$vc * quantities$mu_eta * quantities$derivative_matrix
  )
  
  c(gamma_score, sum(quantities$precision_score_term))
}