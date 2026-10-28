# Dispatch to the recursion implementation used by the package.
karma2d_recursion <- function(y_link, alpha, ar, ma, order,
                              derivatives = FALSE) {
  karma2d_recursion_r(y_link, alpha, ar, ma, order, derivatives)
}
