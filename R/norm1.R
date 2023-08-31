#' Normalizes vector(s) in the L1 norm
#'
#' @description This function normalizes vector(s) in the L1 norm.
#'
#' @param a Numeric, a vector or a matrix
#'
#' @return A vector or a matrix, normalized in the L1 norm
norm1 <- function(a) {
  if(is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sum(x))
  } else {
    a <- a / sum(a)
  }
  return(a)
}