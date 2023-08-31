#' Normalizes vector(s) in the L2 norm
#'
#' @description This function normalizes vector(s) in the L2 norm.
#'
#' @param a Numeric, a vector or a matrix
#'
#' @return A vector or a matrix, normalized in the L2 norm
norm2 <- function(a) {
  if(is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    a <- a / sqrt(sum(a^2))
  }
  return(a)
}