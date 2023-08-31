#' Check if a matrix is negative definite
#'
#' @param matA Numeric, an interaction matrix
#'
#' @return TRUE if the matrix is negative definite, FALSE otherwise
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, 0, "norm")  ## Generate a random interaction matrix
#' check_negative_definite(matA)  ## Check if the matrix is negative definite
#' matB <- generate_inte_gs(4, 1, 1, 0, "norm")  ## Generate a globally stable random interaction matrix
#' check_negative_definite(matB)  ## Check if the matrix is negative definite
check_negative_definite <- function(matA) {
  all(eigen(matA + t(matA))$values < 0)
}