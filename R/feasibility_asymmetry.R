#' Calculate the asymmetry of the feasibility domain of community S with n species
#'
#' @description The asymmetry value is defined as the standard deviation of \eqn{\{||v_1||, ||v_2||, \cdots, ||v_n||\}}, where \eqn{v_i} is the \eqn{i}-th column of the interaction matrix, i.e. the \eqn{i}-th vertex of the feasibility domain. Note that \eqn{v_i} is not normalized. Think about a triangle as the feasibility domain. The asymmetry measures the variation in the legth of borders. Each of the borders represent the feasbility of the n-1 communities.
#'
#' @param matA Numeric, a nxn interaction matrix A
#'
#' @return the asymmetry value
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm") ## Generate a random matrix of n=4 species following a normal distribution with mean=1 and sd=1
#' feasibility_asymmetry(matA) ## This measures the irregularity of the feasiblity domain of the community A with n=4 species. The larger the outcome, the larger the assymetry
feasibility_asymmetry <- function(matA) {
  # standard deviation of all column lengths
  sd(apply(matA, 2, function(x) sqrt(sum(x^2))))
}
