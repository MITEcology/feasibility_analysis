#' Calculate the asymmetry of the feasibility domain
#'
#' @description Calculate the asymmetry of the feasibility domain of community \eqn{\mathcal{S}}.
#'
#' @param matA Numeric, an SxS interaction matrix A
#'
#' @return the asymmetry value
#'
#' @note the asymmetry value is defined as standard deviation of \eqn{\{||v_1||, ||v_2||, \cdots, ||v_S||\}}, where \eqn{v_i} is the \eqn{i}-th column of the interaction matrix, i.e. the \eqn{i}-th vertex of the feasibility domain \eqn{D(\mathbf{S})}. Note that \eqn{v_i} is not normalized.
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_asymmetry(matA)
feasibility_asymmetry <- function(matA) {
  # standard deviation of all column lengths
  sd(apply(matA, 2, function(x) sqrt(sum(x^2))))
}
