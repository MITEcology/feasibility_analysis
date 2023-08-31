#' Calculate the asymmetry of the feasibility domain of community S
#'
#' @param matA Numeric, an SxS interaction matrix A
#'
#' @return Asymmetry or standard deviation of column vectors
#'
#' @note standard deviation of $\{||v_1||, ||v_2||, \cdots, ||v_S||\}$, where $v_i$ is the $i$-th column of the interaction matrix, i.e. the $i$-th vertex of the feasibility domain $D(\mathbf{S})$. Note that $v_i$ is not normalized.
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

