#' Calculate the recovery to abundance perturbations
#'
#' @description recovery rate of community \eqn{\mathcal{S}} to abundance perturbations
#'
#' @param matA the interaction matrix
#' @param r the intrinsic growth rate
#' @param type the type of recovery rate, either "full" or "part"
#'
#' @note recovery rate is defined as the real parts of eigenvalues of the Jacobian matrix \eqn{\mathbf{A}\mathrm{diag}(N^{*})}, the largest eigenvalue \eqn{\lambda_1} is full recovery rate, the second smallest eigenvalue \eqn{\lambda_{S-1}} is partial recovery rate. (largest and smallest in terms of absolute values)
#'
#' @return the recovery rate with corresponding type
#'
#' @examples
#' matA <- generate_inte_rand(4,1,1,"norm")
#' r <- matA %*% c(1,2,1,1)
#' feasibility_recovery(matA, r, "full")
#' feasibility_recovery(matA, r, "part")
#'
#' @export
feasibility_recovery <- function(matA, r, type) {
  N <- solve(matA) %*% r
  matJ <- matA %*% diag(c(N))
  if (type == "full") {
    return(Re(eigen(matJ)$values)[1])
  }
  if (type == "part") {
    return(Re(eigen(matJ)$values)[ncol(matA) - 1])
  }
}
# cite: according to \cite{medeiros2021merging}