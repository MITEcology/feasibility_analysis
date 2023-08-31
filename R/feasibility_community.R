#' Calculate the feasibility of a community S
#'
#' @param matA Numeric, an SxS interaction matrix A
#'
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to 30
#' @param raw Logical, whether to return the raw feasibility or the normalized feasibility, default to TRUE
#'
#' @return Feasibility \eqn{(\Omega)} of all the \eqn{|\mathcal{S}|} species in the community (i.e., the size of \eqn{D(\mathcal{S})}).
#'
#' @note This measure typically decreases with dimension \eqn{|\mathcal{S}|}. If the matrix has positive and negative values then \eqn{\Omega \in [0,0.5]}; otherwise \eqn{\Omega \in [0,1/2^n]}---these bounds are important if the user aims to transform feasibility into a probability measure that assumes a uniform distribution of directions in parameter space.
#'
#' Note Inside the function raw and nt can be changed.
#' @importFrom mvtnorm pmvnorm
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_community(matA)
feasibility_community <- function(matA, nt = 30, raw = TRUE) {
  S <- nrow(matA)
  omega <- function(S, Sigma) {
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    return(d[1])
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = TRUE)) == "matrix"
  if (all(f(matA) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(matA) %*% matA)
    return(replicate(nt, omega(S, Sigma)) %>% mean())
  }
}
