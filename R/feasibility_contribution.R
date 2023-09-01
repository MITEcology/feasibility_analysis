#' Calculate contribution to the feasibility of a community S by a species i
#'
#' @description The ratio between the feasibility with n species and with n-1 species. How much the feasibility of a community S is positively or negatively affected by a species i. Outcomes above (resp. below) 1 mean a positive (resp. negative) contribution.
#' @param matA Numeric, a nxn interaction matrix A
#'
#' @param sp Numeric, index of the focal species
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to 30
#'
#' @return Contribution of species i to the feasibility of a community
#'
#' @importFrom mvtnorm pmvnorm
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a random interaction matrix of nA=4 species
#' feasibility_contribution(matA, sp = 1) ## Contribution of species 1 to the feasibiltiy of the entire community of nA=4 species.
feasibility_contribution <- function(matA, sp, nt = 30) {
  S <- ncol(matA)

  sub_matA <- matA[-sp, -sp]
  f1 <- feasibility_community(sub_matA, nt)
  f2 <- sum(feasibility_partition(matA, nt)[c(2^S, 2^S - sp)])

  long_term_eff <- f2 / f1
  return(long_term_eff)
}
