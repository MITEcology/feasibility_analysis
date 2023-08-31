#' Calculate the feasibility of individual species
#' 
#' @description calculate the feasibility of species \eqn{sp}'s survival inside a community \eqn{\mathcal{S}}
#' 
#' @param matA the interaction matrix of the community
#' @param sp the index of the focal species 
#' 
#' @importFrom magrittr %>%
#' 
#' @return the feasibility value of species \eqn{sp}
#' 
#' @note the feasibility of each individual species \eqn{sp} considering that no species in the community is already at equilibrium. This is the union of feasibility regions related to each species \eqn{sp}. For example, the feasibility of species \eqn{1} is the summation of the size of the feasibility regions \eqn{D(1,2,3)}, \eqn{D(1,2,\overline{3})} \eqn{D(1,\overline{2},3)} and \eqn{D(1,\overline{2},\overline{3})}. Note that \eqn{\overline{\cdot}} corresponds to the absence of a given species.
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 0.8, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_species(matA, 2)
feasibility_species <- function(matA, sp, nt = 30) {
  # add up all the conditional probabilities
  S <- ncol(matA)
  if(sp > S) stop("sp is out of range")
  feas_regions <- which(get_compo(S, 0)[, sp] == 1)
  return(feasibility_partition(matA, nt = nt)[feas_regions] %>% sum())
}