#' Calculate contribution to the feasibility of a community by a species
#'
#' @description The ratio between the feasibility with |S| species and with |S-1| species.
#' @param matA Numeric, an |S|x|S| interaction matrix A
#'
#' @param st Numeric, id of species
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to 30
#'
#' @return How much the feasibility of a community S is positively or negatively affected by a species i. Outcomes above (resp. below) 1 mean a positive (resp. negative) contribution.
#'
#' Note Inside the function nt can be changed.
#' @importFrom mvtnorm pmvnorm
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_contribution(matA, sp = 1)
feasibility_contribution <- function(matA, sp, nt = 30) {
  S <- ncol(matA)
  remaining_sp <- setdiff(1:S, sp)
  record0 <- get_compo(S,0)

  sub_matA <- matA[-sp, -sp]
  f1 <- feasibility_community(sub_matA, nt = nt)
  f2 <- sum(feasibility_partition(matA, nt = nt)[c(2^S, 2^S-sp)])

  long_term_eff <- f2 / f1
  return(long_term_eff)
}

