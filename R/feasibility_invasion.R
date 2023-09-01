#' Calculate feasibility of invasion to a community by a species i
#'
#' @description Estimates the feasibility that a species i can invade a community of n-1 species. Feasibility of the domain of positive intrinsic growth rate (denoted by IGR_i>0). This considers that the rest of the species (\eqn{j \neq i \in \mathcal{S}}) are already at equilibrium. Recall that this can be transformed into a probability.
#' @param matA Numeric, a nxn interaction matrix A
#' @param inv Numeric, index of the invading species
#' @param nt Numeric, number of sampling points, default to 5000
#'
#' @return Feasibility of invasion of species i to a community of n-1 species
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a random interaction matrix of nA=4 species
#' feasibility_invasion(matA, inv = 1) ## Estimates the feasibility of invasion by species 1 to the nA-1 species.
feasibility_invasion <- function(matA, inv, nt = 5000) {
  S <- nrow(matA)
  res <- c(1:S)[-inv]

  # sample r uniformly from the unit sphere
  sample_r <- runif_on_sphere(n = nt, d = S, r = 1)
  r <-  split(sample_r, row(sample_r)) %>% unname()
  spanA <- norm2(matA)

  # sample r from the geometric region of igr > 0 (assume global stability)
  check_igr <- map_dbl(r, ~inside_detection(cbind(diag(S)[, inv], spanA[, res]), .))
  res_only <- map_dbl(r, ~inside_detection(cbind(-diag(S)[, inv], spanA[, res]), .))

  # list of r inside the region
  r_inside <- r[sort(c(which(check_igr == 1), which(res_only == 1)))]

  # number of r inside the region
  n_point_inside <- length(r_inside)

  # calculate the colonization probability with positive igr
  igr <- sum(check_igr) / n_point_inside

  # output
  return(igr)
}
