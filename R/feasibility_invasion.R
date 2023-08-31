#' Calculate feasibility of invasion to a community by a species
#'
#' @param matA Numeric, an SxS interaction matrix A
#'
#' @param invader Numeric, id of invader species
#' @param nt Numeric, number of sampling points, default to 5000
#'
#' @return Feasibility of the domain of positive intrinsic growth rate (denoted by IGR}_i>0). This considers that the rest of the species (\eqn{j \neq i \in \mathcal{S}}) are already at equilibrium.
#'
#' Note Inside the function nt can be changed.
#' @importFrom mvtnorm pmvnorm
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_invasion(matA, invader = 1)
feasibility_invasion <- function(matA, invader, nt = 5000) {
  S <- nrow(matA)
  res <- c(1:S)[-invader]
  
  # sample r uniformly from the unit sphere
  r <- runif_on_sphere(n = nt, d = S, r = 1) %>% split(. , row(.)) %>% unname()
  spanA <- norm2(matA)

  # sample r from the geometric region of igr > 0 (assume global stability)
  check_igr <- map_dbl(r, ~inside_detection(cbind(diag(S)[, invader], spanA[, res]), .))
  res_only <- map_dbl(r, ~inside_detection(cbind(-diag(S)[, invader], spanA[, res]), .))

  # list of r inside the region
  r_inside <- r[sort(c(which(check_igr == 1), which(res_only == 1)))]

  # number of r inside the region
  n_point_inside <- length(r_inside) 
  
  # calculate the colonization probability with positive igr 
  igr <- sum(check_igr)/n_point_inside

  # output
  return(igr)
}
