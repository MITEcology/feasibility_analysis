#' Calculate feasibility of invasion to a community by a species
#'
#' @param matA Numeric, an SxS interaction matrix A
#' @param inv Numeric, index of the invading species
#' @param nt Numeric, number of sampling points, default to 5000
#'
#' @return Feasibility of the domain of positive intrinsic growth rate (denoted by \eqn{IGR_i>0}). This considers that the rest of the species (\eqn{j \neq i \in \mathcal{S}}) are already at equilibrium.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_invasion(matA, inv = 1)
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
