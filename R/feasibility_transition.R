#' Calculate the transition matrix for switching behavior given interaction matrix
#'
#' @description calculate the feasibility partition of all possible 2^n compositions of a pool of n species assuming Lotka-Volterra dynamics. Then calculate the transition matrix based on the regions vector and distances matrix with specified structual parameters.
#'
#' @param inte Numeric, interaction matrix defining the dynamics of the system
#' @param beta Numeric, structual parameter for the transition matrix, default to 1
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom pracma eye
#' @return a row-normalized transition matrix among compositions
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(3, 1, 1, "norm")  ## Generate a random interaction matrix with nA=5 species
#' feasibility_transition(matA, 1) ## Transition matrix among the 2^3=8 possible compositions of the 3-taxa system
#'
feasibility_transition <- function(inte, beta) {
  S <- ncol(inte)
  get_region <- function(inte) {
    l <- 0
    A <- inte
    B <- pracma::eye(nrow(inte), ncol(inte))
    record0 <- get_compo(S, 0)
    region <- list()
    for (l in 1 : 2^S){
      region[[l]] <- A %*% diag(record0[l, ]) + B %*% diag(1 - (record0[l, ]))
    }
    return(region)
  }
  get_centroids <- function(region) {
    if (is.list(region)) {
      purrr::map(region, get_centroids)
    } else {
      center <- region %*% c(rep(1, ncol(region)))
      centrd <- as.vector(center / sqrt(sum(center^2)))
      return(centrd)
    }
  }
  regions <- get_region(inte)
  omegas <- feasibility_partition(inte)
  centroids <- get_centroids(regions)
  gravi_mat <- matrix(0, 2^S, 2^S)
  for (li in 1:2^S) {
    for (lj in 1:2^S) {
      if (li == lj) {
        dij <- 0
      } else {
        dij <- acos(centroids[[li]] %*% centroids[[lj]])
      }
      gravi_mat[li, lj] <- omegas[li] * omegas[lj] * exp(-beta * dij)
    }
    gravi_mat[li, ] <- norm1(gravi_mat[li, ])
  }
  return(gravi_mat)
}
