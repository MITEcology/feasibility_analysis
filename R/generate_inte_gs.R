#' Generate random interaction matrix (globally stable)
#'
#' @description Generate globally stable random interaction matrices A of a community \eqn{\mathcal{S}}
#'
#' @param S Numeric, number of species
#' @param sigma Numeric, standard deviation of interaction
#' @param conne Numeric, connectivity of interaction, default to 1
#' @param dist Character, distribution of interaction, default to "norm"
#' @param mu Numeric, mean interaction, default to 0
#' @param rep Numeric, number of matrices wanted, default to 1
#'
#' @return a single matrix or a list of globally stable SxS interaction matrices A.
#'
#' @note Note that if one chooses a half-normal distribution, the interaction matrix can be interpreted as a purely competition system.
#'
#' As default, the diagonal values are all set to a_ii = -1 (e.g., self-regulation).
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a globally stable random interaction matrix
#' matA
#'
#' @export
generate_inte_gs <- function(S, sigma, conne = 1, dist = "norm", mu = 0, rep = 1) {
  i <- 1; j <- 1
  list_matA <- list()

  while (i <= 10000 * rep) {
    matA <- generate_inte_rand(S, sigma, conne, dist, mu)
    if (check_negative_definite(matA)) {
      list_matA[[j]] <- matA
      j <- j + 1
    }
    if (j == rep + 1) break
    i <-  i + 1
    if (i == 10000 * rep) stop("Cannot generate a globally stable matrix")
  }
  if (rep == 1) return(list_matA[[1]])
  else return(list_matA)
}
