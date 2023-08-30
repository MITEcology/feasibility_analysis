#' Generate random interaction matrix that satisfies global dynamical stability
#'
#' Generate a globally stable random interaction matrix A of a community S
#'
#' @param S Numeric, number of species
#' @param sigma Numeric, standard deviation of interaction
#' @param conne Numeric, connectivity of interaction, default to 1
#' @param dist Character, distribution of interaction, default to "norm"
#' @param mu Numeric, mean interaction, default to 0
#' @param trail Numeric, starting number of trails to generate negative definite matrix, default to 0
#'
#' @return a globally stable SxS interaction matrix A.
#'
#' @note Note that if one chooses a half-normal distribution, the interaction matrix can be interpreted as a purely competition system. As default, the diagonal values are all set to a_ii = -1 (e.g., self-regulation).
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 1, 1, "norm")  ## Generate a globally stable random interaction matrix
#' matA
generate_inte_gs <- function(S, sigma, conne = 1, dist = "norm", mu = -0.1 * sigma, trail = 0) {
  if (dist == "norm") {
    matA <- rnorm(S * S, mean = mu, sd = sigma)
  } else if (dist == "lnorm") {
    matA <- -rlnorm(S * S, mean = mu, sd = sigma)
  }
  zeroes <- sample(
    c(rep.int(1, floor(S * S * conne)), rep.int(0, (S * S - floor(S * S * conne))))
  )
  matA[which(zeroes == 0)] <- 0
  matA <- matrix(matA, ncol = S, nrow = S)
  diag(matA) <- -1
  if (check_negative_definite(matA)) {
    return(matA)
  } else if (trail < 100) {
    trail <- trail + 1
    generate_inte_gs(S, sigma, conne, dist, mu, trail)
  } else {
    stop("Error: cannot generate negative definite matrix within 100 trails")
  }
}
