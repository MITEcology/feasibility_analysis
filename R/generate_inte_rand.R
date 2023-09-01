#' Generate random interaction matrix
#'
#' @description Generate random interaction matrix A of a community \eqn{\mathcal{S}}
#'
#' @param S Numeric, number of species
#' @param sigma Numeric, standard deviation of interaction
#' @param conne Numeric, connectivity of interaction, default to 1
#' @param dist Character, distribution of interaction (normal, log-normal, Half-Normal, or Uniform). Default to "norm"
#' @param mu Numeric, mean interaction, default to 0
#'
#' @return an SxS interaction matrix A.
#'
#' @note Note that if one chooses a half-normal distribution, the interaction matrix can be interpreted as a purely competition system. As default, the diagonal values are all set to a_ii = -1 (e.g., self-regulation).
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' matA
generate_inte_rand <- function(S, sigma, conne = 1, dist = "norm", mu = 0) {
  if(!is.numeric(S)){stop("S must be numerical")}
  if(!is.vector(S)){stop("S must be a single number")}
  if(length(S) > 1){stop("S must be a single number")}
  if(!is.numeric(sigma)){stop("sigma must be numerical")}
  if(!is.vector(sigma)){stop("sigma must be a single number")}
  if(length(sigma) > 1){stop("sigma must be a single number")}
  if(!is.numeric(conne)){stop("conne must be numerical")}
  if(!is.vector(conne)){stop("conne must be a single number")}
  if(length(conne) > 1){stop("conne must be a single number")}
  if(!dist %in% c("norm", "lnorm")){stop("dist must be one either norm or lnorm")}
  #Help says half normal and uniform are allowed, but I don't see this implemented in code. Check and add if necessary.
  if(!is.numeric(mu)){stop("mu must be numerical")}
  if(!is.vector(mu)){stop("mu must be a single number")}
  if(length(mu) > 1){stop("mu must be a single number")}
  
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
  return(matA)
}