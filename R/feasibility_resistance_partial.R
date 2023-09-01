#' Calculate the partial resistance of a community to parameter perturbations
#'
#' @description partial resistance of community to parameter perturbations. Partial resistance is defined as the closest distance between the intrinsic growth rate vector to any of the vertices of the feasibility domain.
#'
#' @param matA the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm the norm to embed feasibility domain, either "l1" or "l2"
#' 
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#' 
#' @return a numeric value of the partial (some but not all species) resistance
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm") ## Generate random interaction matrix of nA=4 species
#' r <- matA %*% c(runif(4, 0, 1)) ## Generate random location of the community inside the feasibility domain
#' feasibility_resistance_partial(matA, r) ## The max resistance to parameter perturbations before losing nA-1 species.
#' @export
#' @seealso feasibility_resistance_full
feasibility_resistance_partial <- function(matA, r, norm = "l2") {
  euclidean_distance <- function(a, b) {
    sqrt(sum((a - b)^2))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }

  if (norm == "l1") {
    r <- norm1(r); matA <- norm1(matA)
    distances <- 1:nrow(matA) %>% 
    map_dbl(~euclidean_distance(r, matA[, .x]))
  }
  if (norm == "l2") {
    r <- norm2(r); matA <- norm2(matA)
    distances <- 1:nrow(matA) %>% 
    map_dbl(~arc_length(r, matA[, .x]))
  }
  names(distances) <- 1:nrow(matA)
  return(distances)
}
