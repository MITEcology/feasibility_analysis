#' Calculate the partial resistance to parameter perturbations
#'
#' @description partial resistance of community $\mathcal{S}$ to parameter perturbations
#'
#' @param matA the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm the norm to embed feasibility domain, either "l1" or "l2"
#'
#' @note full resistance is the distance between an intrinsic growth rate vector inside the feasibility domain and a border of \eqn{D(\mathcal{S})}. It can be calculated by the nearest arc length (i.e. L2 norm) from \eqn{r} to a border as: \eqn{d_b = \arccos <r(N^*), r(\mathrm{border})>}; partial resistance is similarly defined as distance between the intrinsic growth rate vector to the vertices of feasibility domain.
#' 
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#' 
#' @return a numeric value of the partial resistance
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")
#' r <- matA %*% c(runif(4, 0, 1))
#' feasibility_resistance_partial(matA, r)
#' @export
#' @seealso \code{\link{feasibility_resistance_full}}
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