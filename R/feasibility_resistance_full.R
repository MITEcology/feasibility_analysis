#' Calculate the full resistance to parameter perturbations
#'
#' @description full resistance (all species remain extanct) of a community to parameter perturbations. Full resistance is the distance between an intrinsic growth rate vector inside the feasibility domain and a border of the domain. It can be calculated by the nearest arc length (i.e. L2 norm) from \eqn{r} to a border as: \eqn{d_b = \arccos <r(N^*), r(\mathrm{border})>}.
#'
#' @param matA the interaction matrix
#' @param r the intrinsic growth rate
#' @param norm the norm to embed feasibility domain, either "l1" or "l2"
#' @param nt the number of sampling points on the border of feasibility domain
#' @param all whether to return all resistance values from the samples or taking the minimal value as the estimated resistance
#'
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#' @importFrom utils combn
#'
#' @return a numeric value of the full resistance without having species going exrtinct.
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm") ## Generate random interaction matrix with nA=4 species
#' r <- matA %*% c(runif(4, 0, 1)) ## Generate random location inside the feasibility domain
#' feasibility_resistance_full(matA, r) ## Extent of parameter perturbations without leading to species extinctions.
#' @export
#' @seealso feasibility_resistance_partial
feasibility_resistance_full <- function(matA, r, norm = "l2", nt = 100, all = FALSE) {
  simplex_sampling <- function(m, n) {
    r <- list()
    for (j in 1:m) {
      dist <- c(sort(runif(n - 1, 0, 1)), 1)
      r[[j]] <- c(dist[1], diff(dist))
    }
    return(r)
  }
  euclidean_distance <- function(a, b) {
    sqrt(sum((a - b)^2))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }

  vertices <- combn(seq_len(nrow(matA)), nrow(matA) - 1, simplify = FALSE)
  distances_list <- list()
  for (i in seq_len(length(vertices))) {
    vertex <- vertices[[i]]
    if (norm == "l1") {
      r <- norm1(r); matA <- norm1(matA)
      distances_list[[i]] <- 1:nt %>%
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA) - 1))
          border_point <- c(matA[, vertex] %*% t)
          euclidean_distance(r, border_point)
        })
    }
    if (norm == "l2") {
      r <- norm2(r); matA <- norm2(matA)
      distances_list[[i]] <- 1:nt %>%
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA) - 1))
          border_point <- c(matA[, vertex] %*% t)
          border_point_norm <- border_point / sqrt(sum(border_point^2))
          arc_length(r, border_point_norm)
        })
    }
  }
  names(distances_list) <- unlist(lapply(vertices, paste, collapse = "_"))
  if (all) {
    return(distances_list)
  } else {
    return(sapply(distances_list, min))
  }
}
