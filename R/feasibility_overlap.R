#' Feasibility overlap between communities \eqn{\mathcal{S}_1} and \eqn{\mathcal{S}_2}
#'
#' @param matA Numeric, an SxS interaction matrix A
#' @param matB Numeric, an SxS interaction matrix B
#' @param raw Logical, whether to return the raw feasibility (TRUE) or the normalized feasibility (FALSE), default to TRUE
#' @param nsamples Numeric, number of MonteCarlo sampling points, default to 3000
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to TRUE
#'
#' @return Overlap: feasiblity region of overlap
#'
#' @note interaction matrices need to have same dimension.
#'
#' Note Inside the function raw and nt can be changed.
#' @importFrom mvtnorm pmvnorm
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix A
#' matB <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix B
#' feasibility_overlap(matA,matB)
feasibility_overlap <- function(matA, matB, raw = TRUE, nsamples = 3000, nt = 10) {
  num <- nrow(matA)

  overlap_vertex <- vertex_detection(matA, matB) %>%
    cbind(vertex_detection(matB, matA)) %>%
    unique(MARGIN = 2)

  if (qr(overlap_vertex)$rank < num) {
    volume_overlap <- 0
  } else {
    volume_overlap <- tryCatch(
      {
        replicate(nt, calculate_omega(overlap_vertex, raw, nsamples)) %>% mean()
      },
      error = function(cond) {
        0
      }
    )
  }
  return(volume_overlap)
}
