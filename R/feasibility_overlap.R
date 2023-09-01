#' Calculate the feasibility overlap between two communities
#'
#' @description Feasibility overlap between two communities of the same dimension.
#'
#' @param matA Numeric, an SxS interaction matrix A
#' @param matB Numeric, an SxS interaction matrix B
#' @param raw Logical, whether to return the raw feasibility (TRUE) or the normalized feasibility (FALSE), default to TRUE
#' @param nsamples Numeric, number of MonteCarlo sampling points, default to 3000
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to TRUE
#'
#' @return A numeric value of the feasibility overlap between these two communities
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix A of dimension nA=4 species
#' matB <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix B of dimension nB=4 species
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
