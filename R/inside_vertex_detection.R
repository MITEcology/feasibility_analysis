#' Detect inside vertex
#'
#' @description function that computes all the extreme points that belong to original vertexes
#'
#' @param matA matrix of the original vertexes
#' @param matB matrix of the extreme points
#'
#' @return a list with the extreme points that belong to original vertexes
inside_vertex_detection <- function(matA, matB) {
  SpanA <- norm2(matA)
  SpanB <- norm2(matB)

  inside_vertex <- list()
  l <- 1
  for (i in seq_len(ncol(matB))) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in seq_len(ncol(matA))) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}
