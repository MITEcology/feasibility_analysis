# function that computes all the extreme points that belong to original vertexes
inside_vertex_detection <- function(A, B) {
  SpanA <- norm2(A)
  SpanB <- norm2(B)
  # to determine whether a vertex of one cone is inside another cone or not.

  inside_vertex <- list()
  l <- 1
  for (i in 1:ncol(B)) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in 1:ncol(A)) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}
