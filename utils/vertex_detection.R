# function that computes all the extreme points
vertex_detection <- function(A, B) {
  num <- ncol(A)
  inside_vertex <- inside_vertex_detection(A, B)
  intersection_vertex <- intersection_vertex_detection(A, B)

  # combine the two vertex lists
  if (length(inside_vertex) > 0) {
    vertex <- matrix(unlist(inside_vertex), nrow = num, byrow = FALSE)
  } else {
    vertex <- matrix(0, nrow = num, ncol = 2)
  }
  if (length(intersection_vertex) > 0) {
    vertex <- cbind(vertex, matrix(unlist(intersection_vertex), nrow = num, byrow = FALSE))
  }

  # delete the points that are nonzero due to numerical error
  delete_zeroes <- c()
  for (i in 1:ncol(vertex)) {
    if (near(sum(vertex[, i]^2), 0)) {
      delete_zeroes <- c(delete_zeroes, i)
    }
  }
  if (length(delete_zeroes) > 0) vertex <- vertex[, -delete_zeroes]


  # delete the same ones
  if (length(vertex) > num) {
    for (test in 1:ncol(vertex)) {
      vertex[, test] <- norm2(vertex[, test])
    }
    delete_duplicates <- c()
    for (i in 1:(ncol(vertex) - 1)) {
      for (j in (i + 1):ncol(vertex)) {
        if (sum(near(vertex[, i], vertex[, j])) == nrow(vertex)) {
          delete_duplicates <- c(delete_duplicates, j)
        }
      }
    }
    if (length(delete_duplicates) > 0) vertex <- vertex[, -unique(delete_duplicates)]
  }
  return(vertex)
}
