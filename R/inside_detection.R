#' function that detects if a point is inside a simplex
#'
#' @description This function detects if a point is inside a simplex
#'
#' @param Span matrix of the simplex
#' @param vector vector to be tested
#'
#' @return 1 if the point is inside the simplex, 0 otherwise
inside_detection <- function(Span, vector) {
  lambda <- solve(Span, vector)
  if (sum(lambda >= -1e-10) == length(lambda)) {
    return(1)
  } else {
    return(0)
  }
}