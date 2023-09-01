#' @title Uniform random points on a sphere
#' @description Generate uniform random points on a sphere
#'
#' @param n Integer, the number of points to generate
#' @param d Integer, the dimension of the sphere
#' @param r Numeric, the radius of the sphere
#'
#' @return A matrix of size n * d, each row is a point on the sphere
#'
#' @export
#'
#' @examples
#' runif_on_sphere(20, 4, 1)
runif_on_sphere <- function(n, d, r = 1) {
  sims <- matrix(rnorm(n * d), nrow = n, ncol = d)
  r * sims / sqrt(apply(sims, 1L, crossprod))
}