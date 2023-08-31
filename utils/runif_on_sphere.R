runif_on_sphere <- function(n, d, r = 1) {
  sims <- matrix(rnorm(n * d), nrow = n, ncol = d)
  r * sims / sqrt(apply(sims, 1L, crossprod))
}