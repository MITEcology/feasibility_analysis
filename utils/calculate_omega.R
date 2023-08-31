# function that computes the normalized feasibility from an interaction matrix
calculate_omega <- function(vertex, raw = FALSE, nsamples = 100,
                            method = "convex_hull") {
  num <- nrow(vertex)
  vertex <- norm2(vertex)
  
  if (method == "convex_hull") {
    set.seed(1010)
    vertex <- cbind(
      vertex,
      vertex %*% t(abs(runif_on_sphere(n = nsamples, d = ncol(vertex), r = 1)))
    )
    if (num < 5) {
      vertex <- norm2(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .05 * (num - 2))
        )
      )
    } else {
      vertex <- norm2(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .1 * (num - 2))
        )
      )
    }
    
    vertex <- cbind(vertex, rep(0, num))
    
    vol_ori <- (geometry::convhulln(t(vertex), output.options = TRUE)$vol)
    vol_ball <- (pi^(num / 2) / gamma(num / 2 + 1))
    # vol_ball <- calculate_omega(diag(num), nsamples = nsamples)
    
    omega <- ifelse(raw == FALSE,
                    (vol_ori / vol_ball)^(1 / num),
                    vol_ori / vol_ball
    )
  }
  if (method == "sphere") {
    m <- matrix(0, num, 1)
    a <- matrix(0, num, 1)
    b <- matrix(Inf, num, 1)
    d <- pmvnorm(
      lower = rep(0, num),
      upper = rep(Inf, num),
      mean = rep(0, num), sigma = solve(t(vertex) %*% vertex)
    )
    omega <- ifelse(raw == FALSE,
                    d[1]^(1 / num),
                    d[1]
    )
  }
  return(omega)
}
