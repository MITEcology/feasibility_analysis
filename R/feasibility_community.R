#' Calculate the feasibility of a community
#'
#' Calculate the feasibility of a community S
#'
#' @param A Numeric, an SxS interaction matrix A
#'
#' @return Feasibility ($\Omega$) of all the $|\mathcal{S}|$ species in the community (i.e., the size of $D(\mathcal{S})$). This measure typically decreases with dimension $|\mathcal{S}|$. If the matrix has positive and negative values then $\Omega \in [0,0.5]$; otherwise $\Omega \in [0,1/2^n]$---these bounds are important if the user aims to transform feasibility into a probability measure that assumes a uniform distribution of directions in parameter space.
#'
#' @note Note Inside the function nt can be changed to specify the number of replications to reduce numerical instabilities (currently set to 30).
#'
#' @export
#'
#' @examples
#' matA <- generate_inte_rand(4, 1, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_community(matA)
feasibility_community <- function(A, nt = 30, raw = TRUE) {
  S <- nrow(A)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    # out <- d[1]^(1 / S) # species level
    out <- d[1] # community level
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (all(f(A) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(A) %*% A)
    return(replicate(nt, omega(S, Sigma)) %>% mean())
  }
}
