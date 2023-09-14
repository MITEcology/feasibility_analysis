#' Generate random symmetric interaction matrix with speficied eigenvalues
#'
#' @param S Numeric, number of species
#' @param lamba Numberic, either specifing the vector of eigenvalues, or the method to sample eigenvalues
#' @param lmean Numeric, center of the eigenvalues range
#' @param loffset Numeric, half length of the eigenvalues range
#' @param rep Numeric, number of matrices wanted, default to 1
#'
#' @return a single matrix or a list of SxS interaction matrices A.
#' 
#' @note  This approach decompose a symmetric interaction matrix to a product of random orthogonal matrix and specified eignevalues. Thus, it can be used to explictly assign properties or the exact eigenvalues of the interaction matrix. However, the interaction matrix has to be symmetric.
#'
#' @examples
#' matA <- generate_inte_ev(5, c(-1, -2, -3, -4, -5)) #assigning eigenvalues
#' check_negative_definite(matA)
#'
#' matA <- generate_inte_ev(4, "unif", 0.5, 0.5) #uniformly sampling eigenvalues
#' check_negative_definite(matA)
#'
#' matA <- generate_inte_ev(4, "lognorm", 0.5, 0.5) #lognormal sampling eigenvalues
#' check_negative_definite(matA)
#'
#' @export
#'
generate_inte_ev <- function(S, lambda, lmean = -0.5, loffset = 0.5, rep = 1) {
  list_matA <- list()
  if (is.vector(lambda) && length(lambda) == S) {
    for (i in 1:rep) {
      matL <- diag(lambda)
      matG <- matrix(rnorm(S * S), S, S) %>% qr() %>% qr.Q()
      list_matA[[i]] <- matG %*% matL %*% t(matG)
    }
  } else if (lambda == "unif") {
    for (i in 1:rep) {
      matL <- diag(runif(S, lmean - loffset, lmean + loffset))
      matG <- matrix(rnorm(S * S), S, S) %>% qr() %>% qr.Q()
      list_matA[[i]] <- matG %*% matL %*% t(matG)
    }
  } else if (lambda == "norm") {
    for (i in 1:rep) {
      matL <- diag(rnorm(S, lmean, loffset))
      matG <- matrix(rnorm(S * S), S, S) %>% qr() %>% qr.Q()
      list_matA[[i]] <- matG %*% matL %*% t(matG)
    }
  } else if (lambda == "halfnorm") {
    for (i in 1:rep) {
      matL <- diag(-abs(rnorm(S, lmean, loffset)))
      matG <- matrix(rnorm(S * S), S, S) %>% qr() %>% qr.Q()
      list_matA[[i]] <- matG %*% matL %*% t(matG)
    }
  } else if (lambda == "lognorm") {
    for (i in 1:rep) {
      matL <- diag(-rlnorm(S, lmean, loffset))
      matG <- matrix(rnorm(S * S), S, S) %>% qr() %>% qr.Q()
      list_matA[[i]] <- matG %*% matL %*% t(matG)
    }
  } else {
    stop("lambda must be a vector of length S or one of the following sampling methods: unif, norm, halfnorm, lognorm")
  }

  if (rep == 1) return(list_matA[[1]])
  else return(list_matA)
}
