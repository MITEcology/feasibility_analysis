#' Get effective competition matrix
#' @description TBA
#'
#' @param matA matrix of bipartite network (?)
#'
#' @return effective competition matrix
#'
#' @examples
#' # TBA
#'
#' @export
get_effective_competition <- function(matA) {
  gram <- t(matA) %*% matA
  for (i in seq_len(ncol(gram))) {
    gram[, i] <- gram[, i] / sum(gram[, i])
  }
  diag(gram) <- 1
  return(gram)
}