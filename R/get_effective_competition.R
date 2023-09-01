#' Get effective competition matrix
#' @description TBA
#' 
#' @param A matrix of bipartite network (?)
#' 
#' @return effective competition matrix
#' 
#' @examples
#' TBA
#'
#' @export
get_effective_competition <- function(A){
  gram <- t(A) %*% A
  for(i in 1:ncol(gram)){
    gram[,i] <- gram[,i]/sum(gram[,i])
  }
  diag(gram) <- 1
  return(gram)
}