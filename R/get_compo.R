#' Get compositions of a system
#'
#' @description function that generates a table of presence/absence combinations
#'
#' @param num Numeric, number of species in the system
#' @param nv Numeric, value to fill in the empty table and represents absence, default to 0
#'
#' @return the table of all possible communities (as a matrix)
get_compo <- function(num, nv = 0) {
  record <- matrix(nv, nrow = 2^num, ncol = num)
  k <- 2
  for (s in 1:num){
    for (i in 1:choose(num, s)){
      record[k, utils::combn(num, s)[, i]] <- 1
      k <- k + 1
    }
  }
  return(record)
}