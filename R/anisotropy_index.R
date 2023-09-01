#' Anisotropy index
#'
#' Given a vector whose elements are the average probabilities of extinction of
#' the species in a given community, this function calculates its anisotropy 
#' index (i.e., the relative Shannon diversity index of such probabilities).
#' 
#' @param probability_vector vector with length S whose elements are the 
#' average probabilities of extinction of the species in a given community, 
#' where S is the number of species.
#'
#' @return A non-negative number between 0 and 1 (isotropic community).
#'    
#' @examples
#' probability_vector <- c(.25,.15,.60)
#' anisotropy_index(probability_vector)
#'
#' @export
anisotropy_index <- function(probability_vector){
  
  # check the input parameter: probability_vector
  
  if(!all(class(probability_vector) %in% c("numeric"))){
    cat("Probabilities should be numbers.\n")
  }else if(all(class(probability_vector) %in% c("matrix","array"))){
    cat("The input parameter should be vector.\n")
  }else if(!all(probability_vector >=0)){
    cat("Probabilities should be non-negaive numbers.\n")
  }else if(round(sum(probability_vector),2) != 1){
    cat("The sum of all probabilities should be equal to one.\n")
  }else{
    
    # If the vector with probabilities is OK, then the index is calculated.
    
    number_species <- length(probability_vector)
    
    probability_vector_positive <- probability_vector[probability_vector > 0]
    
    log_prob <- log(probability_vector_positive)
    
    return(-sum(probability_vector_positive*log_prob)/log(number_species))
    
  }
  
}
