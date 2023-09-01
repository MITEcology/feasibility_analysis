#' Parameters of the species closest to exclusion
#'
#' Given an interaction matrix and an initial vector of intrinsic growing rates,
#' denoted by center, this function computes the arc distance from each species
#' to the feasibility domain's border where only such species is excluded.
#'
#' @param center an S vector that contains the initial intrinsic growth 
#'    rates that are used to computed the minimal isotropic variation that
#'    excludes one or more species from the community.
#' @param A_int SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom magrittr %>%
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom tibble tibble
#'
#' @return A tibble data frame with S rows and 2 columns. Each row contains the 
#' species name and its great circle distance to the border where it goes
#' excluded (i.e, that distance is the absolute value of the angle between
#' the center and the tangent point in a given "edge" of the feasibility domain).
#'
#' @examples
#' A_int <- (-1) * diag(c(3,5,7))
#' incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
#' I <- incenter_inradius_isoprob[[1]]
#' center <- I + runif(nrow(A_int))
#' sp_distance_to_exclusion(center, A_int)
#'
#' @export

sp_distance_to_exclusion <- function(center, A_int){
  
  test_result_A_int <- check_A_int(A_int)
  test_result_center <- check_center(center, A_int)
  
  if(test_result_A_int[[1]]==T & test_result_center[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    if(test_result_center[[2]] == "Turning the initial growth rates into a unit vector."){
      
      center <- center/sqrt(sum(center*center))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    dimensions <- ncol(A_int)
    curves_data <- cone_vertices_director_vertices(A_int)
    
    col_number <- ncol(curves_data)
    
    dir_ver_mat <- as.matrix(curves_data[,(col_number-(dimensions-1)):col_number])
    
    # we check the direction of the director vectors
    correction_direction <- NULL
    
    for(i in 1: dimensions){
      
      if(all((dir_ver_mat %*% as.numeric(curves_data[i,(1:dimensions)]) %>% round(12))>=0)){
        correction_direction <- c(correction_direction,1)
      }else{
        correction_direction <- c(correction_direction,-1)
      }
      
    }
    
    acr_distance <- NULL
    
    for(i in 1:nrow(curves_data)){
      
      normal_vector_i <- correction_direction[i] * curves_data[i,c((col_number-(dimensions-1)):col_number)]
      sin_theta <- sum(normal_vector_i*center)
      acr_distance <- c(acr_distance,asin(sin_theta))
      
    }
    
    result_rank <- tibble(species = colnames(A_int),
                          arc_distance = abs(acr_distance))
    
    return(result_rank)
    
    
  }else{
    
    if(test_result_A_int[[2]]!=""){
      error_messages <- paste0(test_result_A_int[[2]],"\n",
                               test_result_center[[2]],"\n")
    }else{
      error_messages <- paste0(test_result_center[[2]],"\n")
    }
    
    cat(error_messages)
    
  }
  
}
