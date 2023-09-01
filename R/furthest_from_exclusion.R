#' Parameters of the species furthest away from exclusion
#'
#' Given an interaction matrix and an initial vector of intrinsic growing rates,
#' denoted by center, this function returns the column number/s of the 
#' interaction matrix with the species that is/are further away from exclusion.
#'
#' @param center a S vector that contains the initial intrinsic growth 
#'    rates that are used to computed the minimal isotropic variation that
#'    excludes one or more species from the community.
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#'
#' @return A vector whose length varies between 1 and S. It element contains the 
#' column index of the species that is/are furthest away from exclusion.
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
#' I <- incenter_inradius_isoprob[[1]]
#' center <- I + runif(nrow(A_int))
#' furthest_from_exclusion(center, A_int)
#'
#' @export


furthest_from_exclusion <- function(center, A_int){
  
  
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
    acr_distance_center_vertex <- NULL
    
    for(i in 1:nrow(curves_data)){
      
      vertex_i <- curves_data[i,c(1:dimensions)]
      acr_distance_center_vertex <- c(acr_distance_center_vertex,acos(sum(vertex_i*center)))
      
    }
    
    max_arc_dist <- max(acr_distance_center_vertex)
    vertex_max_arc_dist_index <- which(acr_distance_center_vertex == max_arc_dist)
    
    return(vertex_max_arc_dist_index)
    
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
