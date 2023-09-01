#' Parameters of the species closest to exclusion
#'
#' Given an interaction matrix and an initial vector of intrinsic growing rates,
#' denoted by center, this function computes the intrinsic growing
#' rates and abundances of the species that is/are closest to exclusion.
#'
#' @param center an S vector that contains the initial intrinsic growth 
#'    rates that are used to computed the minimal isotropic variation that
#'    excludes one or more species from the community.
#' @param A_int SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom magrittr %>%
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#'
#' @return A matrix with 2xS columns. The number of rows can take values
#' between 1 and S. Each row contains the information of one of species that
#' is closest to exclusion. Specifically, it shows:
#'    \itemize{
#'    \item The intrinsic growth rate of the species when it is excluded 
#'    (r_T_i).
#'    \item The abundance of species when it is excluded (N_T_i).
#'  }
#'
#' @examples
#' A_int <- (-1) * diag(c(3,5,7))
#' incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
#' I <- incenter_inradius_isoprob[[1]]
#' center <- I + runif(nrow(A_int))
#' closest_to_exclusion(center, A_int)
#'
#' @export

closest_to_exclusion <- function(center, A_int){
  
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
    
    acr_distance_center_side <- NULL
    
    for(i in 1:nrow(curves_data)){
      
      normal_vector_i <- correction_direction[i] * curves_data[i,c((col_number-(dimensions-1)):col_number)]
      sin_theta <- sum(normal_vector_i*center)
      acr_distance_center_side <- c(acr_distance_center_side,asin(sin_theta))
      
    }
    
    min_arc_dist <- min(acr_distance_center_side)
    vertex_min_arc_dist_index <- which(round(acr_distance_center_side,12) == round(min_arc_dist,12))
    
    tangent_points_matrix <- matrix(rep(0,dimensions*length(vertex_min_arc_dist_index)),
                                    nrow = length(vertex_min_arc_dist_index),
                                    ncol = dimensions)
    
    colnames(tangent_points_matrix) <- paste0("r_T_",1:dimensions)
    
    for(i in 1:length(vertex_min_arc_dist_index)){
      
      index <- vertex_min_arc_dist_index[i]
      normal_vector_i <- correction_direction[index] * curves_data[index,c((col_number-(dimensions-1)):col_number)]
      sin_theta <- sum(normal_vector_i*center)
      tangent_aux <- center - sin_theta*normal_vector_i %>% as.numeric()
      module_tangent_aux <- sqrt(sum(tangent_aux*tangent_aux))
      tangent_points_matrix[i,] <- tangent_aux / module_tangent_aux
      
    }
    
    N_tangent_points_matrix <- matrix(rep(0,dimensions*length(vertex_min_arc_dist_index)),
                                      nrow = length(vertex_min_arc_dist_index),
                                      ncol = dimensions)
    
    colnames(N_tangent_points_matrix) <- paste0("N_T_",1:dimensions)
    
    for(i in 1:nrow(N_tangent_points_matrix)){
      
      r_T_i <- tangent_points_matrix[i,] %>% as.numeric()
      r_T_i_column <- matrix(r_T_i, nrow = dimensions, ncol = 1)
      N_T_i_colum <- (-1)*matlib::inv(A_int) %*% r_T_i_column %>% as.numeric()
      N_tangent_points_matrix[i,] <- N_T_i_colum
      
    }
    
    return(cbind(tangent_points_matrix,N_tangent_points_matrix))
    
    
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
