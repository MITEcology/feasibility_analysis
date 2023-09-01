#' Vertices of the feasibility domain
#'
#' This function computes the vertices of the simplicial cone of the 
#' feasibility domain, when it intersects the unit ball, as well as the normal 
#' vector of the great circle that contains all the S-1 remaining vertices and 
#' the origin of the unit ball (i.e., a point where all the intrinsic growth 
#' rates are equal to zero).
#'
#' @param A_int SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom pracma crossn
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#'
#' @return Dataframe with S rows and 2xS+2 columns. Each row contains the vertex
#'     information for each species, when considering the intersection
#'     of the feasibility domain and the unit ball. Specifically, each row
#'     shows:
#'    \itemize{
#'    \item The intrinsic growth rate of species i at the given vertex (r_i).
#'    \item The vertex ID.
#'    \item The vertex name.
#'    \item the normal vector of the great circle that contains all the S-1 
#'     remaining vertices and the origin of the unit ball.
#'  }
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' cone_vertices_director_vertices(A_int)
#'
#' @export

cone_vertices_director_vertices <- function(A_int){
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    dimensions <- ncol(A_int)
    all_vertices_data <- vertices_unit_ball(A_int) %>% dplyr::mutate(V_name=NA)
    all_vertices_data$V_name <- 1:dimensions
    all_vertices_r <- all_vertices_data[,(dimensions+1):ncol(all_vertices_data)]
    
    new_columns_names <- paste0("dir_vec_",1:dimensions)
    all_vertices_r[,new_columns_names] <- NA
    
    for(i in 1:dimensions){
      
      other_vertex_data <- all_vertices_r[-i,1:dimensions]
      
      director_vector_i_aux <- pracma::crossn(as.matrix(other_vertex_data))
      module_director_vector_i_aux <- sqrt(sum(director_vector_i_aux*director_vector_i_aux))
      director_vector_i <- director_vector_i_aux/module_director_vector_i_aux
      all_vertices_r[i,new_columns_names] <- director_vector_i
      
    }
    
    return(all_vertices_r)
    
  }else{
    
    cat(test_result_A_int[[2]],"\n")
    
  }
  
  
  
  
}
