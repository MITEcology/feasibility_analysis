
#' Vertices of the feasibility domain
#'
#' This function computes the vertices of the simplicial cone of the 
#' feasibility domain when it intersects the unit ball.
#'
#' @param A_int SxS interaction matrix, where S is the number of species.
#'
#' @return Dataframe with S rows and 2xS+1 columns. Each row contains the vertex
#'     information for each species, when considering the intersection
#'     of the feasibility domain and the unit ball. Specifically, each row
#'     shows:
#'   \itemize{
#'    \item The abundance of species i at the given vertex (N_i).
#'    \item The intrinsic growth rate of species i at the given vertex (r_i).
#'    \item The vertex ID.
#'  }
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' vertices_unit_ball(A_int)
#'
#' @export

vertices_unit_ball <- function(A_int){
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    dimensions <- ncol(A_int)
    
    num_species <- 1:dimensions
    
    r_values <- data.frame(matrix(ncol = (2*dimensions), nrow = dimensions))
    r_names <- c(paste0("N_", num_species),paste0("r_", num_species)) 
    colnames(r_values) <- r_names
    r_values$vertex_ID <- NA
    
    for(i in 1:dimensions){
      
      sup_A_col_i <- A_int[,i]
      
      Nv_i_aux <- sqrt(1/sum( sup_A_col_i* sup_A_col_i))
      
      # We add a zero in x_aux at position i
      if(i==1){
        Nv_i <- c(Nv_i_aux,rep(0,dimensions-1))
      }else if(i==dimensions){
        Nv_i <- c(rep(0,dimensions-1),Nv_i_aux)
      }else{
        Nv_i <- c(rep(0,i-1),Nv_i_aux,rep(0,dimensions-i))
      }
      
      Nv_i_mat <- matrix(Nv_i,nrow = dimensions, ncol = 1)
      rv_i_mat <- (-1)*A_int %*%  Nv_i_mat
      rv_i <- as.numeric(rv_i_mat)
      
      r_values[i,1:dimensions] <- Nv_i
      r_values[i,(dimensions+1):(2*dimensions)] <- rv_i
      r_values$vertex_ID[i] <- paste0("vertex ",i)
      
    }
    
    return(r_values)
    
  }else{
    
    cat(test_result_A_int[[2]],"\n")
  
  }
  
}
