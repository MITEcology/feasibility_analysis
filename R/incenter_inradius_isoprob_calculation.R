#' Main parameters of the maximum isotropic cap
#'
#' Given an interaction matrix, this function computes the main parameters of 
#' its maximum isotropic cap.
#'
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' 
#' @importFrom magrittr %>%
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#'
#' @return A list with 3 elements:
#' \itemize{
#'    \item The first element is a vector that contains the intrinsic growth 
#'    rates of the feasibility domain's incenter, I.
#'    \item The second element is a number that represents the great circle 
#'    distance from I to the border of the maximum  isotropic cap, theta.
#'    \item The third element is a number that shows area of the maximum 
#'    isotropic cap, Xi.
#'  }
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' incenter_inradius_isoprob_calculation(A_int)
#'
#' @export


incenter_inradius_isoprob_calculation <- function(A_int){
  
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
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
    
    inv_dir_ver_mat <- matlib::inv(dir_ver_mat)
    I_aux_mat <- inv_dir_ver_mat %*% matrix(correction_direction,nrow = dimensions, ncol = 1)
    I_aux <- as.numeric(I_aux_mat)
    sin_theta <- sqrt(1/sum(I_aux*I_aux)) # After imposing normalization
    # cos_theta <- sqrt(1-sin_theta*sin_theta)
    I_aux2 <- I_aux*sin_theta
    I <- I_aux2/sqrt(sum(I_aux2*I_aux2))
    
    dir_ver_mat_zero <- as.matrix(curves_data[,(dimensions+3):ncol(curves_data)])
    
    dir_ver_mat_zero %*% I
    
    # probability of the isotropic area: Area of the spherical cap with angle theta/ Area of full circle
    a_par <- 0.5*(dimensions-1)
    b_par <- 0.5
    probability_sp_cap <- 0.5*zipfR::Rbeta(sin_theta*sin_theta,a_par,b_par)
    
    return(list(I,asin(sin_theta),probability_sp_cap))
    
  }else{
    
    cat(test_result_A_int[[2]],"\n")
    
  }
  
}
