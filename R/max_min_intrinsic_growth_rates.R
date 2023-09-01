#' Maximum and minimum (feasible) intrinsic growth rates for the community's 
#' species
#'
#' This function computes a rough estimation of the species range of feasible
#' intrinsic growth rates in a community with S species, under Lotka-Volterra 
#' dynamics. To do so, for each species, we randomly select 1,000 intrinsic 
#' growth rates vectors in the border where a given species goes extinct, and 
#' then we extract the the minimum and maximum values for each component of the 
#' 1,000xS selected vectors.
#'
#' @param A_int SxS interaction matrix, where S is the number of species.
#' @param number_points_per_species number of intrinsic growth rates vectors to
#' be randomly sampled in the border of the feasibility domain where a given
#' species goes excluded.
#'
#' @return Dataframe with S rows and 3 columns. Each row contains the 
#' information for each species. Specifically, each row shows:
#'   \itemize{
#'    \item The species name.
#'    \item The minimum intrinsic growth rate estimated for such species.
#'    \item The miaximum intrinsic growth rate estimated for such species.
#'  }
#'
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' max_min_intrinsic_growth_rates(A_int, number_points_per_species = 1000)
#'
#' @export

max_min_intrinsic_growth_rates <- function(A_int,
                                           number_points_per_species = 1000){
  
  test_result_A_int <- check_A_int(A_int)
  
  if(test_result_A_int[[1]]==T){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    limit_points <- NULL
    number_species <- ncol(A_int)
    
    limits_feasibility_domain <- function(A_int, species_limiting, 
                                          number_points_per_species){
      
      num_species <- 1:(ncol(A_int))
      
      r_values <- data.frame(matrix(ncol = (2*ncol(A_int)), 
                                    nrow = number_points_per_species))
      r_names <- c(paste0("N_", num_species),paste0("r_", num_species)) 
      colnames(r_values) <- r_names
      
      for(j in 1:ncol(A_int)){
        
        r_values[,j] <- stats::runif(number_points_per_species) # uniform on [-1, 1]
        
      }
      
      r_values[,species_limiting] <- 0
      
      for(i in 1:nrow(r_values)){
        
        N_column <- matrix(as.numeric(r_values[i,1:ncol(A_int)]),
                           ncol = 1, nrow = ncol(A_int) )
        
        r_column <- (-1)*(A_int %*% as.matrix(N_column))
        
        sq_sum <- sum(r_column*r_column)
        
        r_values[i,c((ncol(A_int)+1):(2*ncol(A_int)))] <- r_column/sqrt(sq_sum)
        
      }
      
      return(r_values)
      
    }
    
    
    for(i in 1:nrow(A_int)){
      
      limit_FD_i <- limits_feasibility_domain(A_int, species_limiting=i, 
                                              number_points_per_species=1000)
      
      limit_points <- dplyr::bind_rows(limit_points, limit_FD_i)
      
    }
    
    species_names <- colnames(A_int)
    
    if(is.null(species_names)){
      species_names <- paste0("sp_",1:number_species)
    }
    
    min_intrinsic_growth_rates <- apply(limit_points,2,min)
    max_intrinsic_growth_rates <- apply(limit_points,2,max)
    
    solution <- tibble(
      species = species_names,
      min_R = (min_intrinsic_growth_rates[(number_species+1):(2*number_species)] %>% as.numeric()),
      max_R = (max_intrinsic_growth_rates[(number_species+1):(2*number_species)]  %>% as.numeric())
    )
    
    return(solution) 
    
  }else{
    
    cat(test_result_A_int[[2]],"\n")
    
  }
  
  
}
