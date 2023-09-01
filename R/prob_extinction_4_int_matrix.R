
#' Probabilities of exclusion
#'
#' Given an interaction matrix, this function uses the output of the function
#' boot_prob_excl_Omega_raw to estimate the mean values and 95% confidence
#' intervals of the species' probabilities of exclusion.
#' 
#'
#' @param A_int a SxS interaction matrix, where S is the number of species.
#' @param number_Omega_replicates a number that specifies how many estimations 
#'    of Omega will be calculated by the quasi-Monte Carlo method in function
#'    boot_prob_excl_Omega_raw. By default, this parameter is set to 1,000 
#'    replicates.
#' @param number_boot_replicates a number that specifies how many bootstrap
#'    estimations of the probability of exclusion and Omega will be calculated
#'    for each species in function boot_prob_excl_Omega_raw. By default, this 
#'    parameter is set to 1,000 replicates.
#' @param use_chol_decomp boolean that specifies if the QR decomposition should
#' be used to compute Omega.
#' 
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom stats quantile
#' @import foreach
#' @import doParallel
#' @importFrom mvtnorm pmvnorm
#' @importFrom boot boot
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#'
#' @return A tibble with S rows and 4 columns, that contain the following data:
#'    \itemize{
#'    \item The species ID.
#'    \item The mean value of the species' probability of exclusion.
#'    \item The lower bound of the 95% confidence interval for the species' 
#'    probability of exclusion.
#'    \item The upper bound of the 95% confidence interval for the species' 
#'    probability of exclusion.
#'    } 
#'    
#' @examples
#' numCores <- detectCores()
#' numCores
#' registerDoParallel(2)
#' A_int <- -1*diag(c(3,5,7))
#' prob_extinction_4_int_matrix(A_int)
#' stopImplicitCluster()
#' 
#' @references \url{https://doi.org/10.1016/j.jtbi.2018.04.030}
#'
#' @export


prob_extinction_4_int_matrix <- function(A_int, number_Omega_replicates = 1e3,
                                         number_boot_replicates = 1e3,
                                         use_chol_decomp = F){
  
  # check the input parameters
  
  is.wholenumber <- function(x){
    integer_part <- round(x,0)
    decimal_part <- x - integer_part
    return(decimal_part == 0)
  }
  
  test_result_A_int <- check_A_int(A_int)
  test_omega_rep <- is.wholenumber(number_Omega_replicates)
  test_boot_rep <- is.wholenumber(number_boot_replicates)
  test_chol_decomp <- is.logical(use_chol_decomp)
  
  if((test_result_A_int[[1]]==T) & test_omega_rep & test_boot_rep & 
     test_chol_decomp){
    
    if(test_result_A_int[[2]] == "Adding species labels to interaction matrix."){
      
      colnames(A_int) <- paste0("sp_",1:nrow(A_int))
      rownames(A_int) <- paste0("sp_",1:nrow(A_int))
      
      cat(test_result_A_int[[2]],"\n")
      
    }
    
    prob_excl_Omega_df <- boot_prob_excl_Omega_raw(A_int, number_Omega_replicates,
                                                   number_boot_replicates,
                                                   use_chol_decomp = F)
    
    prob_excl_agg_info <- tibble::tibble(species = colnames(A_int))
    prob_excl_agg_info$prob_excl_mean <- NA
    prob_excl_agg_info$prob_excl_lowerCI <- NA
    prob_excl_agg_info$prob_excl_upperCI <- NA
    
    prob_excl_df <- prob_excl_Omega_df[,2:(number_boot_replicates-1)]
    
    for(i in 1:nrow(prob_excl_df)){
      
      probability_excl_Sp <- prob_excl_df[i, ] %>% as.numeric()
      
      probability_excl_CI <- stats::quantile(probability_excl_Sp, prob=c(.025,.975), na.rm = T) %>% as.numeric()
      
      prob_excl_agg_info$prob_excl_mean[i] <- mean(probability_excl_Sp, na.rm = T)
      prob_excl_agg_info$prob_excl_lowerCI[i] <- probability_excl_CI[1]
      prob_excl_agg_info$prob_excl_upperCI[i] <- probability_excl_CI[2]
      
    }
    
    return(prob_excl_agg_info)
    
  }else{
    
    if(test_result_A_int[[1]] != T){
      cat(test_result_A_int[[2]],"\n")
    }
    if(test_omega_rep != T){
      cat("number_Omega_replicates should be a whole number.","\n")
    }
    if(test_boot_rep != T){
      cat("number_boot_replicates should be a whole number.","\n")
    }
    if(test_chol_decomp != T){
      cat("use_chol_decomp should be logical, that is, TRUE or FALSE.","\n")
    }
    
  }
  
}
