#' Anisotropy metrics for a given interaction matrix
#'
#' Given an interaction matrix, this function uses the output of the function
#' boot_prob_excl_Omega_raw to estimate the mean values and 95% confidence
#' intervals of the following metrics: small omega (i.e., the S root of capital
#' Omega, that is, the proportion of the feasible parameter space inside the 
#' unit sphere) and the anisotropy index (i.e., the relative Shannon diversity
#' index of the species' probabilities of exclusion).
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
#' @import foreach
#' @import doParallel
#' @importFrom mvtnorm pmvnorm
#' @importFrom boot boot
#' @importFrom pracma crossn
#' @importFrom matlib inv
#' @importFrom zipfR Rbeta
#' @importFrom stats quantile
#' @importFrom tibble tibble
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom stats quantile
#' 
#' @return A tibble with 1 row and 6 columns, that contain the following data:
#'  \itemize{
#'    \item The mean value of small Omega.
#'    \item The lower bound of the 95% confidence interval for small Omega.
#'    \item The upper bound of the 95% confidence interval for small Omega.
#'    \item The mean value of the anisotropy index.
#'    \item The lower bound of the 95% confidence interval for the anisotropy 
#'    index.
#'    \item The upper bound of the 95% confidence interval for the anisotropy 
#'    index.
#'  } 
#'    
#' @examples
#' numCores <- detectCores()
#' numCores
#' registerDoParallel(2)
#' A_int <- -1*diag(c(3,5,7))
#' anisotropy_metrics_4_int_matrix(A_int)
#' stopImplicitCluster()
#' 
#' @references \url{https://doi.org/10.1016/j.jtbi.2018.04.030}
#'
#' @export
anisotropy_metrics_4_int_matrix <- function(A_int, 
                                            number_Omega_replicates = 1e3,
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
    
    anisotropy_agg_info <- tibble::as_tibble(data.frame(matrix(1)))
    colnames(anisotropy_agg_info) <- "interaction_matrix"
    
    small_omega_mean <- NA
    small_omega_lowerCI <- NA
    small_omega_upperCI <- NA
    anisotropy_index_mean <- NA
    anisotropy_index_lowerCI <- NA
    anisotropy_index_upperCI <- NA
    
    
    Omega_rep <- prob_excl_Omega_df[1,(1+number_boot_replicates):(2*number_boot_replicates)] %>%
      as.numeric()
    small_omega_rep <- Omega_rep^(1/nrow(A_int))
    anisotropy_rep <- NULL
    
    for(i in 1:number_boot_replicates){
      
      probability_vector <- prob_excl_Omega_df[,i+1] %>% dplyr::pull()
      
      anisotropy_rep <-  c(anisotropy_rep, anisotropy_index(probability_vector))
      
    }
    
    small_omega_rep_CI <- stats::quantile(small_omega_rep, prob=c(.025,.975)) %>% as.numeric()
    anisotropy_rep_CI <- stats::quantile(anisotropy_rep, prob=c(.025,.975)) %>% as.numeric()
    
    anisotropy_agg_info$small_omega_mean <- mean(small_omega_rep, na.rm = T)
    anisotropy_agg_info$small_omega_lowerCI <- small_omega_rep_CI[1]
    anisotropy_agg_info$small_omega_upperCI <- small_omega_rep_CI[2]
    anisotropy_agg_info$anisotropy_index_mean <- mean(anisotropy_rep, na.rm = T)
    anisotropy_agg_info$anisotropy_index_lowerCI <- anisotropy_rep_CI[1]
    anisotropy_agg_info$anisotropy_index_upperCI <- anisotropy_rep_CI[2]
    
    return(anisotropy_agg_info[,2:ncol(anisotropy_agg_info)])
    
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
