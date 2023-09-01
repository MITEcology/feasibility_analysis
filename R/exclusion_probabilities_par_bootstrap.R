
exclusion_probabilities_par_bootstrap <- function(A_int, replicates = 1e3){
  
  dimensions <- ncol(A_int)
  
  incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
  
  I <- incenter_inradius_isoprob[[1]]
  
  all_vertices_data_aux <- vertices_unit_ball(A_int)
  all_vertices_data <- all_vertices_data_aux[,c((dimensions+1):(2*dimensions))]
  
  feasibility_parts <- matrix(rep(0,nrow(A_int)*replicates),
                              nrow = nrow(A_int), ncol = replicates)
  
  names_sp <- colnames(A_int)
  
  for (i in 1:nrow(A_int)) {
    
    cat(names_sp[i],"\n")
    
    A_int_mod <- A_int
    
    A_int_mod[,i] <- -I
    
    feasibility_parts[i, ] <- Omega_bootstrap(A_int_mod, replicates)
    
  }
  
  return(feasibility_parts)
  
}
