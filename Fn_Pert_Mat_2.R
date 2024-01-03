#============
# 18 Dec.2023
#============
# Aim:
  # 1. since the perturb magnitude is viable as long as
  # the it's relative small compared to the large
  # matrix elements

  # 2. want to modify Pert_Mat function to allow
    # "larger" perturb quantity beyond 0.1
    # from 0.1, 0.2, ..., 1


Pert_Mat_2 <- function(M){
  source("Fn_I_sparse.R")
  
  if(all(eigen(M, symmetric = T, only.values = T)$val > 0) == F){
    
    smallest_pert <- Inf # initialize the smallest pert
    for(pert in c(10^seq(-7, -1, by = 1), seq(0.2, 2, by = 0.1))){
      
      M_pert <- M + I_sparse(size = nrow(M), value = pert)
      
      if(all(eigen(M_pert, symmetric = T, only.values = T)$val > 0)){
        cat("smallest pert:", pert, "\n")
        smallest_pert <- pert
        
        break
      } 
    }
    
    if (smallest_pert != Inf){
      return(M_pert)
    } else {
      cat("No suitable pert found.", "\n")
      cat("Min & Max singular value:", range(svd(M)$d), "\n" )
      cat("Condition number is:",kappa(M), "\n" )
      return(M)}
    
  } else {
    cat("No need to perturb.", "\n")
    return(M)}
  
}


