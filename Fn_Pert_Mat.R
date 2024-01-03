#============
# 27 Nov.2023
#============

# Function:
  # This function fisrt check if an input matrix is p.d.
  # if not, it will choose the smallest sufficient perturbation
  # to adjust the matrix's p.d due to numerical instablility

# Arg:
  # M: matrix to be tested p.d. and if not perturbed

# Output:
  # if M is tested non-p.d.
  # a smallest sufficient pertub from 1e-07 1e-06 1e-05 1e-04 1e-03 1e-02 1e-01
  # will be added to the M, M + perturb will be returned
  
  # it will also return origianl M if none of the perturb 
    # can make M p.d.

Pert_Mat <- function(M){
  source("Fn_I_sparse.R")
  
  if(all(eigen(M, symmetric = T, only.values = T)$val > 0) == F){
    
    smallest_pert <- Inf # initialize the smallest pert
    for(pert in 10^seq(-7, -1, by = 1)){
      
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








