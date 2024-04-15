#=============
# 15 Apr. 2024
#=============

# Aim:
  # function: GPUmatrix version of forceSymmetric

forceSym_gpu <- function(gpu_mat) {
  
  sym_gpu_mat <- (gpu_mat + t(gpu_mat))/2
  
  return(sym_gpu_mat)
}


