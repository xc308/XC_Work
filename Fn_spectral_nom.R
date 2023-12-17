#=============
# 17 Dec. 2023
#=============

# Function:
  # Aim: spectral normalize an input matrix
  # Arg: a matrix to be spectral normalized
  

spectral_norm <- function(A_mat) {
  eigen_val <- eigen(A_mat, only.values = T)$val
  max_eigen_val <- max(Mod(eigen_val))
  
  scale_factor <- 1/max_eigen_val
  
  # element-wise scaling
  normed_A <- scale_factor * A_mat
  return(normed_A)
}




