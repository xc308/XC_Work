#======================
# 4 Jan., 19 Feb. 2024
#=====================

# Aim:
# Try lower Reg number to see if SIGMA_inv p.d.
# is still robust


# Function:
# Aim: modify B function to supress its condition number

# Method:
# 1. check all zero-element row and input
# 2. spectral normalization
# 3. replace diag element to reg_num

# Arg: A_mat

check_set_SpNorm_Reg <- function(A_mat, reg_num) {
  
  require(Matrix)
  source("Fn_I_sparse.R")
  # check each row of the matrix
  # if a row has all zero elements
  # set the last element of such a row 
  # with the same position element of the last row
  # to ensure B is not low-rank or singular
  
  for(i in 1:nrow(A_mat)){
    if(all(A_mat[i, ] == 0)){
      A_mat[i, ncol(A_mat)] <- A_mat[i-1, ncol(A_mat)]
    }
  }
  
  # spectral normalization
  singular_val <- svd(A_mat)$d
  max_sg_val <- max(singular_val)
  scale_factor <- 1/max_sg_val
  
  # element-wise scaling
  normed_A <- scale_factor * A_mat
  
  # diag(norm_A) set to a small regurlaize number
  Reg_diag <- I_sparse(size = nrow(normed_A), value = reg_num)
  
  normed_A <- normed_A + Reg_diag
  
  return(normed_A)
  
}

