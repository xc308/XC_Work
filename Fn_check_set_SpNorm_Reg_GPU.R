#============
# 8 Apr. 2024
#============

# Aim:
  # To modify check_set_SpNorm_Reg function for GPU matrix


# Arg:
  # A_mat_pgu: a gpu matrix to be Spectral normalized
  # reg_num: the regularize number
  

check_set_SpNorm_Reg_gpu <- function(A_mat_gpu, reg_num) {
  
  #require(Matrix)
  #source("Fn_I_sparse.R")
  # check each row of the matrix
  # if a row has all zero elements
  # set the last element of such a row 
  # with the same position element of the last row
  # to ensure B is not low-rank or singular
  
  #for(i in 1:nrow(A_mat)){
  # if(all(A_mat[i, ] == 0)){
  #   A_mat[i, ncol(A_mat)] <- A_mat[i-1, ncol(A_mat)]
  # }
  #}
  
  
  # check if a matrix is all zero, if yes, assign an identity matrix
  if (all(A_mat_gpu == 0)) {
    #A_mat <- I_sparse(size = nrow(A_mat), value = 1)
    A_mat_gpu <- as.gpu.matrix(diag(1, nrow(A_mat_gpu), ncol(A_mat_gpu)))
    
  }
  
  
  # spectral normalization
  singular_val <- svd(A_mat_gpu)$d
  max_sg_val <- max(singular_val)
  scale_factor <- 1/max_sg_val
  
  # element-wise scaling
  normed_A <- scale_factor * A_mat_gpu
  
  # diag(norm_A) set to a small regurlaize number
  #Reg_diag <- I_sparse(size = nrow(normed_A), value = reg_num)
  Reg_diag <- as.gpu.matrix(diag(reg_num, nrow(normed_A), ncol(normed_A)))
  
  normed_A <- normed_A + Reg_diag
  
  return(normed_A)
  
}


#class(a.gpu) # [1] "gpu.matrix.torch"
#all(a.gpu == 0)
# [1] FALSE

#2 * a.gpu
#a.gpu + a.gpu


#====================
# Test on gpu matrix
#====================

#check_set_SpNorm_Reg_gpu(A_mat_gpu = a.gpu, reg_num = 1e-9)

# GPUmatrix
#torch_tensor
#-0.3337  0.8498  0.2597
#0.0978  0.1755  0.3933
#-0.4451 -0.4371  0.3067
#[ CPUDoubleType{3,3} ]







