#============
# 9 Apr. 2024
#============

# Aim:
  # funtion to construct inverse of a GPU matrix using cholesky

# Arg:
  # A_gpu: a matrix of GPUmatrix class, real, symmetry, p.d.

# Return:
  # an inversed GPU matrix obtained using cholesky solve


chol_inv_gpu <- function(A_gpu) {
  require(GPUmatrix)
  
  I <- diag(1, nrow(A_gpu), nrow(A_gpu))
  A_gpu_inv <- chol_solve(t(chol(A_gpu)), I)
  
  return(A_gpu_inv)
}


#======================
# Test on a GPU matrix
#======================

#a.gpu <- as.gpu.matrix(c(4, 1, 1, 3), 2, 2)
#Tst_sym_pd_gpu(a.gpu)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


#a.gpu %*% chol_inv_gpu(a.gpu)
# GPUmatrix
#torch_tensor
#1.0000  0.0000
#0.0000  1.0000
#[ CPUDoubleType{2,2} ]

#str(nrow(a.gpu))




