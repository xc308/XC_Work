#==============
# 29 Mar. 2024
#==============

# Aim:
  # To install "torch" package for the use of GPUmatrix

# Requirements:
  #  NVIDIA CUDA Toolkit versions 11.6 or 11.7
  # cuDNN (a version compatible with your CUDA version)


install.packages("torch")
library(torch)
Yes


install.packages("GPUmatrix")
library(GPUmatrix)


#==============
# installTorch
#==============
# check if torch has installed correctly

installTorch()
# [1] TRUE


#============
# gpu.matrix 
#============

# create and store a matrix in the GPU

# Returns a GPUmatrix object that can be either "gpu.matrix.tensorflow" 
# gpu.matrix.torch".

# If the gpu.matrix-class object is not sparse it will show on the console the matrix as it is. If the
#gpu.matrix is sparse, it will return to the console the position where there are number different from
#zero.

gpu_a <- gpu.matrix(1:9,nrow=3,ncol=3)
class(gpu_a)
# [1] "gpu.matrix.torch"
# attr(,"package")
gpu_a@gm$is_cuda
# [1] FALSE
gpu_a@gm
# torch_tensor
#1  4  7
#2  5  8
#3  6  9
#[ CPUIntType{3,3} ]

