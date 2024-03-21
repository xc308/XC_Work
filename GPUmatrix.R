#=============
# 20 Mar. 2024
#=============

# Aim:
  # Learn GPUmatrix package 
  # to turn ordinary matrix into GPU-accelerated ones
  # so as to utilize the GPU resources provided 


install.packages("torch")
library(torch)
Yes

install.packages("tensorflow")
Yes
library(tensorflow)


install.packages("GPUmatrix")
library(GPUmatrix)
# Not torch or tensorflow installed


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

#gpu_b <- gpu.matrix(1:9,nrow=3,ncol=3, type = "tensorflow")
#class(gpu_b)
# r went aborted

gpu_a@gm$is_cuda

a <- gpu.matrix(1:9,nrow=3,ncol=3, device="cpu")
class(a)
a@gm$is_cuda





