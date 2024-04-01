#==============
# 29 Mar. 2024
#==============

# Aim:
  # To install "torch" package for the use of GPUmatrix

# Requirements:
  #  NVIDIA CUDA Toolkit versions 11.6 or 11.7
  # cuDNN (a version compatible with your CUDA version)


install.packages("torch", lib="/bask/projects/v/vjgo8416-xchen", force = TRUE)
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(torch)

#Sys.setenv(LD_LIBRARY_PATH = "/bask/projects/v/vjgo8416-xchen/torch/lib") # share objs in lib
#Sys.setenv(LD_LIBRARY_PATH = "/bask/projects/v/vjgo8416-xchen/torch") # share objs in lib

#Sys.setenv(TORCH_INSTALL = 1) #automatic installation of LibTorch and LibLantern in non-interactive 
torch::install_torch(timeout = 6000)


install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
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


