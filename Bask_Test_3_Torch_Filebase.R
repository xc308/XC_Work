#=============
# 2 Apr. 2024
#=============

# Aim: try the file-based download

# Requirements:
#  NVIDIA CUDA Toolkit versions 11.6 or 11.7
# cuDNN (a version compatible with your CUDA version)



# Method:
  # first upload the two files to the hpb torch-lantern folder

#getwd()
#/Users/xchen/Library/CloudStorage/OneDrive-UniversityofExeter/XC_PhD/HPC/lantern-0.12.0+cpu+x86_64-Linux.zip
#/Users/xchen/Library/CloudStorage/OneDrive-UniversityofExeter/XC_PhD/HPC/libtorch-cxx11-abi-shared-with-deps-2.0.1+cpu.zip



#install.packages("torch", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(torch)

Sys.setenv(TORCH_URL="/bask/projects/v/vjgo8416-xchen/torch-lantern/libtorch-cxx11-abi-shared-with-deps-2.0.1+cpu.zip")
Sys.setenv(LANTERN_URL="/bask/projects/v/vjgo8416-xchen/torch-lantern/lantern-0.12.0+cpu+x86_64-Linux.zip")
Sys.setenv(LANTERN_BASE_URL = "/bask/projects/v/vjgo8416-xchen/torch-lantern")

torch::install_torch()



#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
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
