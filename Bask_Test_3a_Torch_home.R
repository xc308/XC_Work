#=============
# 1 Apr. 2024
#=============

# Aim:
  # try torch installation

# TORCH_HOME only include the torch additional dependencies
Sys.setenv(TORCH_HOME = "/bask/projects/v/vjgo8416-xchen/torch-lantern")

#Install the torch package
install.packages("torch", lib="/bask/projects/v/vjgo8416-xchen", force = TRUE)

# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)


install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)


#==============
# installTorch
#==============
# check if torch has installed correctly

installTorch()
# [1] TRUE


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






