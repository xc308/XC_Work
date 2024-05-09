#===========
# 9 May 2024
#===========

# Aim:
# compare the speed of matrix multiplication using GPU
# 200*5 by 200*5 vs 4000*5 vs 4000*5


#==============
# GPU settings
#==============
#------
# torch
#------
# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)


#-----------
# GPUmatrix
#------------

#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

system("nvidia-smi")


C <- matrix(rnorm(4000*5), 20000, 20000)
C_gpu <- as.gpu.matrix(C, device = "cuda")
D <- matrix(rnorm(4000*5), 20000, 20000)
D_gpu <- as.gpu.matrix(D, device = "cuda")


CD <- C_gpu %*% D_gpu







