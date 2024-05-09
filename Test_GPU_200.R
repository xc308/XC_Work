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

as.gpu.matrix(tau2_mat, device = "cuda")

A <- matrix(rnorm(200*5), 1000, 1000)
A_gpu <- as.gpu.matrix(A, device = "cuda")
B <- matrix(rnorm(200*5), 1000, 1000)
B_gpu <- as.gpu.matrix(B, device = "cuda")


AB <- A_gpu %*% B_gpu





