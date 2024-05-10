#============
# 10 May 2024
#============

# Aim:
  # Want to know 6600*5 matrix on GPU computing time




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


C <- matrix(rnorm(6600*5), 33000, 33000)
C_gpu <- as.gpu.matrix(C, device = "cuda")
D <- matrix(rnorm(6600*5), 33000, 33000)
D_gpu <- as.gpu.matrix(D, device = "cuda")


CD <- C_gpu %*% D_gpu








