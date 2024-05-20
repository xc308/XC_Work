#============
# 20 May 2024
#============

# Aim:
  # Want to know the time of construction of SG_inv on GPU
  # involving 10 pairs of (r-1)n by (r-1)n, r sum up from 2 to 5



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


#torch_set_num_interop_threads(2)
#torch_set_num_threads(2)

cat("inter threads:", "\n")
torch_get_num_interop_threads()

cat("intra threads:", "\n")
torch_get_num_threads()


#-------------------------------
# check BLAS and OPENBLAS info: only 36 works
#-------------------------------
#install.packages("RhpcBLASctl")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(RhpcBLASctl)


#cat("Check Current BLAS Library", "\n")
#sessionInfo()

cat("Check the current number of BLAS threads", "\n")
blas_get_num_procs()

#blas_set_num_threads(48)

#cat("Updated BLAS threads:", "\n")
#blas_get_num_procs()


#-----------
# GPUmatrix
#------------

#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

system("nvidia-smi")


#======================================
# total matrices multiplications on GPU
#======================================

#----------------------------------------------
# 10 times of 3793 by 3793 multiply 3793 by 3793
#----------------------------------------------
C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(3793*3793), 3793, 3793)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(3793*3793), 3793, 3793)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu




#------------------------------
# 10 times of 3793*2 by 3793*2
#------------------------------
C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu

C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*2)^2), 3793*2, 3793*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu



#------------------------------
# 10 times of 3793*3 by 3793*3
#------------------------------
C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")
            
D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*3)^2), 3793*3, 3793*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


#------------------------------
# 10 times of 3793*4 by 3793*4
#------------------------------
C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4) # 1.8G
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((3793*4)^2), 3793*4, 3793*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


