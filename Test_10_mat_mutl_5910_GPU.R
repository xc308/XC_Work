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
# 10 times of 5910 by 5910 multiply 5910 by 5910
#----------------------------------------------
C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm(5910*5910), 5910, 5910)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm(5910*5910), 5910, 5910)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu




#------------------------------
# 10 times of 5910*2 by 5910*2
#------------------------------
C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu

C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*2)^2), 5910*2, 5910*2)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu



#------------------------------
# 10 times of 5910*3 by 5910*3
#------------------------------
C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*3)^2), 5910*3, 5910*3)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


#------------------------------
# 10 times of 5910*4 by 5910*4
#------------------------------
C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4) # 1.8G
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu


C <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
C_gpu <- as.gpu.matrix(C, device = "cuda")

D <- matrix(rnorm((5910*4)^2), 5910*4, 5910*4)
D_gpu <- as.gpu.matrix(D, device = "cuda")
CD <- C_gpu %*% D_gpu















