#==============
# 29 Mar. 2024
#==============

# Aim:
  # To install "torch" package for the use of GPUmatrix

# Requirements:
  #  NVIDIA CUDA Toolkit versions 11.6 or 11.7
  # cuDNN (a version compatible with your CUDA version)


#install.packages("torch", lib="/bask/projects/v/vjgo8416-xchen", force = TRUE)
#.libPaths("/bask/projects/v/vjgo8416-xchen")
#library(torch)

#Sys.setenv(LD_LIBRARY_PATH = "/bask/projects/v/vjgo8416-xchen/torch/lib") # share objs in lib (not necessary)

#Sys.setenv(TORCH_INSTALL = 1) #automatic installation of LibTorch and LibLantern in non-interactive 


#Sys.setenv(TORCH_URL="/bask/projects/v/vjgo8416-xchen/torch/https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip")
#Sys.setenv(LANTERN_URL="/bask/projects/v/vjgo8416-xchen/torch/https://storage.googleapis.com/torch-lantern-builds/binaries/refs/heads/cran/v0.12.0/latest/lantern-0.12.0+cpu+x86_64-Linux.zip")
#torch::install_torch(timeout = 6000, 
 #                    LANTERN_BASE_URL = "/bask/projects/v/vjgo8416-xchen/torch",
 #                    TORCH_INSTALL_DEBUG = 1)


# Set TORCH_HOME environment variable to the desired path
#Sys.setenv(TORCH_HOME = "/bask/projects/v/vjgo8416-xchen")
#Sys.setenv(TORCH_INSTALL_DEBUG = 1)

system("nvidia-smi")


## try binary versioin
options(timeout = 6000) # increasing timeout is recommended since we will be downloading a 2GB file.
# For Windows and Linux: "cpu", "cu117" are the only currently supported
# For MacOS the supported are: "cpu-intel" or "cpu-m1"
kind <- "cu118"
version <- available.packages()["torch","Version"]
options(repos = c(
  torch = sprintf("https://storage.googleapis.com/torch-lantern-builds/packages/%s/%s/", kind, version),
  CRAN = "https://cloud.r-project.org" # or any other from which you want to install the other R dependencies.
))

#Install the torch package
install.packages("torch", lib="/bask/projects/v/vjgo8416-xchen", force = TRUE)

# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)

# dependencies
#torch::install_torch(timeout = 6000)


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


