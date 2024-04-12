#==============
# 12 Apr. 2024
#==============

# Aim:
  # addressing the error encoutered when submit 059 to HPC
  
  # Error encountered:" Error in (function (self, other)  : 
    #Expected all tensors to be on the same device, 
    # but found at least two devices, cuda:0 and cpu!"
  
  

# Method:
  # 1. add device detection for tensor, and ensure all tensores on 1 GPU
  # 2. move tensor from GPU to CPU and vice versa

# Ref:
  # https://blogs.rstudio.com/ai/posts/2020-10-01-torch-network-from-scratch/#running-on-gpu





#=============
# GPU settings
#=============
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


#==============
# Testing code
#==============

device <- torch_device("cuda")
cat("device:", device, "\n")

a <- matrix(rnorm(9), 3, 3)
#a.gpu <- as.gpu.matrix(a, device = device)
#a.gpu

#--------
# Error
#--------
# Error in cpp_torch_device(type, index) : 
#Expecting a single string value: [type=externalptr; extent=1].
#Calls: as.gpu.matrix ... <Anonymous> -> <Anonymous> -> <Anonymous> -> cpp_torch_device
#Execution halted

#------------------
# Meaning and solo
#------------------
# in gpu.matrix package
# device: It indicates the device to load cuda. If not indicated, ’device’ will be set to ’cuda’
#if it is available.

# it cannot be assigned to the device value that we obtained
  # from GPU runing, which is  CUDA Version: 12.2 

#------------------
# Has to use Tensor
#------------------

a_tensor <- torch_tensor(a, device = device)

a_tensor

# torch_tensor
#0.8204  0.4109  0.9122
#0.1312 -1.0588  0.3119
#0.8899 -0.4120 -1.0433
#[ CUDAFloatType{3,3} ]


#----------------------------
# move tensor from GPU to CPU
#----------------------------

a_tsr_cpu <- a_tensor$to("cpu")

a_mat <- as.matrix(a_tsr_cpu)

b_mat <- matrix(rnorm(9), 3, 3)
ab <- a_mat %*% b_mat

print(ab)

#---------
# overhead
#---------


