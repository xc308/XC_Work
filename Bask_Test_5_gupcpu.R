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
#cpu <- torch_device("cpu")

cat("device:", as.character(device), "\n")

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


a_tsr_cpu <- a_tensor$cpu()
a_tsr_cpu_arry <- as.array(a_tsr_cpu)

a_mat <- as.matrix(a_tsr_cpu_arry)

b_mat <- matrix(rnorm(9), 3, 3)
ab <- a_mat %*% b_mat

print(ab)

#---------
# overhead
#---------

# there is overhead involved in moving data between CPU and GPU memory. 
# When you create a normal matrix a on the CPU and then convert it to 
# a tensor using torch_tensor() with a specified GPU device, the data 
# needs to be transferred from CPU memory to GPU memory. 
# This transfer operation incurs overhead, which includes the time 
# it takes to copy the data and the resources used during the transfer process.


# The overhead of transferring data between CPU and GPU depends on factors 
# such as the size of the data, the bandwidth of the system's PCI Express bus, 
# and the processing capabilities of the GPU. Generally, larger data transfers 
# incur more overhead, 

# To minimize overhead when working with tensors on the GPU, 
# it's recommended to perform as much computation as possible 
# directly on the GPU without unnecessary data transfers between 
# CPU and GPU memory. This often involves creating tensors directly 
# on the GPU whenever possible and performing operations on them 
# without moving them back to CPU memory unless required for further 
# processing or analysis.


#=========================================
# tensor on GPU in combination of lapply
#=========================================

R_mat <- matrix(rnorm(16), 4, 4)
R <- torch_tensor(R_mat, device = device)

t <- c(1, 2)
n <- 2
Subset_cols <- function(t) {
  start_col <- (t - 1) * n + 1
  end_col <- t * n
  
  result <- R[, start_col:end_col]
}

result_lst <- lapply(t, FUN = Subset_cols)
#R_subset <- do.call(cbind, result_lst)
R_subset <- torch_cat(result_lst, dim = 2)

cat("R_subset:", "\n")
R_subset
# torch_tensor
#1.0098  1.6499 -0.3109  1.5242
#-0.0963 -0.3132  1.4742 -0.8559
#-0.6514 -0.4540 -0.8400  1.4369
#-0.6229  0.0098 -1.2988 -0.6238
#[ CUDAFloatType{4,4} ]

# -2.0243  0.4468 -1.4094 -0.5827
# 0.2032 -1.6359 -0.0236 -0.6563
# 1.5153 -0.8843  0.4664  1.0171
# 0.1104 -0.0247 -2.6884 -2.3827
# [ CUDAFloatType{4,4} ]


#=================
# Tensor transpose
#=================

#R_subset_t <- torch_t(R_subset)
#cat("R_subset_t:", "\n")

#R_subset_t

# torch_tensor
# -2.0243  0.2032  1.5153  0.1104
# 0.4468 -1.6359 -0.8843 -0.0247
# -1.4094 -0.0236  0.4664 -2.6884
# -0.5827 -0.6563  1.0171 -2.3827
# [ CUDAFloatType{4,4} ]


#-----------------------
# try: torch$transpose() 
#-----------------------
cat("R_subset_t_v2:", "\n")

torch$transpose(R_subset, 1, 2)


#==========================
# forceSymmetric for Tensor
#==========================





