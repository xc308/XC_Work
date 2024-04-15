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

R_subset_t <- torch_t(R_subset)
cat("R_subset_t:", "\n")

R_subset_t

# torch_tensor
# -2.0243  0.2032  1.5153  0.1104
# 0.4468 -1.6359 -0.8843 -0.0247
# -1.4094 -0.0236  0.4664 -2.6884
# -0.5827 -0.6563  1.0171 -2.3827
# [ CUDAFloatType{4,4} ]


#-----------------------
# try: torch$transpose() (Fail)
#-----------------------
#cat("R_subset_t_v2:", "\n")

#R_subset_t <- R_subset$transpose()

#R_subset_t

#==========================
# forceSymmetric for Tensor
#==========================

source("Fn_forceSym_Tsr.R")

M <- matrix(c(4, 1.003, 1.002, 1), 2, 2)
M_tsr <- torch_tensor(M, device = device)
cat("M_tsr:", "\n")
M_tsr

# torch_tensor
#4.0000  1.0020
#1.0030  1.0000
#[ CUDAFloatType{2,2} ]


M_tsr_sym <- forceSym_Tsr(M_tsr)
cat("M_tsr_sym:", "\n")
M_tsr_sym 
# torch_tensor
#4.0000  1.0025
#1.0025  1.0000
#[ CUDAFloatType{2,2} ]


#==========================
# random number from N(0, 1)
#==========================

torch_randn(c(4))
# a vector of 4 elements
# torch_tensor
#1.4655
#0.1317
#0.8498
#1.4178
#[ CPUFloatType{4} ]

torch_randn(c(2, 3))
# a matrix of 2 rows 3 cols
#torch_tensor
#-1.2619 -0.4416 -0.7054
#-1.2961 -1.6491  0.5342
#[ CPUFloatType{2,3} ]



#===================
# Repeat_interleave (error encountered)
#===================

#y = torch_tensor(matrix(c(1, 2, 3, 4), ncol = 2, byrow=TRUE), device = device)
#torch_repeat_interleave(y, 2)


#========
# reshape
#========

# torch_reshape(self, shape)

b <- torch_tensor(matrix(c(0, 1, 2, 3), ncol = 2, byrow=TRUE), device = device)
torch_reshape(b, list(-1))


# -1 means reshaping the tensor b into a one-dimensional tensor, 
# but you're not explicitly specifying the size of that dimension. 
# letting PyTorch infer it based on the total number of elements 
# in the original tensor.


#============
# torch_round
#============
# torch_round(self, decimals)
# Returns a new tensor with each of the elements of input rounded 
# to the closest integer (default deci = 0)
# 


a = torch_randn(c(4))
a
# torch_tensor
#-0.5267
#0.0806
#-0.6318
#-1.8366
#[ CPUFloatType{4} ]
torch_round(a)
# -1
#0
#-1
#-2


torch_round(a, 4)
# torch_tensor
#-0.5267
#0.0806
#-0.6318
#-1.8366
#[ CPUFloatType{4} ]



#===========
# torch_sum
#===========
# torch_sum(self, dim, keepdim = FALSE, dtype = NULL)
# Returns the sum of each row of the input tensor in the given dimension dim.

a <- torch_randn(c(4, 4))
a
# torch_tensor
#0.6186  0.1394  1.6781 -0.1519
#0.7125 -0.7671 -0.3556 -1.3737
#0.8013  1.2749 -2.0088 -0.1229
#-1.4502  0.6271 -2.2641 -1.3831
#[ CPUFloatType{4,4} ]
 

cat("dim = 1", "\n")
torch_sum(a, 1) # sum across different rows
# torch_tensor
#0.6821 = 0.6186 + 0.7125 + 0.8013 + (-1.4502)
#1.2743
#-2.9504
#-3.0316
#[ CPUFloatType{4} ]



cat("dim = 2", "\n")
torch_sum(a, 2) # sum across different cols
#  2.2841 = 0.618 + 0.1394 + 1.6781 -0.1519
#-1.7839
#-0.0554
#-4.4704
#[ CPUFloatType{4} ]


#===========
# torch_svd
#===========

# torch_svd(self, some = TRUE, compute_uv = TRUE)

# singular value decomposition of a input real matrix
# input = U  diag(S)  t(V) .

# compute_uv (bool, optional) option whether to compute U and V or not
# If compute_uv is FALSE, the returned U and V matrices will be zero matrices of 
# shape (m by m) and
# (n by n) respectively. some will be ignored he

# The singular values S are returned in descending order.

a <- torch_randn(c(5, 3), device = device)

a_svd <- torch_svd(a)
#a_svd
#[[1]]
#torch_tensor
#-0.2079 -0.0886  0.2968
#0.7005 -0.5243 -0.3773
#0.3274  0.6936 -0.4111
#0.1188 -0.3923  0.2472
#-0.5871 -0.2868 -0.7345
#[ CUDAFloatType{5,3} ]

#[[2]]
#torch_tensor
#3.7027
#2.0156
#1.3863
#[ CUDAFloatType{3} ]

#[[3]]
#torch_tensor
#0.9956  0.0809  0.0466
#-0.0191 -0.3129  0.9496
#0.0914 -0.9463 -0.3100
#[ CUDAFloatType{3,3} ]

s <- a_svd[[2]]
#print(s)
# torch_tensor
#2.7430
#1.4159
#1.0300
#[ CUDAFloatType{3} ]

u <- a_svd[[1]]
v <- a_svd[[3]]

# verify

#cat("distance between a and svd:", "\n")
#torch_dist(a, torch_mm(torch_mm(u, torch_diag(s)), v$t()))

#distance between a and svd: 
#  torch_tensor
# 1.36946e-06


#==============
# torch_arrange
#==============
# torch_arange(
#start,
#end,
#step = 1,
#dtype = NULL,
#layout = NULL,
#device = NULL,
#requires_grad = FALSE
#)

# start: (Number) the starting value for the set of points. Default: 0.
# end: (Number) the ending value for the set of points
# step: (Number) the gap between each pair of adjacent points. Default: 1.
# device

torch_arange(start = 0, end = 5)
# [ CPUFloatType{6} ]

torch_arange(1, 10, step = 1.5, device = device)
# [ CUDAFloatType{7} ]



#=================
# torch_threshold_
#=================

# torch_threshold_(self, threshold, value)

torch_threshold_(a, 0.5, 0)
# torch_tensor
#0.0000  0.0000  0.0000
#1.2704  0.7653  1.2864
#0.0000  0.0000  1.3933
#0.0000  0.0000  0.0000
#0.5665  0.0000  0.0000


#==========================
# torch_repeat_interleave (still error)
#==========================
#y = torch_tensor(matrix(c(1, 2, 3, 4), ncol = 2, byrow=TRUE), device = device)
#torch_repeat_interleave(y, 2)


#===========
# torch_tril
#===========
# torch_tril(self, diagonal = 0L)
# Returns the lower triangular part of the matrix (2-D tensor)
# The lower triangular part of the matrix is defined as the elements on and below the diagonal

# If diagonal = 0, all elements on and below the main diagonal are retained.

# A positive value includes just as many diagonals above the main diagonal,
# a negative value excludes just as many diagonals below the main diagonal.

a <- torch_randn(c(4, 4), device = device)

cat("torch_tril(a):", "\n")
torch_tril(a)
# torch_tensor
#0.4209  0.0000  0.0000  0.0000
#0.5231 -1.0766  0.0000  0.0000
#-0.3282 -0.2651 -0.2431  0.0000
#0.5711  0.9287  0.6600 -0.7264
#[ CUDAFloatType{4,4} ]



cat("torch_tril(a, diagonal = 1", "\n")
torch_tril(a, diagonal = 1)
# torch_tensor
#0.4209 -0.4726  0.0000  0.0000
#0.5231 -1.0766  0.5269  0.0000
#-0.3282 -0.2651 -0.2431  0.8449
#0.5711  0.9287  0.6600 -0.7264
#[ CUDAFloatType{4,4} ]


torch_tril(a, diagonal = -1)
# torch_tensor
# 0.0000  0.0000  0.0000  0.0000
# 0.5231  0.0000  0.0000  0.0000
#-0.3282 -0.2651  0.0000  0.0000
# 0.5711  0.9287  0.6600  0.0000
#[ CUDAFloatType{4,4} ]
## pure low tri 


#============
# torch_trunc
#============

# Returns a new tensor with the truncated integer values of the elements of input

a <- torch_randn(c(4), device = device)
cat("trunc a:", "\n")
torch_trunc(a)


#=============
#torch_unbind (similar to unlist)
#=============

# Removes a tensor dimension.

#a <- torch_randn(c(4, 4), device = device)
a <- torch_tensor(matrix(1:9, nrow = 3, byrow = T), device = device)
#torch_unbind(torch_tensor(matrix(1:9, ncol = 3, byrow=TRUE)))
a
# torch_tensor
#1  2  3
#4  5  6
#7  8  9


cat("unbind:", "\n")
torch_unbind(a)

#unbind: 
#  [[1]]
#torch_tensor
#1
#2
#3
#[ CUDALongType{3} ]

#[[2]]
#torch_tensor
#4
#5
#6
#[ CUDALongType{3} ]

#[[3]]
#torch_tensor
#7
#8
#9
#[ CUDALongType{3} ]


#===========
# torch_var
#===========
# torch_var(self, dim, unbiased = TRUE, keepdim = FALSE)

# var(input, unbiased=TRUE)
  # Returns the variance of all elements in the input tenso

# var(input, dim, keepdim=False, unbiased=TRUE, out=NULL)
  # Returns the variance of each row of the input tensor in the given dimension dim

#a <- torch_randn(c(4, 4), device = device)
#torch_var(a)
# 0.925266
# [ CUDAFloatType{} ]

#cat("var at dim = 1:", "\n")
#torch_var(a, dim = 1)

#var at dim = 1: 
#  torch_tensor
#2.1873
#1.1289
#0.1728
#0.7682
#[ CUDAFloatType{4} ]


#===============
# torch_var_mean
#===============
# torch_var_mean(self, dim, unbiased = TRUE, keepdim = FALSE)


#==============
# torch_vstack
#==============

# Stack tensors in sequence vertically (row wise).
# torch_vstack(tensors)
  # tensors: (sequence of Tensors) sequence of tensors to concatenate

a <- torch_tensor(c(1, 3, 4), device = device) # a vector tensor
b <- torch_tensor(c(2, 3, 2), device = device)
#cat("str(b):", "\n")
#str(b): 
#  Float [1:3]

#cat("vstack a, b:", "\n")
#torch_vstack(list(a, b))
# torch_tensor
#1  3  4
#2  3  2
# [ CUDAFloatType{2,3} ]

# same effect as rbind(vec1, vec2) on CPU




#---------------
# Ty rbind tsr a,b: 
#---------------
#Error in rbind(a, b) : 
#  cannot coerce type 'externalptr' to vector of type 'list'
#Execution halted


## cannot directly combine two torch tensors using rbind() 
# because rbind() is designed to work with standard R objects like 
# vectors, matrices, or data frames. 

# NEED to use torch_cat(list(tsr1, tsr2), dim = 1) # by row: dim = 1


c <- torch_tensor(rbind(1, 2, 3), device = device)
d <- torch_tensor(rbind(4, 5, 6), device = device)
cat("str(d):", "\n")
str(d)
# str(d): 
#Float [1:3, 1:1]

e <- torch_tensor(matrix(1:9, nrow = 3), device = device)
cat("str(e):", "\n")
str(e)
# Long [1:3, 1:3], integer values ranging from 1 to 9,
        # doesn't require floating-point precision.


#cat("vstack c, d:", "\n")
#torch_vstack(list(c, d))


#str(c(1, 2, 3))  # num [1:3] 1 2 3, a vector
#str(rbind(1, 2, 3)) # num [1:3, 1] 1 2 3, a matrix

## for vstack:
  # if tensor is vector, then vstack will first transpose the vector, then row-wise stack
  # if tensor is matrix, then vstack will just row-wisely stack them



# Create two vectors
#vec1 <- c(1, 2, 3)
#vec2 <- c(4, 5, 6)

# Combine vectors by rows using rbind()
#str(rbind(vec1, vec2))
# num [1:2, 1:3] 1 4 2 5 3 6

#str(rbind(1, 2, 3)) # num [1:3, 1]


#============
# torch_where
#============

# torch_where(condition, self = NULL, other = NULL)
  # condition: (BoolTensor) When TRUE (nonzero), yield x, otherwise yield y
  # self: (Tensor) values selected at indices where condition is T
  # other: (Tensor) values selected at indices where condition is FALSE

# where(condition, x, y)
# Return a tensor of elements selected from either x or y, depending on condition


x <- torch_randn(c(3, 2), device = device) # condition T
y <- torch_ones(c(3, 2), device = device)

cat("condition x > 0:", "\n")
x > 0
# torch_tensor
#0  1
#0  0
#1  0
#[ CUDABoolType{3,2} ]


torch_where(x > 0, x, y)
# torch_tensor
#1.0000  1.9086
#1.0000  1.0000
#0.7810  1.0000
#[ CUDAFloatType{3,2} ]


#============
# torch_zeros
#============

#zeros(*size, out=NULL, dtype=NULL, layout=torch.strided, device=NULL, requires_grad=False)

torch_zeros(c(2, 3), device = device)
# torch_tensor
#0  0  0
#0  0  0
#[ CUDAFloatType{2,3} ]


torch_zeros(c(3), device = device)
# torch_tensor
#0
#0
#0
#[ CUDAFloatType{3} ]


#===============
# torch_nonzero
#===============

# Nonzero elements of tensors.
# torch_nonzero(self, as_list = FALSE)
  # If FALSE, the output tensor containing indices
  # If TRUE, one 1-D tensor for each dimension, 
    # containing the indices of each nonzero element along that dimension.

a <- torch_tensor(c(1, 2, 0, 0, 1, 3), device = device)
torch_nonzero(a)
# torch_tensor
#1
#2
#5
#6
#[ CUDALongType{4,1} ]


#======================================
# check symmetricity of a tensor matrix
#======================================

a <- torch_tensor(matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 1), 3, 3), device = device)
#cat("dim a:", "\n")
#dim(a)
# [1] 3 3
#str(dim(a))
#cat("length:", "\n")
#length(dim(a))
# 2

#length(dim(a)) != 2
# [1] FALSE

a_t <- torch_t(a)

#cat("torch_eq:", "\n")
a_at_eq <- torch_eq(a, a_t)


#cat("all:", "\n")
#torch_prod(a_at_eq) == 1 
# if all elements in a_at_eq product to be 1, then no element is 0
# then symmetric
  
# torch_tensor
#1
#[ CUDABoolType{1} ] 


eq_prod <- torch_prod(a_at_eq) 

#torch_eq(eq_prod, 1)
# torch_tensor
#1
#[ CUDABoolType{} ]

eq_prod_cpu <- eq_prod$cpu()

if (as.arry(eq_prod_cpu) == 1) {
  print("sym: Yes")
  
} else{
  print("sym: No")
}


#eq_prod > 0
# torch_tensor
#1
#[ CUDABoolType{1} ]

#torch_logical_not(torch_eq(eq_prod, 1), dtype = torch_int())

#eq_res <- torch_eq(eq_prod, 1)
#eq_res <- (eq_prod > 0)
#eq_res[1]
# torch_tensor
#1
#[ CUDABoolType{} ]


#torch_where(torch_eq(eq_prod, 1), print("sym: Yes"), print("sym: NO"))
#torch_where(eq_prod > 0, print("sym: Yes"), print("sym: NO"))


#source("Fn_isSymmetric_Tsr.R")
#isSymmetric_Tsr(a)


#============
# as.array()
#============

# see if as.array can be directly used on a torch tensor
#cat("as.arrray:", "\n")
#a_cpu <- a$cpu()

#a_ary_cpu <- as.array(a_cpu)

#isSymmetric(a_ary_cpu)
#[1] TRUE


#if(isSymmetric(a_ary_cpu) == T) {
#  print("sym: Yes")
#} else {
#  print("sym: No")
#}
# [1] "sym: Yes"


#a_tsr_cpu <- a_tensor$cpu()
#a_tsr_cpu_arry <- as.array(a_tsr_cpu)

#a_mat <- as.matrix(a_tsr_cpu_arry)






