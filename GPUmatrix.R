#=============
# 20 Mar. 2024
#=============

# Aim:
  # Learn GPUmatrix package 
  # to turn ordinary matrix into GPU-accelerated ones
  # so as to utilize the GPU resources provided 


install.packages("torch")
library(torch)

Yes

install.packages("tensorflow")
Yes
library(tensorflow)


install.packages("GPUmatrix")
library(GPUmatrix)
# Not torch or tensorflow installed


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


#gpu_b <- gpu.matrix(1:9,nrow=3,ncol=3, type = "tensorflow")
#class(gpu_b)
# r went aborted

gpu_a@gm$is_cuda

a <- gpu.matrix(1:9,nrow=3,ncol=3, device="cpu")
class(a)
a@gm$is_cuda
#[1] FALSE


mb <- matrix(rnorm(4), nrow = 2, ncol = 2)
gpu.mb <- gpu.matrix(mb, nrow = nrow(mb), ncol = ncol(mb))
class(gpu.mb)
# [1] "gpu.matrix.torch"
#attr(,"package")
#[1] "GPUmatrix"

gpu.mb@rownames
gpu.mb@colnames
gpu.mb@gm
# torch_tensor
#0.5366 -0.8550
#-0.3407 -0.5913
#[ CPUDoubleType{2,2} ]
#> mb
#       [,1]       [,2]
#[1,]  0.5365908 -0.8549864
#[2,] -0.3407146 -0.5913342
gpu.mb@type
# [1] "torch"
gpu.mb@sparse
# [1] FALSE


gpu.mb.f64 <- gpu.matrix(mb, nrow = nrow(mb), ncol = ncol(mb),
                     dtype = "float64")

gpu.mb.f64@gm
# # torch_tensor
#0.5366 -0.8550
#-0.3407 -0.5913
#[ CPUDoubleType{2,2} ]


#create a sparse
gpu.mat.sparse <- gpu.matrix(data=c(0,0,1,0,0,0),nrow = 3,ncol = 2, 
                             sparse = T)
gpu.mat.sparse@sparse



#====================================
# transpose of a gpu.matrix-class obj
#====================================

gpu.mb.t <- t(gpu.mb)
gpu.mb.t@gm
# torch_tensor
#0.5366 -0.3407
#-0.8550 -0.5913
#[ CPUDoubleType{2,2} ]


#=======
# apply
#=======
# do: Apply Functions over ’gpu.matrix-class’ margins
# return: a vector or a list of values obtained by 
          # applying a function to margins of a GPUmatrix.

# apply will return a list of length dim(x)[MARGIN].
# margin: 1 row; 2 col

apply(gpu.mb, 1, mean, simplify = F) # no simplify, return a list
# [[1]]
#[1] -0.1591978

#[[2]]
#[1] -0.4660244

apply(gpu.mb, 1, mean, simplify = T)
# [1] -0.1591978 -0.4660244


#============
# as_methods
#============
# do: as.matrix turn its argument into a matrix
    # as.list turn its argument into a list
    # as.numeric turns the input gpu.matrix into a numeric vector
    # as.vector turns the input gpu.matrix into a vector
    # is.numeric returns TRUE or FAALSE if input can be interpretable as numbers

mb_mat <- as.matrix(gpu.mb)
str(mb_mat)
# num [1:2, 1:2]
class(mb_mat)
# [1] "matrix" "array" 

is.numeric(gpu.mb) # [1] TRUE


#====================
# cbind_rbind_methods
#====================
# Usage:
  # cbind2(x, y)
  # rbind2(x, y)
  # cbind(x, y)
  # rbind(x, y)

# x, y: a gpu.matrix object or any other matrix class
  # if one of the input values is a gpu.matrix-class object, 
  # then the output will also be a gpu.matrix-class object. 


# new row
nwr <- c(1, 2)
gpu.mb.new <- rbind(gpu.mb, nwr)
class(gpu.mb.new)

gpu.mb.new@gm
# torch_tensor
#0.5366 -0.8550
#-0.3407 -0.5913
#1.0000  2.0000
#[ CPUDoubleType{3,2} ]


nwc <- c(3, 4, 5)
gpu.mb.new2 <- cbind(gpu.mb.new, nwc)
gpu.mb.new2@gm
# torch_tensor
#0.5366 -0.8550  3.0000
#-0.3407 -0.5913  4.0000
#1.0000  2.0000  5.0000
#[ CPUDoubleType{3,3} ]


#=======
# c: concatenate_gpu.matrix
#=======
# do: "combines its arguments to form a vector
# All arguments are coerced to a common type which is the type of the
    # returned value.the returned object is of type ’numeric’.

mb.vec <- c(gpu.mb)
class(mb.vec) # [1] "numeric"
str(mb.vec) # num [1:4] 0.537 -0.341 -0.855 -0.591


#========
# cor_cov
#========

# Correlation, Variance and Covariance for ’GPUmatrix’ objects

# If x and y are matrices 
  # then the covariances (or correlations) between 
  # the columns of x and the columns of y are computed.

a <- gpu.matrix(rnorm(10)) # a vector
a@gm
# [ CPUDoubleType{10,1} ]

cor(gpu.mb)
gpu.mb@gm

cor(mb_mat)


#=====
# det
#=====
# det Calculate the Determinant of a ’GPUMatrix’

# The function det and determinant internally call the corresponding function 
  # of the library torch or tensorflow 
  # (depending on the type of input gpu.matrix-class).

#If the input gpu.matrix-class object(s) are stored on the GPU, 
  # then the operations will be performed on the GPU.

# Value:
  #det returns the same output corresponding to the base function det, 
  # which is the determinant of x.
  #The returned value is a object of class numeric sotred in the cpu.

det(gpu.mb) # [1] -0.6086109
class(det(gpu.mb))
# [1] "numeric"


#=======
# diag
#=======

# "extract or replace the diagonal of a matrix
# or constructs a diagonal matrix."

# extract
diag(gpu.mb)
# [1]  0.5365908 -0.5913342

diag(gpu.mb) <- c(1, -1)
gpu.mb@gm
# torch_tensor
#1.0000 -0.8550
#-0.3407 -1.0000
#[ CPUDoubleType{2,2} ]


#=====
# dim
#=====

dim(gpu.mb)
# [1] 2 2

rownames(gpu.mb)
# NULL

colnames(gpu.mb)
# NULL

length(gpu.mb)
# [1] 4

nrow(gpu.mb)
# [1] 2


#=====
# dist
#=====
# do: Distance Matrix Computation with GPU
    # to compute the distances between the rows of a data matrix

# Value:
  # The function returns a gpu.matrix-class object 
  # with the corresponding distances between the rows
  # of the input gpu.matrix object.


mat_x <- matrix(rnorm(100), ncol = 2)
DST_mat_x <- as.matrix(dist(mat_x,diag = TRUE,upper = TRUE,method = "euclidean"))
str(DST_mat_x)
# num [1:50, 1:50]


gpu_DST_mat_x <- dist(as.gpu.matrix(mat_x),method = "euclidean")
gpu_DST_mat_x@gm
# [ CPUDoubleType{50,50} ]
nrow(gpu_DST_mat_x)
# [1] 50
ncol(gpu_DST_mat_x)
# [1] 50

Adj_m <- as.numeric(abs(gpu_DST_mat_x < 1))
class(Adj_m) # [1] "numeric"
str(Adj_m) # [1] "numeric"

gpu.Adj_m <- as.gpu.matrix(Adj_m, nrow(gpu_DST_mat_x), ncol(gpu_DST_mat_x))
class(gpu.Adj_m)
gpu.Adj_m@gm
# [ CPUDoubleType{50,50} ]

diag(gpu.Adj_m) <- 0
gpu.Adj_m@gm
# [ CPUDoubleType{50,50} ]


#==================================
# Extract elements from gpu.matrix
#==================================

# a: a gpu matrix
# a[i, j]: specifying elements to extract or replace

a_mat <- gpu.matrix(1:9,nrow=3,ncol=3)
rownames(a_mat) <- c("R1","R2","R3")
colnames(a_mat) <- c("C1","C2","C3")

a_mat[1, 3]
# GPUmatrix
#torch_tensor
#7
#[ CPUIntType{1,1} ]
#rownames: R1 
#colnames: C3 

a_mat[c(1, 2), ]
# GPUmatrix
#torch_tensor
#1  4  7
#2  5  8
#[ CPUIntType{2,3} ]
#rownames: R1 R2 
#colnames: C1 C2 C3 

a_mat[3,3] <- 100 # replace the element 3,3
a_mat[3, 3]
# GPUmatrix
#torch_tensor
#100
#[ CPUIntType{1,1} ]
#rownames: R3 
#colnames: C3 

a_mat[1,] <- c(1,2,1) # replace the first row
a_mat[1, ]
# GPUmatrix
#torch_tensor
#1  2  1
#[ CPUIntType{1,3} ]


#==============
# installTorch
#==============
# check if torch has installed correctly

installTorch()
# [1] FALSE








