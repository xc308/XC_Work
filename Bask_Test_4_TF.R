#==============
# 29 Mar. 2024
#==============

# Aim:
  # to install "tensorflow" package for use of GPUmatrix

# Requirements:
  # python 3.8 - 3.11 
  # CUDA 11.8
  # cuDNN 8.6



#=====================
# install tensorflow
#=====================

# first get python 3.10 installed

#install.packages("devtools")
#devtools::install_github("rstudio/tensorflow")

#install.packages("remotes")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(remotes)


install.packages("Rcpp", lib = "/bask/projects/v/vjgo8416-xchen/TF")
.libPaths("/bask/projects/v/vjgo8416-xchen/TF")
library(Rcpp)


# Use the loaded Python module
reticulate::use_python("Python/3.9.5-GCCcore-10.3.0")


#remotes::install_github("rstudio/tensorflow")
remotes::install_github("rstudio/tensorflow", lib = "/bask/projects/v/vjgo8416-xchen/TF")


.libPaths("/bask/projects/v/vjgo8416-xchen/TF")
library(tensorflow)
#install_tensorflow(envname = "/bask/projects/v/vjgo8416-xchen/TF")


hello <- tf$constant("Hello TensorFlow!")
print(hello)


#===================
# install GPUmatrix
#===================
install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

# Torch tensors allowed
#Tensorflow tensors allowed
#Your Torch installation does not have CUDA tensors available. Please check the Torch requirements and installation if you want to use CUDA tensors.

a_tf <- gpu.matrix(1:9,nrow=3,ncol=3,type="tensorflow")
class(a_tf)
# [1] "gpu.matrix.tensorflow"
#attr(,"package")
#[1] "GPUmatrix"


a_tf@gm$device
# [1] "/job:localhost/replica:0/task:0/device:CPU:0"


#============
# Test on GPU
#============

tf$config$list_physical_devices("GPU")

# Expected output:
#[[1]]
#PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')

# If instead you see an empty list(), then TensorFlow is not using a GPU.









#============
# Alternative
#============
#install.packages("devtools")
#Yes
#library(devtools)
#devtools::install_github("rstudio/tensorflow")
#1




# The current release version of TensorFlow requires 
# a Python version between 3.8 and 3.11. 
# Python versions >=3.12 are not supported.
#reticulate::install_python('3.10:latest')
# or mannually install python 

