#============
# 22 May 2024
#============

# Aim:
  # understand the Pytorch memory allocation for big matrices of size 
  # 4*6855 by 4*6855


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


install.packages("sys")
library(sys)



# Function to calculate tensor size in MB
tensor_size_mb <- function(tensor) {
  # Get the dimensions of the tensor
  tensor_dims <- dim(tensor)
  # Calculate the total number of elements
  numel <- prod(tensor_dims)
  # Get the size of each element in bytes
  element_size <- torch_dtype(tensor)$itemsize()
  # Calculate the total size of the tensor in MB
  size_mb <- numel * element_size / 1024^2
  return(size_mb)
}

# Function to monitor memory allocation behavior
observe_memory_allocation <- function() {
  for (i in 1:4) {
    tensor <- torch_randn(c(i*6855, i*6855), device = "cuda")
    tensor_size_mb <- tensor_size_mb(tensor)
    cat("Iteration", i, ": Allocated tensor size:", tensor_size_mb, "MB\n")
    # Free the tensor to avoid memory leakage
    rm(tensor)
  }
}

# Observe memory allocation behavior
observe_memory_allocation()


# Print value of PYTORCH_CUDA_ALLOC_CONF if set
if (!is.null(Sys.getenv("PYTORCH_CUDA_ALLOC_CONF"))) {
  cat("PYTORCH_CUDA_ALLOC_CONF:", Sys.getenv("PYTORCH_CUDA_ALLOC_CONF"), "\n")
} else {
  cat("PYTORCH_CUDA_ALLOC_CONF not set\n")
}



