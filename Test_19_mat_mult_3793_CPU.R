#============
# 19 May 2024
#============

# Aim:
  # mircro benchmark the matrices multiplication time on CPU
  # with BLAS set to 36 threads

  # want to know how time scale up for pairs of matrices multipications
  # when the matrices size
  # grows from 3793 by 3793 to 6855 by 6855



# Slurm JD:
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --gpus-per-task=1
    #SBATCH --cpus-per-gpu=36

# When apptainer exec --nv --env OPENBLAS_NUM_THREADS=1 ../cuda_4.3.3.sif Rscript
  # inter threads: 
  #[1] 72
  # intra threads: 
  #[1] 36
  # Check the current number of BLAS threads 
  #[1] 1


# When apptainer exec --nv ../cuda_4.3.3.sif Rscript
  #Check the current number of BLAS threads 
  #[1] 36

# Means:
  # the BLAS will automatically parallelize over 36 threads

  # but when the code contain GPU matrix
  # the torch code has already parallelize over intra threads 36
  # so does not need BLAS to be parallelize over 36 again,
  # so set the env variable OPENBLAS threads = 1. 





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





A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B



A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B



A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
B <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B

