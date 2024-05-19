#============
# 19 May 2024
#============

# Aim:
  # mircro benchmark the matrices multiplication time on CPU
  # with BLAS set to 36 threads

  # want to know how time scale up for pairs of matrices multipications
  # when the matrices size
  # grows from 3793 by 3793 to 6855 by 6855


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


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B



A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B



A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B


A <- matrix(rnorm(3793*3793), 3793, 3793)
# dim(A)
A <- matrix(rnorm(3793*3793), 3793, 3793)
AB <- A %*% B

