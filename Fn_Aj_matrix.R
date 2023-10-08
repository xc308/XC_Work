#===========================
# Adjacency matrix function
#===========================

# Aim:
  # This function finds adjacency matrix corresponding to Neighbour
  # input: Neighbour indications; logic vector
  # output: Adjancency matrix (n by n), n is the # of spatial location
          # for each row i, if sj in N(si), a_ij = 1
          # otherwise, a_ij = 0



#----------------------------
# Find Neighouour indications
#----------------------------

# Data matrix n * p, each row is an observation, in R^p

# distance of between pairs of rows of data matrix

# set a Neighbour selection radius

# obtain Neighbour indices logical vector


#--------------------#
# Adjancency Matrix
#--------------------#

AJ_Mat <- function(Neighbour_index, data_matrix) {
  
  n <- nrow(data_matrix)
  Aj_mat <- matrix(0, n, n)
  
  pair_indx <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (Neighbour_index[pair_indx]) {
        Aj_mat[i, j] <- 1
        Aj_mat[j, i] <- 1
      }
      
      pair_indx <- pair_indx + 1
    }
  }
  
  return(Aj_mat)
}


#---------------
# Generate Data
#---------------

## Data 1: rnorm
set.seed(08-10-23)

n <- 100
p <- 2

data_mat <- matrix(rnorm(n*p), n, p)


## Data 2: separation lag
ds <- 0.01
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

DIST <- abs(t(outer(s, s,  FUN = "-")))
all(DIST == D_vec) # [1] TRUE



#----------------------------------
# distance among each pair of rows
#----------------------------------
Dist <- dist(data_mat) # Euclidean
quantile(Dist)

#        0%        25%        50%      75%       100% 
# 0.03889918 1.02634659 1.58693156 2.24224756 4.76265507


# for DIST
quantile(DIST)
# 0%  25%  50%  75% 100% 
# 0.0  0.3  0.6  1.0  1.9 

   
#-----------------
# Select Neighbour
#-----------------

# data1: set Neighbour radius 1.5

N_indx <- as.matrix(Dist <= 1.5)
str(N_indx) # logi [1:4950, 1]


# data2: set Neighbour radius 0.5
N_indx2 <- as.matrix(DIST <= 0.5)


#-----------------------------------
# Generate AJ matrix using function
#-----------------------------------
Aj_mat <- AJ_Mat(Neighbour_index = N_indx, data_matrix = data_mat)
str(Aj_mat)
Aj_mat[1:10, 1:10]


Aj_mat2 <- AJ_Mat(Neighbour_index = N_indx2, data_matrix = DIST)

DIST

