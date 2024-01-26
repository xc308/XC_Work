#=============
# 23 Jan. 2024
#=============

# Aim:
  # To select and iterate over Tst_indx and Fit_indx for C.V.
  # and Obs_indx for co-kring


# Create a sequence of numbers from 1 to 200
n1 <- nrow(df)
sequence <- 1:n1

# Define the number of folds
num_folds <- 4

# Use cut to create folds based on equal intervals
folds <- cut(sequence, breaks = seq(1, max(sequence), length.out = num_folds + 1), 
             labels = F, include.lowest = TRUE)
folds 
# 50 1's, 50 2's, 50 3's, 50 4's


fold_indx <- which(folds == 1)
# [1]  1  2  3  4  5  6  7  8  9
# [10] 10 11 12 13 14 15 16 17 18
# [19] 19 20 21 22 23 24 25 26 27
# [28] 28 29 30 31 32 33 34 35 36
# [37] 37 38 39 40 41 42 43 44 45
# [46] 46 47 48 49 50
## can be Tst_indx

# so the remaining would be Fit_indx
seq(1:200)[-fold_indx]
# [1]  51  52  53  54  55  56  57
# [8]  58  59  60  61 ...
# [141] 191 192 193 194 195 196 197
# [148] 198 199 200


#------------------------------
# iterate over different folds
#------------------------------

# Aim:
  # so each time, Tst_indx will be shift from 
    # 1:50; 51:100; 101:150; 151:200
  
  # and the remaining Fit_indx can be used in 
    # construction neg_logL, optim

  # still need to think how to connect different Fit_indx
    # with obs_indx to allow for co-kring index selection

n1 <- nrow(df)

Tst_indx <- list()
Fit_indx <- list()
Obs_indx <- list()
for (i in 1:num_folds){
  Tst_indx[[i]] <- which(folds == i)
  Fit_indx[[i]] <- seq(1:n1)[-Tst_indx[[i]]]
  Obs_indx[[i]] <- c(Fit_indx[[i]], n1 + (1:n1), 2*n1 + (1:n1))
}

Obs_indx[[1]]
Obs_indx[[2]]


#=================
# Update Obs_indx 
#=================
# Aim:
  # to expand the Obs_indx as p grows
n1

Tst_indx <- list()
Fit_index <- list()
Obs_indx <- list()
for (i in 1:num_folds){
  Tst_indx[[i]] <- which(folds == i)
  Fit_indx[[i]] <- seq(1:n1)[-Tst_indx[[i]]]
  
  obs_indx <- list(Fit_indx[[i]]) # append new things on Fit_index use list is efficient
  for(l in 1:(p-1)){
    obs_indx[[l+1]] <- ((l*n1 + 1): ((l+1)*n1))
  }
  
  Obs_indx[[i]] <- unlist(obs_indx)
}


Obs_indx[[1]]















