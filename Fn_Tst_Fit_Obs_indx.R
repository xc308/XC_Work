#============
# 6 Mar. 2024
#============

# Aim: 
  # Function to generate Tst, Fit, Obs indx

# Arg:
  # df: the data frame contains spatial locations to be splitted into Test and Fit
        # and the all the variate observations spanning these locations (for Obs indx)
        
        # Note: for df structure see the bottom of this script

  # num_folds: the number of folds to create for C.V.

# Ref:
  #"042_Tst_Fit_Obs_indx.R"


#==========
# Function
#==========

Fn_Tst_Fit_Obs_indx <- function(df, num_folds){
  
  # Create a sequence of the number of rows of df
  n1 <- nrow(df) # 200
  sequence <- 1:n1
  
  # Define the number of folds
  num_folds <- num_folds
  
  # Use cut to create folds based on equal intervals
  # 4 folds require 5 cuts
  folds <- cut(sequence, breaks = seq(1, max(sequence), length.out = num_folds + 1), 
               labels = F, include.lowest = TRUE)
  
  
  
  Tst_indx <- list()
  Fit_indx <- list()
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
  
  
  return(list(Tst_indx = Tst_indx, Fit_indx = Fit_indx, Obs_indx = Obs_indx))
}



#==========
# Test
#==========

# data generated using 046b
#df_TW

#Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

#Tst_indx <- Tst_Fit_Obs_indx$Tst_indx
#str(Tst_indx) # List of 4
# List of 4
#$ : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:50] 51 52 53 54 55 56 57 58 59 60 ...
#$ : int [1:50] 101 102 103 104 105 106 107 108 109 110 ...
#$ : int [1:50] 151 152 153 154 155 156 157 158 159 160 ...

#Fit_indx <- Tst_Fit_Obs_indx$Fit_indx
#str(Fit_indx)
# List of 4
#$ : int [1:150] 51 52 53 54 55 56 57 58 59 60 ...
#$ : int [1:150] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:150] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:150] 1 2 3 4 5 6 7 8 9 10 ...


#Obs_indx <- Tst_Fit_Obs_indx$Obs_indx
#str(Obs_indx)
# List of 4
#$ : int [1:1150] 51 52 53 54 55 56 57 58 59 60 ...
#$ : int [1:1150] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:1150] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:1150] 1 2 3 4 5 6 7 8 9 10 ...




#=============
# df structure
#=============
# Ref: 046b
#df_TW[1:2,]
#s    smp_Y1    smp_Y2    smp_Y3     smp_Y4     smp_Y5    smp_Y6
#1 -9.95 0.1935268 -0.976847 0.5577696 0.33677533 -0.7798481 -1.374773
#2 -9.85 0.7690593 -1.201465 0.6862664 0.05829913 -0.5546527  1.494220
#Z1         Z2         Z3        Z4         Z5         Z6
#1 -0.6844408 -0.9371033  0.5554426 0.6663564 -0.3062566 -0.1402088
#2  1.0181379 -1.1651346 -0.1464337 0.1724847 -0.2614616  1.7834590

