#=============
# 17 Jan. 2024
#=============

#-------
# Aim1 :
#-------
  # set initial values in all_par_lst list of par matrices
    # with pre-defined values for each type of parameters
    # and collect into a vector for optim function par argument


all_pars_lst <- All_paras(p = 5, data = hierarchy_data)


# sub Aim:
# want to set the initial values for A, dlt, sig2, kappa
# in all_pars_lst,
# want the result is a vector, with all the A's, all the dlt's
# sig2's and kappa's to have their own initial values
# such a vector containing all the intial values
# for optim function par arg

sum(is.na(all_pars_lst[[1]]))
# [1] 6


# now AIM:
# want to replace NA's with initial values and
# collect them in one vector for optim function


Vals <- c()
ini_vals <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
for (i in 1:length(all_pars_lst)){
  value <- rep(ini_vals[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


# [1] 1.0 1.0 1.0 1.0 1.0 1.0
# [7] 0.1 0.1 0.1 0.1 0.1 0.1
# [13] 1.0 1.0 1.0 1.0 1.0 2.0
# [19] 2.0 2.0 2.0 2.0



#--------
# Aim2:
#--------
  # turn the optimized pars vector into different matrices A_mat, dlt_mat,
    #  sig2_mat, kappa_mat

  
vec_2_mat <- function(vector, all_pars_lst) {
  
  vec_indx <- 1
  
  for(lst in 1:length(all_pars_lst)){
    for (i in 1:nrow(all_pars_lst[[lst]])){
      for (j in 1:ncol(all_pars_lst[[lst]])){
        if (is.na(all_pars_lst[[lst]][i, j])){
          all_pars_lst[[lst]][i, j] <- vector[vec_indx]
          vec_indx <- vec_indx + 1
        }  
      }
    }
  }
  return(all_pars_lst)
}


vec_2_mat(Vals) 

#$A_mat
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    0    0    0    0    0
#[2,]    1    0    0    0    0
#[3,]    1    1    0    0    0
#[4,]    0    1    1    0    0
#[5,]    0    0    0    1    0

#$dlt_mat
#.    [,1] [,2] [,3] [,4] [,5]
#[1,]  0.0  0.0  0.0  0.0    0
#[2,]  0.1  0.0  0.0  0.0    0
#[3,]  0.1  0.1  0.0  0.0    0
#[4,]  0.0  0.1  0.1  0.0    0
#[5,]  0.0  0.0  0.0  0.1    0

#$sig2_mat
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    1    0    0    0    0
#[2,]    0    1    0    0    0
#[3,]    0    0    1    0    0
#[4,]    0    0    0    1    0
#[5,]    0    0    0    0    1


#-------------------------------
# Additional parameters in theta
#-------------------------------
# Aim:
  # obs level measurement error for each process tau2
  # are also parameters on top of those in all_pars_lst

# Method:
  # 1st count the total NA parameters
  # then count on top of that
SUM <- 0
for(i in 1:length(all_pars_lst)){
  s <- sum(is.na(all_pars_lst[[i]]))
  SUM <- SUM + s
  
}

# [1] 22

diag(5)








