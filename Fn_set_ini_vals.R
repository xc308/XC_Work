#=============
# 26 Oct. 2023
#=============

# This is a function to replace the initial values
  # of given matrix

# Args:
  # pars_mat: parameter matrix constructed using 
              source("Fn_para_mat_construct.R")
  # ini_vals: the initial vals you want set for simulatioins


Fn_set_ini_vals <- function(pars_mat, ini_vals) {
  
  for(i in 1:nrow(pars_mat)){
    for(j in 1:ncol(pars_mat)){
      if(is.na(pars_mat[i, j])){
        pars_mat[i, j] <- ini_vals
      }
    }
  }
  pars_mat
}




