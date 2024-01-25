#=============
# 25 Jan.2024
#=============

# Aim:
  # This function turns a vector of optimized pars to a list of matrices
  # different parameters

  # Arg:
    # vector: optimzed pars vector
    # all_pars_lst: the pars list with NA's reflecting data_str


vec_2_mat <- function(vector, all_pars_lst) {
  
  if(!all(is.na(all_pars_lst[[1]]) | all_pars_lst[[1]] == 0)){
    stop("all_pars_lst must consist of NA's and 0")
  }
  
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


#all_pars_lst[[1]][2, 1] <- 1
#vec_2_mat(rep(1, 4), all_pars_lst)


