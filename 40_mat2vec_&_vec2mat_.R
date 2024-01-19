#=============
# 17 Jan. 2024
#=============

# Aim:
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





