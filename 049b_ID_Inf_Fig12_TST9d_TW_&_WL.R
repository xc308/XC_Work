#============
# 9 Mar. 2024
#============

# Aim:
  # use neg_logL_Matern function in 049 to do optimization
  

# Method:
source("049_ID_Inf_Fig12_TST9d_neg_logL_Matern.R")


#==========
# Settings
#==========

ini <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


## lower bound for each parameters, 
# NA: no lower bound
# kappa > 0 
lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
                 rep(0.05, sum(is.na(all_pars_lst[[2]]))),
                 rep(0.001, sum(is.na(all_pars_lst[[3]]))),
                 rep(0.001, sum(is.na(all_pars_lst[[4]]))),
                 rep(0.001, p))

# [1]    NA    NA    NA 0.050
# [5] 0.050 0.050 0.001 0.001
# [9] 0.001 0.001 0.001 0.001
# [13] 0.001 0.001 0.001



#================
# Optim: Tir-Wave
#================

optm_pars_Matern_TW <- optim(par = c(Vals, rep(1, p)), # ini guess
                             fn = neg_logL_Matern,
                             p = p, data_str = hierarchy_data6, 
                             all_pars_lst = all_pars_lst_Matern_6, df = df_TW, 
                             fit_indx = fit_indx,
                             method = "L-BFGS-B",
                             lower = lower_bound,
                             control = list(trace = 1, 
                                            maxit = 2,
                                            pgtol = 1e-5))



# final  value 1717.052764 
#stopped after 3 iterations
optm_pars_Matern_TW$message
# [1] "NEW_X"
optm_pars_Matern_TW$counts
# function gradient 
# 12       12 

# During a line search, the objective function may be evaluated multiple times 
# to find the step size that minimizes the objective function along the search direction

# it may evaluate the gradient multiple times within each iteration 
# to determine the search direction

# The optimization algorithm may perform additional function evaluations 
# or convergence criteria, such as changes in the objective function values 
# or parameter value

# Q: what's the difference with iteration?
# a simplified overview of what typically happens in an iteration 
# of an optimization algorithm

# 1. Current Solution
# 2. Compute Gradient: If the algorithm uses gradients (first-order derivatives)
# to guide the search, it computes the gradient of the objective function 
# with respect to the parameters at the current solution point. 
# This gradient provides information about the direction of steepest ascent 
# (for maximization) or descent (for minimization) of the objective function.


# 3. Update Parameters: Based on the gradient information the algorithm updates 
# the parameters (or decision variables) in the direction that is expected to improve the objective function value.


# 4. Convergence Check

# 5. Repeat or Terminate: If the convergence criteria are not met, 
# the algorithm repeats the process by starting a new iteration 
# with the updated parameters. Otherwise, if the convergence criteria 
# are satisfied, the algorithm terminates, and the current solution 
# is considered the optimal or near-optimal solution.


##================================
# Need to find a good inital start
##================================
# Grid search 


#================
# Optim: Wendland
#================

optm_pars_Matern_TW <- optim(par = c(Vals, rep(1, p)), # ini guess
                             fn = neg_logL_Matern,
                             p = p, data_str = hierarchy_data6, 
                             all_pars_lst = all_pars_lst_Matern_6, df = df_WL, 
                             fit_indx = fit_indx,
                             method = "L-BFGS-B",
                             lower = lower_bound,
                             control = list(trace = 1, 
                                            maxit = 1,
                                            pgtol = 1e-5))




