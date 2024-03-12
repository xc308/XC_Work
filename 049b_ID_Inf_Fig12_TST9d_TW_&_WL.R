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

ini <- c(0.5, 0.5, 1, 2) # A, dlt, sig2, kappa
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
                             fit_indx = fit_indx, b = "Tri-Wave",
                             method = "L-BFGS-B",
                             lower = lower_bound,
                             control = list(trace = 0, 
                                            maxit = 200,
                                            pgtol = 1e-4))


# 1st run
# final  value 1717.052764 
#stopped after 3 iterations
optm_pars_Matern_TW$message
# [1] "NEW_X"
optm_pars_Matern_TW$counts
# function gradient 
# 12       12 

# 2nd run 
#final  value 1563.972442 
#stopped after 85 iterations
optm_pars_Matern_TW$message
# [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"
# initial value

optm_pars_Matern_TW$counts
# function gradient 
# 256      256 

optm_pars_Matern_TW$par
# [1] 0.48628902 0.54297915
#[3] 0.55838092 0.53285438
#[5] 0.49051338 0.50782285
#[7] 0.50001489 0.49953224
#[9] 0.49947753 0.08658082
#[11] 0.42727072 2.22067814
#[13] 1.05183670 0.42125199
#[15] 0.47794978 0.19150125
#[17] 0.56442973 0.65038259
#[19] 1.63223177 1.15320121
#[21] 2.07270914 1.55084153
#[23] 1.40549876 0.67715591
#[25] 2.57817159 1.36318412
#[27] 2.96140610 2.47055638
#[29] 1.55639986 1.17298566
#[31] 1.24305039 1.42695151
#[33] 1.63231830 1.44185847
#[35] 1.40587505 1.58684220

optm_pars_Matern_TW$convergence
# [1] 52


# 3rd run:
  # use the returned pars as 3rd run initial value
  # meanwhile change pgtol = 1e-4
optm_pars_Matern_TW$message
#[1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"
optm_pars_Matern_TW$convergence 
#[1] 52

optm_pars_Matern_TW$counts
#function gradient 
# 103      103 

optm_pars_Matern_TW$value
# [1] 1648.297

optm_pars_Matern_TW$par

# [1] 0.4829771 0.4803281 0.5261407
#[4] 0.4673359 0.5177298 0.5000071
#[7] 0.4995936 0.5015581 0.4936938
#[10] 2.2549959 0.3976075 2.6997030
#[13] 0.5816490 0.6544571 0.4532069
#[16] 0.5492003 0.2490240 0.7504331
#[19] 1.3748299 1.1990764 1.2761730
#[22] 0.9253619 1.3559191 0.3015382
#[25] 1.3712188 1.2371381 1.2213733
#[28] 1.2601855 1.1369781 1.0350282
#[31] 1.3310886 1.2168976 1.6480171
#[34] 1.7890516 1.7061369 1.8597839


#============
# Conclusion
#============
# For non-cross MRF, Matern + DAG
  # it's very difficulty to get the converged optimization results
  # due to the kappa in Matern is highly unlikely to converge



## from this point of view, cross-MRF evaude the kappa
  # and has a set parameters that are much easier to fit.





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




