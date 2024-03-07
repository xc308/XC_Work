#=============
# 27 Feb.2024
#=============

# Aim:
  # 1D inference for 6 fields in Fig12 
  # neg_logL function including fit_indx to allow for denoising, 
      # filling the gap and C.V.


# Method:
  # fit_indx obtained from source("Fn_Tst_Fit_Obs_indx.R")

  # 1. vector theta to mat for Tst10c function
  # 2. df_ft, H_ft, H_adj_ft, phi
  # 3. Tst10c for SG_Y_ft, SG_Y_ft_inv
  # 4. SG_Ng_ft
  # 5. SG_Z_ft
  # 6. SG_Z_ft_inv
  # 7. neg_logL

# df_TW: obtained from 046b_generate 


#=========
# settings
#=========

#---------
# fit indx
#---------
# Ref: source("042_Tst_Fit_Obs_indx.R")
source("Fn_Tst_Fit_Obs_indx.R") # generate indx

Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

Fit_indx <- Tst_Fit_Obs_indx$Fit_indx # lst
fit_indx <- Fit_indx[[1]]


#-----
# pars
#-----
p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 


#------------------------------
# ini value to run the neg_logL
#------------------------------

# A, del, sig2
ini <- c(1, 0.1, 1)

# expand to a vector
Val <- c()
for (i in 1:length(ini)){
    vals <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
    Val <- c(Val, vals)
  }
  

Val
vals_ini_vec <- c(Val, rep(1, p)) # together with tau2
theta <- vals_ini_vec


#==================
# neg_logL function
#==================
# Ref: source("034c_1D_simu_CAR_SpNReg_Thres_b_choice.R")

neg_logL_CAR <- function(theta, ..., p, data_str, all_pars_lst, b = "Tri-Wave", 
         df, Nb_radius = 0.4, fit_indx){
  #print("neg_logL_CAR function called") # for parallel
  
  source("Fn_TST10c_SpNReg_Thres_SG_SGInv_b_choice.R") # for TST
  
  ## theta vector to matrix for TST10c
  theta_indx <- 1
  for (lst in 1:length(all_pars_lst)) {
    for (i in 1:nrow(all_pars_lst[[lst]])){
      for (j in 1:ncol(all_pars_lst[[lst]])){
        if (is.na(all_pars_lst[[lst]][i, j])){
          all_pars_lst[[lst]][i, j] <- theta[theta_indx]
          theta_indx <- theta_indx + 1
        }
      }
    }
  }
  
  
  ## df_ft, H_ft, H_adj_ft, phi
  df_ft <- df[fit_indx, ]
  H_ft <- t(outer(df_ft$s, df_ft$s, "-"))
  
  H_adj_ft <- matrix(as.numeric(abs(H_ft) < Nb_radius), nrow(H_ft), nrow(H_ft))
  diag(H_adj_ft) <- 0
  #str(H_adj_ft) # num [1:150, 1:150]
  
  eg_H_adj_ft <- eigen(H_adj_ft, symmetric = T, only.values = T)$val
  spec <- 1/max(abs(eg_H_adj_ft)) # 0.1412107
  phi <- trunc(spec * 100)/100 # 0.14
  
  Z_ft <- c(df_ft$Z1, df_ft$Z2, df_ft$Z3, df_ft$Z4, df_ft$Z5, df_ft$Z6)
  #str(Z_ft) # num [1:900]
  
  
  ## Construct SG_Y_ft, SG_Y_ft_inv
  SG_SG_inv_Y_ft <- TST10c_SpNReg_Thres_SG_SGInv(p = p, data = data_str, 
                               A_mat = all_pars_lst[[1]],
                               dlt_mat = all_pars_lst[[2]],
                               sig2_mat = all_pars_lst[[3]],
                               phi = phi, H_adj = H_adj_ft,
                               h = H_ft, b = b, 
                               reg_ini = 1e-9, thres_ini = 1e-3)
  
  
  
  SG_Y_ft <- SG_SG_inv_Y_ft$SIGMA
  SG_Y_ft_inv <- SG_SG_inv_Y_ft$SIGMA_inv
  
  ## Construct SG_Ng_ft
  # re-assign NA in all_pars_lst
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras_CAR(p = p, data = data_str)
   
  # count the NAs accumulated in A_mat, dlt_mat, sig2_mat
  SUM <- 0
  for (lst in 1:length(all_pars_lst)){
    s <- sum(is.na(all_pars_lst[[lst]]))
    SUM <- SUM + s
  }
  # SUM = 24
  
  THETA <- c()
  for (i in 1:p){
     THETA <- c(THETA, theta[SUM + i])
  }
  
  tau2_mat <- diag(THETA)
  
  n1 <- nrow(df_ft) # 150
  source("Fn_I_sparse.R")
  I_sp_ft <- I_sparse(size = n1, value = 1)
  
  SG_Ng_ft <- kronecker(tau2_mat, I_sp_ft)
  #str(SG_Ng_ft) # @ x: num [1:900]
  
  SG_Ng_ft_inv <- solve(SG_Ng_ft)
  
  
  
  ## SG_Z_ft
  SG_Z_ft <- SG_Y_ft + SG_Ng_ft
  
  ## SG_Z_ft_inv
  SG_Y_Ng_ft <- SG_Y_ft_inv + SG_Ng_ft_inv
  SG_Y_Ng_ft_inv <- chol2inv(chol(SG_Y_Ng_ft))
  
  SG_Z_ft_inv <- SG_Ng_ft_inv - SG_Ng_ft_inv %*% SG_Y_Ng_ft_inv %*% SG_Ng_ft_inv
  
  
  ## neg_logL
  L <- length(Z_ft) # 900 = 150*6
  
  source("Fn_log_det.R")
  SG_Z_ft_chol <- chol(SG_Z_ft)
  logdet_SG_Z_ft <- log_det(SG_Z_ft_chol)
  
  neg_logL <- - (- (L/2) * log(2 * pi) -  1/2 * logdet_SG_Z_ft - 
                   1/2 * t(Z_ft) %*% SG_Z_ft_inv %*% Z_ft)
  
  neg_logL <- as.numeric(neg_logL) # 1707.153
  
  
  ## return
  return(neg_logL)
}



#=============
# Optimization
#=============

#-----
## ini
#-----
# A, del, sig2
#ini <- c(1, 0.1, 1)
ini <- c(0.2, 0.1, 0.5)
Vals <- c()
for(i in 1:length(all_pars_lst)){
  vals <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, vals)
}


# together with tau2's
#all_ini_Vals <- c(Vals, rep(1, p))
all_ini_Vals <- c(Vals, rep(0.2, p))

#---------------
# lower boundary
#---------------

# lower bd:
  # A: NA
  # dlt: 0.05
  # sig2: 0.001
  # tau2: 0.001

lower_bd <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
  rep(0.05, sum(is.na(all_pars_lst[[2]]))),
  rep(0.001, sum(is.na(all_pars_lst[[3]]))),
  rep(0.001, p))




#------------------------
# Without parallelzation
#------------------------
optm_pars <- optim(par = all_ini_Vals,
      fn = neg_logL_CAR,
      p = p, data_str= hierarchy_data6,
      all_pars_lst = all_pars_lst_CAR_6, df = df_TW, 
      fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
      method = "L-BFGS-B",
      lower = lower_bd,
      control = list(trace = 1, 
                     pgtol = 1e-5, 
                     maxit = 1000))


## 1st run problem
  # 1. typo fit_indx
  # 2. dlt lower bound typo 0.005, should be 0.05, cause potential numerical issue

## 2nd run problem
  #final  value 1613.269346 
  #stopped after 29 iterations

# optm_pars
#[1] 1.00000000 1.00000000
#[3] 1.00000000 1.00000000
#[5] 1.00000000 1.00000000
#[7] 1.00000000 1.00000000
#[9] 1.00000000 0.66954274
#[11] 0.13565162 0.08873377
#[13] 0.19409196 0.05307354
#[15] 0.43634767 0.10775537
#[17] 0.10740839 0.80025650
#[19] 0.92498019 0.92742585
#[21] 0.93314032 1.14215374
#[23] 0.94766956 1.77502088
#[25] 0.93006914 0.93736088
#[27] 0.99467995 0.94774043
#[29] 0.94032351 0.94202619

# guess: the initial guess too close to true value? 
  # so quickly stopped?
optm_pars_2nd <- optm_pars # for cv



# 3rd run: adjust ini-guess to enlarge the difference
optm_pars_3rd <- optm_pars

#iter   50 value 1590.905970
#final  value 1590.905970 
#converged

optm_pars_3rd
#  [1] 0.2000000 0.2000000 0.2000000
#[4] 0.2000000 0.2000000 0.2000000
#[7] 0.2000000 0.2000000 0.2000000
#[10] 0.0500000 3.9109734 6.0352246
#[13] 0.4499258 0.1167661 0.4619853
#[16] 0.7618806 0.7896155 1.1601078
#[19] 0.6293863 0.8932284 1.1000549
#[22] 1.0595865 0.7254632 1.1834395
#[25] 0.8371254 0.5700147 1.0913365
#[28] 0.6053106 0.4848475 0.6098196


## 4th run: modify Fn_Waves.R to free A from being constant 1
  # to see if this could better the par optimization of A and dlt

# iter   50 value 1590.905970
#final  value 1590.905970 
#converged

optm_pars$message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
# means that the algorithm has found a solution that satisfies 
  # its convergence criteria, and further iterations are not necessary.

optm_pars_4th <- optm_pars$par
# same as the 3rd


optm_pars$counts

optm_pars$convergence
# indicates that the optimization process converged satisfactorily.


#-----------
# Conclusion
#-----------

# even free A from 1 in the Tri-Wave function, 
# the optimized A's still remain the same as the ini values
# while the optimization still converged,
# due to the optimized function is Gaussian, hence unlikely to be
# multi-modal such that optimization stuck in a local optima
# so the only possible reason is that the for the simulated noisy data
# A is a redudant parameters; so even A is set to initial, it does not
# have any impact on the convergence. 
















#=================
# Parallelization
#=================
install.packages("parallel")
library(parallel)

install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(doParallel)




# detect the number of cores available
num_cores <- detectCores() #10
num_cores <- 3

# Set up a parallel backend using the makeCluster() function, 
  # specifying the number of CPU cores:
cl <- makePSOCKcluster(num_cores)
cl

# Register parallel backend
registerDoParallel(cl)

# Export cluster objects for each cluster to work parallizely

## Note:
b = "Tri-Wave" 
Nb_radius = 0.4 
# has been hard coding in above TST10c function

clusterExport(cl, c('all_ini_Vals', 'neg_logL_CAR',
                    'p', 'hierarchy_data6', 
                    'all_pars_lst_CAR_6', 'df', 
                    'fit_indx','lower_bd'))


foreach(i = 1:3) %dopar% {
  print(paste("Starting on core", i))
  
  tryCatch({
    optm_res <- optim(par = all_ini_Vals,
                      fn = neg_logL_CAR,
                      p = p, data_str= hierarchy_data6,
                      all_pars_lst = all_pars_lst_CAR_6, df = df, 
                      fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
                      method = "L-BFGS-B",
                      lower = lower_bd,
                      control = list(trace = 1, 
                                     pgtol = 1e-5, 
                                     maxit = 100)) 
    
    print("Optimization completed on core ", i)
    return(optm_res)
  }, error = function(e) {
    print(paste("Error occurred on core", i, ": ", conditionMessage(e)))
  })
}

stopCluster(cl)




optm_pars <- parLapply(cl, 1:num_cores, function(i){
  
  print("Starting optimization on core ", i) 
  optm_res <- optim(par = all_ini_Vals,
                     fn = neg_logL_CAR,
                     p = p, data_str= hierarchy_data6,
                     all_pars_lst = all_pars_lst_CAR_6, df = df, 
                     fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
                     method = "L-BFGS-B",
                     lower = lower_bd,
                     control = list(trace = 1, 
                                    pgtol = 1e-5, 
                                    maxit = 100)) 
  
  print("Optimization completed on core ", i)
  return(optm_res)
})

stopCluster(cl)


## Problem
# there's no output when I run the parLapply code; 
# the fan of the laptop starts. 







#-------
# Test
#-------

# Working

# Register parallel backend
registerDoParallel(cl)

foreach(i = 1:3) %dopar% {
  print(paste("Starting on core", i))
  print("Hi")
  print(paste("Test completed on core", i))
}
stopCluster(cl)




## NOT Working
# parallelize optim function on each core
parLapply(cl, 1:3, function(i){
  print("Starting on core ", i)
  print("Hi")
  print("Test completed on core", i)
  
})
stopCluster(cl)




optm_pars <- NULL
tryCatch({
  optm_pars <- parLapply(cl, 1:num_cores, function(i){
    print("Starting optimization on core ", i) 
    optm_res <- optim(par = all_ini_Vals,
                      fn = neg_logL_CAR,
                      p = p, data_str= hierarchy_data6,
                      all_pars_lst = all_pars_lst_CAR_6, df = df, 
                      fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
                      method = "L-BFGS-B",
                      lower = lower_bd,
                      control = list(trace = 1, 
                                     pgtol = 1e-5, 
                                     maxit = 100)) 
    
    print("Optimization completed on core ", i)
    return(optm_res)
  })
}, error = function(e) {
  stop("Error occurred during parallel computation: ", conditionMessage(e))
})





