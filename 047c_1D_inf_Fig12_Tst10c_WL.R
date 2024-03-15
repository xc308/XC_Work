#============
# 3 Mar. 2024
#============

# Aim:
  # use true process and noisy data generated in 046b df_WL 
  # True parameters see 046b
  # to do the optimization use neg_logL in 047b

# Method:
  # neg_logL uses cross-MRF Tst10c with b choice


#===========
# settings
#===========

# Ref:"042_Tst_Fit_Obs_indx.R"
source("Fn_Tst_Fit_Obs_indx.R")
Tst_Fit_Obs_lst <- Fn_Tst_Fit_Obs_indx(df = df_WL, num_folds = 4)
Fit_indx <- Tst_Fit_Obs_lst$Fit_indx
fit_indx <- Fit_indx[[1]]


p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 


#==================
# neg_logL function
#==================

neg_logL_CAR <- function(theta, ..., p, data_str, all_pars_lst, b = "Wendland", 
                         df, Nb_radius = 0.4, fit_indx){
  
  source("Fn_TST10c_SpNReg_Thres_SG_SGInv_b_choice.R") # for TST10c
  
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
ini <- c(0.3, 0.1, 1)
Vals <- c()
for(i in 1:length(all_pars_lst)){
  vals <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, vals)
}

# together with tau2's
all_ini_Vals <- c(Vals, rep(1, p))


#-----------
# boundary
#-----------

# lower bd:
# A: NA
# dlt: 0.05
# sig2: 0.001
# tau2: 0.001

lower_bd <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
              rep(0.05, sum(is.na(all_pars_lst[[2]]))),
              rep(0.001, sum(is.na(all_pars_lst[[3]]))),
              rep(0.001, p))



# upper bd:
# A: NA
# dlt: 20, i.e., dlt <= max(abs(H)); reason see 034d investigate
# sig2: NA
# tau2: NA
#upper_bd <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
#              rep(20, sum(is.na(all_pars_lst[[2]]))),
#              rep(NA, sum(is.na(all_pars_lst[[3]]))),
#              rep(NA, p))




#------
# optm:use df_WL
#------
optm_pars_WL <- optim(par = all_ini_Vals,
                   fn = neg_logL_CAR,
                   p = p, data_str= hierarchy_data6,
                   all_pars_lst = all_pars_lst_CAR_6, df = df_WL, 
                   fit_indx = fit_indx, b = "Wendland", Nb_radius = 0.4,
                   method = "L-BFGS-B",
                   lower = lower_bd,
                   control = list(trace = 1, 
                                  maxit = 500,
                                  pgtol = 1e-5))

optm_pars_WL$message
# [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"
# Indicate poor ini values
optm_pars_WL$convergence
# [1] 52
# a convergence code of 52 suggests that the optimization process has converged successfully according to the specified convergence criterion.
optm_pars_WL_1st <- optm_pars_WL$par
head(df_WL)



#==========================
# To detect early poor ini
#==========================

custom_optim <- function(par, fn, ..., stop_on_line_search_failure = TRUE) {
  result <- optim(par = par, fn = fn, ...)
  
  # Check for line search failure
  if (stop_on_line_search_failure && grepl("ERROR: ABNORMAL_TERMINATION_IN_LNSRCH", result$message)) {
    warning("Line search failure detected. Terminating optimization.")
    return(result)
  }
  
  return(result)
}


optm_pars_WL <- custom_optim(par = all_ini_Vals,
                             fn = neg_logL_CAR,
                             p = p, 
                             data_str = hierarchy_data6,
                             all_pars_lst = all_pars_lst_CAR_6, 
                             df = df_WL, 
                             fit_indx = fit_indx, 
                             b = "Wendland", 
                             Nb_radius = 0.4,
                             method = "L-BFGS-B",
                             lower = lower_bd,
                             control = list(trace = 1, 
                                            pgtol = 1e-5, 
                                            maxit = 500))


# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#final  value 1469.260325 
#converged


optm_pars_WL_converged <- optm_pars_WL$par
optm_pars_WL$convergence #[1] 0
optm_pars_WL$message
optm_pars_WL$counts 
# function gradient 
#
#108      108 


optm_pars_WL_converged
# [1] 0.2999997 0.2999999 0.2999999
#[4] 0.3000001 0.3000001 0.3000000
#[7] 0.3000000 0.3000000 0.3000000
#[10] 0.1908242 0.3056953 0.4228101
#[13] 0.3731729 0.4684762 0.4904502
#[16] 0.1115371 0.5658871 0.5311577
#[19] 0.6307438 0.8018036 0.4783926
#[22] 0.8275809 0.6521290 0.7524731
#[25] 0.4679805 0.6683368 0.6326358
#[28] 0.6003878 0.5977174 0.6497580



#========================
# Try run optim only once
#========================

optm_pars_WL <- optim(par = all_ini_Vals,
                      fn = neg_logL_CAR,
                      p = p, data_str= hierarchy_data6,
                      all_pars_lst = all_pars_lst_CAR_6, df = df_WL, 
                      fit_indx = fit_indx, b = "Wendland", Nb_radius = 0.4,
                      method = "L-BFGS-B",
                      lower = lower_bd,
                      control = list(trace = 1, 
                                     maxit = 1,
                                     pgtol = 1e-5))




