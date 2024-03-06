#============
# 3 Mar. 2024
#============

# Aim:
  # cokring formula 
  # can be used for both df_WL(df)047c optm_pars_WL and df_TW 047b

# 

#==========
# Settings
#==========

ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H) 
str(H) # num [1:200, 1:200]

Nb_radius <- 0.4

H_adj <- matrix(as.numeric(abs(H) < Nb_radius), nrow(H), nrow(H))
diag(H_adj) <- 0
H_adj
str(H_adj)

eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))

phi <- trunc(spec * 100)/100


#---------
# Obs_indx
#---------
# Ref: 042_Tst_Fit_Obs_indx.R"
source("Fn_Tst_Fit_Obs_indx.R") # generate indx

Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

Obs_indx <- Tst_Fit_Obs_indx$Obs_indx #lst
obs_indx <- Obs_indx[[1]]




#========
# Cokrig
#========
cokrig <- function(optm_pars_vec = optm_pars, all_pars_lst, p, data_str, 
                   phi, H_adj, h, b = "Wendland", obs_indx, df){
  
  source("Fn_I_sparse.R")
  source("Fn_vec2mat.R")
  source("Fn_TST10c_SpNReg_Thres_SG_SGInv_b_choice.R") # for TST
  
  optm_pars_lst <- vec_2_mat(vector = optm_pars_vec, all_pars_lst = all_pars_lst)
  
  
  optm_SG_SG_inv <- TST10c_SpNReg_Thres_SG_SGInv(p = p, data = data_str, 
                                                 A_mat = optm_pars_lst[[1]],
                                                 dlt_mat = optm_pars_lst[[2]],
                                                 sig2_mat = optm_pars_lst[[3]],
                                                 phi = phi, H_adj = H_adj,
                                                 h = H, b = b, 
                                                 reg_ini = 1e-9, thres_ini = 1e-3)
  
  
  optm_SG_Y <- optm_SG_SG_inv$SIGMA
  # str(optm_SG_Y) # num [1:600, 1:600]
  optm_SG_Y_obs <- optm_SG_Y[obs_indx, obs_indx]
  #str(optm_SG_Y_obs) # num [1:550, 1:550]
  
  # construct SG_Ng_obs
  tau2_optm <- optm_pars_vec[-(1:length(Vals))]
  tau2_optm_mat <- diag(tau2_optm)
  
  n1 <- nrow(df)
  I_mat <- I_sparse(size = n1, value = 1)
  SG_Ng <- kronecker(tau2_optm_mat, I_mat)
  SG_Ng_obs <- SG_Ng[obs_indx, obs_indx]
  
  
  # construct SG_Z_obs
  SG_Z_obs <- optm_SG_Y_obs + SG_Ng_obs
  
  
  # construct SG_Z_obs_inv
  SG_Z_obs_chol <- chol(SG_Z_obs)
  SG_Z_obs_inv <- chol2inv(SG_Z_obs_chol)
  
  # all joint observations Z
  Z <- c(df$Z1, df$Z2, df$Z3, df$Z4, df$Z5, df$Z6)
  
  # cokring
  all_true_mu <- optm_SG_Y[, obs_indx] %*% SG_Z_obs_inv %*% Z[obs_indx]
  #str(all_true_mu) # Formal class 'dgeMatrix', but can still sub-selection
  all_true_var <- diag(optm_SG_Y - optm_SG_Y[, obs_indx] %*% SG_Z_obs_inv %*% optm_SG_Y[obs_indx, ])
  
  
  # save true mu
  sub_vec <- length(all_true_mu) %/% p 
  each_true_mu_lst <- lapply(1:p, function(i) all_true_mu[((i-1)*sub_vec + 1) : (i*sub_vec)])
  
  for(i in 1:p){
    df[paste0("true_mu", i)] <- each_true_mu_lst[[i]]
  }
  
  # save true var
  sub_vec <- length(all_true_var) %/% p
  each_var_lst <- lapply(1:p, function(i) all_true_var[((i-1)*sub_vec + 1):(i*sub_vec)])
  
  for (i in 1:p){
    df[paste0("true_var", i)] <- each_var_lst[[i]]
  }
  
  
  # return df to preserve the above changes
  return(df)
}


#========
# cokring
#========
# try using Tri-Wave optimzed pars in 047b
#df_cokrig <- cokrig(optm_pars_vec = optm_pars_3rd, all_pars_lst = all_pars_lst_CAR_6,
#       p = 6, data_str = hierarchy_data6, phi = phi, H_adj = H_adj, 
#       h = H, b = "Tri-Wave", obs_indx = obs_indx, df = df)



#Y_pred_TW <- df_cokrig$true_mu1[Tst_indx[[1]]]
#Y_true_TW <- df_cokrig$smp_Y1[Tst_indx[[1]]]

#m <- length(Y_pred_TW)


#======
# C.V.
#======
# MAE
#sum(abs(Y_pred_TW - Y_true_TW)) / m

# RMSE
#sqrt(sum(abs(Y_pred_TW - Y_true_TW)^2) / m)



















