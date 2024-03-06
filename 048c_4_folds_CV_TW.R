#=============
# 6 Mar. 2024
#=============
# Aim:
  # To do 4 folds C.V. for df_TW 

# Method:
  # df_TW: from 046b
  # neg_log_L for optim: 047b
  # cokrig formula: 048_cokrig_fomula
  # SG_inv generation: TST10c, cross-MRF


#=============
# pre-settings
#=============

#---------------------------------------------
# spatial grids, displacement, Adjacent matrix
#---------------------------------------------
ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

H <- t(outer(s, s, FUN = "-"))
str(H) # num [1:200, 1:200]

H_adj <- matrix(as.numeric(abs(H) < 0.4), nrow(H), nrow(H))
diag(H_adj) <- 0


eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1 / max(abs(eig)) # [1] 0.1412107
phi <- trunc(spec * 100)/100 # [1] 0.14


#------
# pars
#------
p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 


#-----
## ini
#-----
# A, del, sig2
ini <- c(0.3, 0.1, 1)
Vals <- c()
for(i in 1:length(all_pars_lst)){
  vals <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, vals)
}

# together with tau2's
all_ini_Vals <- c(Vals, rep(1, p))


#----------------------------------
# lower bound for A, dlt, sig2, tau2
#----------------------------------
# NA: no lower bound
# lower bd:
# A: NA
# dlt: 0.05
# sig2: 0.001
# tau2: 0.001

lower_bd <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
              rep(0.05, sum(is.na(all_pars_lst[[2]]))),
              rep(0.001, sum(is.na(all_pars_lst[[3]]))),
              rep(0.001, p))


#----------
# neg_logL
#----------

# see 047b neg_logL_CAR



#=======
# C.V.
#=======
# Ref:source("042_Tst_Fit_Obs_indx.R") # indx
source("Fn_Tst_Fit_Obs_indx.R") # generate indx
# cokrig: Ref "048_cokrig_formula_Tst10c_bchoice.R" 

#-----------------------------
# Generate Tst, Fit, Obs indx
#-----------------------------
Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)
Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 2)


Tst_indx <- Tst_Fit_Obs_indx$Tst_indx # lst
Fit_indx <- Tst_Fit_Obs_indx$Fit_indx
Obs_indx <- Tst_Fit_Obs_indx$Obs_indx


fit_indx <- Fit_indx[[1]]
tst_indx <- Tst_indx[[1]]
obs_indx <- Obs_indx[[1]]


#-------
# C.V.
#-------

num_folds <- 2
MAE <- c()
RMSE <- c()
for (f in 1:num_folds){
  
  fit_indx <- Fit_indx[[f]]
  obs_indx <- Obs_indx[[f]]
  tst_indx <- Tst_indx[[f]]
  
  ## neg_logL ##
  #neg_logL_CAR 047b
  
  ## optmization ##
  
  optm_pars <- optim(par = all_ini_Vals,
                     fn = neg_logL_CAR,
                     p = p, data_str= hierarchy_data6,
                     all_pars_lst = all_pars_lst_CAR_6, df = df_TW, 
                     fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
                     method = "L-BFGS-B",
                     lower = lower_bd,
                     control = list(trace = 1, 
                                    pgtol = 1e-5, 
                                    maxit = 200))$par
  
  
  
  ## cokrig##
  df_cokrig_TW <- cokrig(optm_pars_vec = optm_pars, all_pars_lst = all_pars_lst_CAR_6,
                         p = 6, data_str = hierarchy_data6, phi = phi, H_adj = H_adj, 
                         h = H, b = "Tri-Wave", obs_indx = obs_indx, df = df_TW)
  
  
  ## C.V. Metric calculation ##
  Y1_pred <- df_cokrig_TW$true_mu1[tst_indx]
  Y1_true <- df_cokrig_TW$smp_Y1[tst_indx]
  
  m <- length(Y1_pred)
  
  # Mean Absolute Error (MAE)
  MAE[f] <- sum(abs(Y1_pred - Y1_true))/m 
  RMSE[f] <- sqrt(sum((Y1_pred - Y1_true)^2) / m)
  
}












