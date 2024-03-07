#============
# 4 Mar. 2024
#============
# Aim:
  # use the optm_par from 047b to do cokrig and C.V.

# Method:
  # df_TW: from 046b
  # optm_par: from 047b
  # cokrig formula: 048_cokrig_fomula


#=============
# pre-settings
#=============

#---------------------------
# spatial grid, displacement
#---------------------------
ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

H <- t(outer(s, s, FUN = "-"))
str(H) # num [1:200, 1:200]

H_adj <- matrix(as.numeric(abs(H) < 0.4), nrow(H), nrow(H))
diag(H_adj) <- 0

eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1 / max(abs(eig)) # [1] 0.06666361
phi <- trunc(spec * 100)/100 # trunc only preserves integer part 6


#-----------------------------
# Generate Tst, Fit, Obs indx
#-----------------------------
# Ref: 042_Tst_Fit_Obs_indx.R"
source("Fn_Tst_Fit_Obs_indx.R") # generate indx

Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

Tst_indx <- Tst_Fit_Obs_indx$Tst_indx # lst
Fit_indx <- Tst_Fit_Obs_indx$Fit_indx
Obs_indx <- Tst_Fit_Obs_indx$Obs_indx


fit_indx <- Fit_indx[[1]]
tst_indx <- Tst_indx[[1]]
obs_indx <- Obs_indx[[1]]



#------
# pars
#------
p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 



#=======
# cokrig
#=======
source("048_cokrig_formula_Tst10c_bchoice.R")

df_cokrig_TW <- cokrig(optm_pars_vec = optm_pars_3rd, all_pars_lst = all_pars_lst_CAR_6,
       p = 6, data_str = hierarchy_data6, phi = phi, H_adj = H_adj, 
       h = H, b = "Tri-Wave", obs_indx = obs_indx, df = df_TW)


Y_pred_TW <- df_cokrig_TW$true_mu1[tst_indx]
Y_true_TW <- df_cokrig_TW$smp_Y1[tst_indx]

m <- length(Y_pred_TW)

#=======
# C.V.
#=======

## MAE
sum(abs(Y_pred_TW - Y_true_TW)) / m # [1] 1.105604

## RMSE
sqrt(sum((Y_pred_TW - Y_true_TW)^2) / m) # [1] 1.380682

