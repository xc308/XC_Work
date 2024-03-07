#============
# 7 Mar. 2024
#============

# Aim:
  # use obtained optimized pars in 047c_WL to do cokrig and 1 fold C.V.


# Method:
  # cokrig: function source("Fn_cokrig")
  # optm_pars:  optm_pars_WL_converged in 047c_WL
  
  # C.V. metric: Mean Absolute Error (MAE); Root Mean Square Error (RMSE)
  


#========
# cokrig
#========

source("Fn_cokrig.R")
df_cokrig <- cokrig(optm_pars_vec = optm_pars_WL_converged, 
                    all_pars_lst = all_pars_lst_CAR_6,
                    p = p, data_str = hierarchy_data6, 
                    phi = phi, H_adj = H_adj, h = H,
                    b = "Wendland", obs_indx = obs_indx, df = df_WL)

# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001     


#============
# 1 fold C.V.
#============

source("Fn_Tst_Fit_Obs_indx.R")
Tst_Fit_Obs_lst <- Fn_Tst_Fit_Obs_indx(df = df_WL, num_folds = 4)

Obs_indx <- Tst_Fit_Obs_lst$Obs_indx
obs_indx <- Obs_indx[[1]]


Tst_indx <- Tst_Fit_Obs_lst$Tst_indx
tst_indx <- Tst_indx[[1]]


Y1_pred <- df_cokrig$true_var1[tst_indx]
Y1_true <- df_cokrig$smp_Y1[tst_indx]

m <- length(Y1_pred)

## MAE
sum(abs(Y1_pred - Y1_true)) / m # [1] 1.002774

## RMSE
sqrt(sum((Y1_pred - Y1_true)^2) / m) # [1] 1.219671
 







