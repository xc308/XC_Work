#============
# 26 Jan.2024
#============

# Aim:
  # C.V. use Fit_indx generated in 042_Tst_Fit_Obs_indx.R

#==============
# pre-settings
#==============
p <- 3
data_str <- hierarchy_data_3

d_vec <- D_vec
str(D_vec) # num [1:40000]

h <- H
str(H) #num [1:200, 1:200] 

head(df)
#     s      smp_Y1      smp_Y2     smp_Y3         Z1          Z2          Z3
# 1 -9.95 0.182526845 -0.07904852  0.5509568  0.1154027 -0.05305937  0.25250491

str(Z) # num [1:600] see 043_generation_true_process_noisy_data


#======
# C.V.
#======
# see 042_Tst_Fit_Obs_indx.R for Fit_Tst_Obs indx list generation


num_folds <- 4
MAE <- c()
RMSE <- c()
for (f in 1:num_folds){
  
  fit_indx <- Fit_indx[[f]]
  obs_indx <- Obs_indx[[f]]
  tst_indx <- Tst_indx[[f]]
  
  ## neg_logL ##
  #neg_logL
  
  ## optmization ##
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras(p = p, data = data_str)
  
  # set ini
  ini <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
  Vals <- c()
  for (i in 1:length(all_pars_lst)){
    value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
    Vals <- c(Vals, value)
  }
  
  # lower bound for each parameters, 
  # NA: no lower bound
  lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
                   rep(0.05, sum(is.na(all_pars_lst[[2]]))),
                   rep(0.001, sum(is.na(all_pars_lst[[3]]))),
                   rep(0.001, sum(is.na(all_pars_lst[[4]]))),
                   rep(0.001, p))
  
  
  
  optm_pars <- optim(par = c(Vals, rep(1, p)), # ini guess
                     fn = neg_logL,
                     p = p, data_str = hierarchy_data_3, 
                     all_pars_lst = all_pars_lst, df = df, 
                     fit_indx = fit_indx,
                     method = "L-BFGS-B",
                     lower = lower_bound,
                     control = list(trace = 1, 
                                    pgtol = 1e-5,
                                    maxit = 3000))$par
  
  
  
  ## cokrig##
  df_cokrig <- cokrig(optm_pars_vec = optm_pars, all_pars_lst = all_pars_lst,
                      p = p, data_str = hierarchy_data_3, d_vec = D_vec, h = H, 
                      chain = F, obs_indx = obs_indx, df, Z)
  
  
  ## C.V. Metric calculation ##
  Y1_pred <- df_cokrig$true_mu1[tst_indx]
  Y1_true <- df_cokrig$smp_Y1[tst_indx]
  
  m <- length(Y1_pred)
  
  # Mean Absolute Error (MAE)
  MAE[f] <- sum(abs(Y1_pred - Y1_true))/m # [1] 0.2646896
  RMSE[f] <- sqrt(sum((Y1_pred - Y1_true)^2) / m)
  
}

# final  value 430.508893 
# converged

# MAE
# [1] 0.2646896 0.1625117
# [3] 0.2178182 0.2922185

# RMSE
# [1] 0.3055528 0.2123753
# [3] 0.2595639 0.3947463

mean(MAE) # [1] 0.2343095
mean(RMSE) # [1] 0.2930596




