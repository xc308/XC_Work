#=============
# 4 Mar. 2024
#=============

# Aim:
  # 4 fold C.V. for Tri.Wave

# Method:
  # data: df_TW 046b
  # optim: 047b
  # cokrig: 048_cokrig_fomula; Fn_cokrig.R

df_TW <- readRDS("df_TW.rds")
#str(df_TW)
# 'data.frame':	200 obs. of  13 variables:


#==============
# pre-settings
#==============

hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)

data_str <- hierarchy_data6


source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 

# set ini
ini <- c(0.2, 0.1, 1) # A, dlt, sig2
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}

# together with tau2's
all_ini_Vals <- c(Vals, rep(1, p))


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



#==============================
# Generate Tst, Fit, Obs indx
#==============================
# Ref: 042_Tst_Fit_Obs_indx.R"
source("Fn_Tst_Fit_Obs_indx.R") # generate indx

Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

Tst_indx <- Tst_Fit_Obs_indx$Tst_indx # lst
Fit_indx <- Tst_Fit_Obs_indx$Fit_indx
Obs_indx <- Tst_Fit_Obs_indx$Obs_indx

#str(Tst_indx)
# List of 4
#$ : int [1:50] 1 2 3 4 5 6 7 8 9 10 ...
#$ : int [1:50] 51 52 53 54 55 56 57 58 59 60 ...
#$ : int [1:50] 101 102 103 104 105 106 107 108 109 110 ...
#$ : int [1:50] 151 152 153 154 155 156 157 158 159 160 ...


#str(Fit_indx)
#List of 4

#str(Obs_indx)
#List of 4


#======================
# neg_logL_CAR function
#======================
neg_logL_CAR <- function(theta, p, data_str, all_pars_lst, b = "Tri-Wave", 
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



#==========
# Settings
#==========

ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
#str(s) # num [1:200]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H) 
#str(H) # num [1:200, 1:200]

Nb_radius <- 0.4

H_adj <- matrix(as.numeric(abs(H) < Nb_radius), nrow(H), nrow(H))
diag(H_adj) <- 0
#H_adj
#str(H_adj)

eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))

phi <- trunc(spec * 100)/100


#========
# C.V.
#========
source("Fn_cokrig.R")


num_folds <- 4
MAE <- numeric(num_folds) # initialize a vector of length 4, w/o re-allocation memory each time as c()
RMSE <- numeric(num_folds)
for (f in 1:num_folds) {
  
  fit_indx <- Fit_indx[[f]]
  obs_indx <- Obs_indx[[f]]
  tst_indx <- Tst_indx[[f]]
  
  ## neg_logL ##
  #neg_logL
  
  ## optmization ##
  optm_pars_TW <- optim(par = all_ini_Vals,
                        fn = neg_logL_CAR,
                        p = p, data_str= hierarchy_data6,
                        all_pars_lst = all_pars_lst_CAR_6, df = df_TW, 
                        fit_indx = fit_indx, b = "Tri-Wave", Nb_radius = 0.4,
                        method = "L-BFGS-B",
                        lower = lower_bd,
                        control = list(pgtol = 1e-3, 
                                       maxit = 250))$par
  
  
  
  df_cokrig <- cokrig(optm_pars_vec = optm_pars_TW, all_pars_lst = all_pars_lst_CAR_6,
                      p = 6, data_str = hierarchy_data6, phi = phi, H_adj = H_adj, 
                      h = H, b = "Tri-Wave", obs_indx = obs_indx, df = df_TW)
  
  
  
  ## C.V. Metric calculation ##
  Y1_pred <- df_cokrig$true_mu1[tst_indx]
  Y1_true <- df_cokrig$smp_Y1[tst_indx]
  
  m <- length(Y1_pred)
  
  
  MAE[f] <- sum(abs(Y1_pred - Y1_true))/m 
  RMSE[f] <- sqrt(sum((Y1_pred - Y1_true)^2) / m)
  
}





