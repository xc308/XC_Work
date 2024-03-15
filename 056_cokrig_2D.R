#=============
# 15 Mar. 2024
#=============

# Aim:
  # cokring all the true processes (pure denoising) for 2D inference

# Method
  # df_2D_TW: 046c; readRDS("df_2D_TW.rds")
  # optm_pars_2D_TW: 055; readRDS("optm_pars_2D_TW.rds")
  # cokrig formula: 048_cokrig_fomula


#=========
# Settings: grid locations, displacement, Adj matrix, phi range
#=========

#---------------
# 2D Grid locations
#---------------
# Aim: for df_2D

ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

coords <- cbind(s, s)
head(coords)
tail(coords)


#----------------
# 2D displacement
#----------------
# Aim: for TST12

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = coords)
str(DSP)
# num [1:200, 1:200, 1:2]

DSP[, , 1] # all Lon DSP
DSP[, , 2] # all Lat DSP


#--------------
# 2D H_adj, phi
#--------------
# Aim: for UniCAR in 2D

Dist <- as.matrix(dist(coords, diag = T, upper = T))
str(Dist) # num [1:200, 1:200]

#Nb_radius <- 0.8 # lag-5
Nb_radius <- 0.6 # lag-4
H_adj <- matrix(as.numeric(Dist < Nb_radius), nrow(Dist), nrow(Dist))
diag(H_adj) <- 0
H_adj


eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))

phi <- trunc(spec * 100)/100
# [1] 0.12


#==========================
# graph structure: 6 fields
#=========================

hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


data_str <- hierarchy_data6


#-----
# pars
#-----
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6_2D <- All_paras_CAR_2D(p = p, data = data_str)


#optm_pars_2D_TW[-(1:SUM)]


#================================
# Cokrig Function: pure denoising
#================================


cokrig_denoise_2D <- function(optm_pars_vec, all_pars_lst, p, data_str, 
                      phi, H_adj, dsp_lon_mat, dsp_lat_mat, b = "Tri-Wave", df){
  
  source("Fn_I_sparse.R")
  source("Fn_vec2mat.R")
  source("Fn_TST12_SG_SGInv_CAR_2D.R") # for TST12
  
  optm_pars_lst <- vec_2_mat(vector = optm_pars_vec, all_pars_lst = all_pars_lst)
  
  
  optm_SG_SG_inv <- TST12_SG_SGInv_CAR_2D(p = p, data = data_str, 
                                          A_mat = optm_pars_lst[[1]],
                                          dsp_lon_mat = dsp_lon_mat,
                                          dsp_lat_mat = dsp_lat_mat,
                                          dlt_lon_mat = optm_pars_lst[[2]],
                                          dlt_lat_mat = optm_pars_lst[[3]],
                                          b = b, phi = phi, H_adj = H_adj, 
                                          sig2_mat = optm_pars_lst[[4]],
                                          reg_ini = 1e-9, thres_ini = 1e-3)
  
  
  optm_SG_Y <- optm_SG_SG_inv$SIGMA
 
  
  # construct SG_Ng_obs
  ## 1st need to know how many non-tau2s 
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras_CAR_2D(p = p, data = data_str)
  # for re-assign NA 
  
  SUM <- 0 
  for (i in 1:length(all_pars_lst)){
    s <- sum(is.na(all_pars_lst[[i]]))
    SUM <- SUM + s
  } # 33
  
  tau2_optm <- optm_pars_vec[-(1:SUM)]
  tau2_optm_mat <- diag(tau2_optm)
  
  n1 <- nrow(df)
  I_mat <- I_sparse(size = n1, value = 1)
  SG_Ng <- kronecker(tau2_optm_mat, I_mat)
  #SG_Ng_obs <- SG_Ng[obs_indx, obs_indx]
  
  
  # construct SG_Z_obs
  SG_Z <- optm_SG_Y + SG_Ng
  
  
  # construct SG_Z_inv
  SG_Z_chol <- chol(SG_Z)
  SG_Z_inv <- chol2inv(SG_Z_chol)
  
  # construct joint Z, stack each Zi in df
  Z <- c()
  for (i in 1:p) {
    Z <- c(Z, df[[paste0("Z", i)]])
  }
  
  

  # cokring pure denoising
  all_true_mu <- optm_SG_Y %*% SG_Z_inv %*% Z
  #str(all_true_mu) # Formal class 'dgeMatrix', but can still sub-selection
  all_true_var <- diag(optm_SG_Y - optm_SG_Y %*% SG_Z_inv %*% optm_SG_Y)
  
  
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
# cokrig
#========
# Aim: cokrig the true processes for each field using optimized parameters

df_cokrig_2D_TW <- cokrig_denoise_2D(optm_pars_vec = optm_pars_2D_TW, 
                  all_pars_lst = all_pars_lst_CAR_6_2D,
                  p = p, data_str = hierarchy_data6,
                  dsp_lon_mat = DSP[, , 1], 
                  dsp_lat_mat = DSP[, , 2], 
                  b = "Tri-Wave", 
                  phi = phi, H_adj = H_adj, 
                  df = df_2D_TW)


head(df_cokrig_2D_TW)
#    lon   lat      smp_Y1      smp_Y2     smp_Y3   smp_Y4
# 1 -9.95 -9.95  0.10083728 -1.56349529  1.0610056 1.462496

#      smp_Y5   smp_Y6          Z1         Z2         Z3        Z4
#1 1.9154356 3.815431  0.60167016 -1.0989157  0.3484638 1.9253850

# Z5       Z6    true_mu1    true_mu2    true_mu3
# 1  1.2178417 4.104056  0.52795320 -0.87518829  0.33254777

#  true_mu4  true_mu5 true_mu6 true_var1 true_var2 true_var3
# 1 2.000137 1.3547079 3.759893 0.1018876 0.2427961 0.1426144

#   true_var4 true_var5 true_var6
# 1 0.3468417 0.3739804 0.3480970
  







