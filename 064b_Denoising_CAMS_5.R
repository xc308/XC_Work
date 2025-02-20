#============
# 22 May 2025
#============

# Aim:
  # denoise the CAMS data

# Method:
  # 056

# Data: 
  # 064a: df_Lon_Strip_1_Sort_new.rds; optm_Lon_Strip_1_GPU_pars.rds


#=========
# Settings: grid locations, displacement, Adj matrix, phi range
#=========


df_Lon_Strp_1_Srt <- readRDS("df_Lon_Strip_1_Sort_new.rds")

head(df_Lon_Strp_1_Srt)

#----------
# 2D coords
#----------
crds <- cbind(df_Lon_Strp_1_Srt$Lon, df_Lon_Strp_1_Srt$Lat)
#head(crds)


#----------------
# 2D displacement
#----------------
# Aim: for TST12

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
str(DSP)
# num [1:3793, 1:3793, 1:2]

DSP[, , 1] # all Lon DSP
DSP[, , 2] # all Lat DSP


#--------------
# 2D H_adj, phi
#--------------
# Aim: for UniCAR in 2D

DIST <- as.matrix(dist(crds, diag = T, upper = T))
str(DIST) # num [1:3793, 1:3793]

Nb_radius <- 1.5 # degree 

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0

#dim(H_adj) # 3793 3793
#length(which(H_adj != 0)) # 27772

spec <- eigen(H_adj, symmetric = T, only.values = T)$val
#max(abs(spec)) # [1] 7.976106

phi <- 1/max(abs(spec)) # [1] 0.1253745
phi <- trunc(phi * 100)/100 # [1] 0.12


#---------
# data str
#---------
hierarchy_data_CAMS <- data.frame(
  node_id = c(1, 2, 3, 4,  5, 5),
  par_id = c(NA, 1, 2, 3, c(4, 1))
)

p = 5
data_str <- hierarchy_data_CAMS


#-----
# pars
#-----
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)
all_pars_lst <- all_pars_lst_CAR_2D_CMS



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
  
  
  # save true mu
  sub_vec <- length(all_true_mu) %/% p 
  each_true_mu_lst <- lapply(1:p, function(i) all_true_mu[((i-1)*sub_vec + 1) : (i*sub_vec)])
  
  for(i in 1:p){
    df[paste0("true_mu", i)] <- each_true_mu_lst[[i]]
  }
  
  # return df to preserve the above changes
  return(df)
}


#========
# cokrig
#========
# Aim: cokrig the true processes for each field using optimized parameters
# optm_Lon_Strip_1_GPU_pars <- readRDS("optm_Lon_Strip_1_GPU_pars.rds")

df_cokrig_2D_CAMS <- cokrig_denoise_2D(optm_pars_vec = optm_Lon_Strip_1_GPU_pars, 
                                     all_pars_lst = all_pars_lst_CAR_2D_CMS,
                                     p = p, data_str = hierarchy_data_CAMS,
                                     dsp_lon_mat = DSP[, , 1], 
                                     dsp_lat_mat = DSP[, , 2], 
                                     b = "Tri-Wave", 
                                     phi = phi, H_adj = H_adj, 
                                     df = df_Lon_Strp_1_Srt)



head(df_cokrig_2D_CAMS)
# 	Lon	Lat	Year	Z1	Z2	Z3	Z4	Z5	PM25_res	Lon_div_indx	true_mu1	true_mu2	true_mu3	true_mu4	true_mu5
#<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>	<dbl>
#  1	-179.25	66.75	2016	0.9560692	0.03698943	-0.14400463	-0.06502712	1.1542898	0.4596799	1	-0.7989188	-0.03317636	0.13852287	-0.01409495	-0.9705909
# 2	-179.25	67.50	2016	1.1864928	0.07368500	-0.04935164	0.06293557	0.9486452	0.4804245	1	-0.9914673	-0.06608915	0.04747299	0.01364159	-0.7976735

# 1. DU, 2. SU, 3. BC, 4. OM, 5. SS









