#==============
# 14 Mar. 2024
#==============

# Aim:
  # write neg_logL function for 2D TST12 
  # without fit_index, due to application of pure denoising, no filling the gap


# Method:
  # TST: source("054_2D_simu_CAR_6_TST12.R") # TST
  # df: 046c_2D_true_process_generation


#=============
# pre-setting
#=============

data_str <- hierarchy_data6


#-----
# pars
#-----
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6_2D <- All_paras_CAR_2D(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6_2D



#-------
# set some ini for theta to run the neg_logL
#-------
Vals <- c(0.7, 0.8, 0.8, 0.65, 0.65, 0.5, rep(0.6, 3), 
               0.2, 0.3, 0.3, 0.4, 0.4, 0.5, rep(0.6, 3),
               rep(0.2, 6), rep(0.1, 3), 
               rep(1, 6)) # w/o tau2s


theta <- c(Vals, rep(1, p)) # with tau2s


dsp_lon_mat <- DSP[, , 1]
dsp_lat_mat <- DSP[, , 2]


# H_adj and phi from 046c
b <- "Tri-Wave"

df <- df_2D_TW


#=========
# neg_logL
#=========


neg_logL_CAR_2D <- function(theta, ..., p, data_str, all_pars_lst, 
                     dsp_lon_mat, dsp_lat_mat, b, phi, H_adj, df){
  
  source("Fn_TST12_SG_SGInv_CAR_2D.R")
  source("Fn_I_sparse.R")
  
  
  # connect each component of theta to all_pars_lst 
  # to incoporate each theta component into the neg log L function
  
  theta_indx <- 1
  for (lst in 1:length(all_pars_lst)){
    for (i in 1:nrow(all_pars_lst[[lst]])){
      for (j in 1:ncol(all_pars_lst[[lst]])){
        if (is.na(all_pars_lst[[lst]][i, j])){
          all_pars_lst[[lst]][i, j] <- theta[theta_indx]
          theta_indx <- theta_indx + 1
        } 
      }
    }
  }
  
  
  
  # construct SIGMA_Y, SIGMA_Y_inv for process Y
  SG_SG_inv_Y <- TST12_SG_SGInv_CAR_2D(p = p, data = data_str, 
                                       A_mat = all_pars_lst[[1]],
                                       dsp_lon_mat = dsp_lon_mat,
                                       dsp_lat_mat = dsp_lat_mat,
                                       dlt_lon_mat = all_pars_lst[[2]],
                                       dlt_lat_mat = all_pars_lst[[3]],
                                       b = b, phi = phi, 
                                       H_adj = H_adj, 
                                       sig2_mat = all_pars_lst[[4]],
                                       reg_ini = 1e-9, thres_ini = 1e-3)
    
  
  SG_Y <- SG_SG_inv_Y$SIGMA
  SG_Y_inv <- SG_SG_inv_Y$SIGMA_inv
  str(SG_Y) # num [1:1200, 1:1200]
  
  # calculate SG_Ng
  ## 1st calcuate the # of parameters accumulated so far, 
  # so can connect theta components on top of current index
  # with measurement error tau2
  
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras_CAR_2D(p = p, data = data_str)
  # for assign NA 
  
  SUM <- 0 
  for (i in 1:length(all_pars_lst)){
    s <- sum(is.na(all_pars_lst[[i]]))
    SUM <- SUM + s
  }
   
  # SUM 33
  
  # tau2 diag matrix
  #tau2_mat <- diag(theta[SUM+1], theta[SUM+2], ..., theta[SUM+p] )
  THETA <- c()
  for(i in 1:p){
    THETA <- c(THETA, theta[SUM + i])
  }
  tau2_mat <- diag(THETA)
  
  #n1 <- length(Z1) # total # of locations of univariate process
  n1 <- nrow(df) # 200
  I_sp_mat <- I_sparse(size = n1, value = 1)
  SG_Ng <- kronecker(tau2_mat, I_sp_mat)
  SG_Ng_inv <- solve(SG_Ng)
  #str(SG_Ng_ft) # num [1:450]
  
  
  # SIGMA, SIGMA_inv for observation Z
  SG_Z = SG_Y + SG_Ng 
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng <- SG_Y_inv + SG_Ng_inv
  SG_Y_Ng_inv <- chol2inv(chol(SG_Y_Ng))
  
  SG_Z_inv <- SG_Ng_inv - SG_Ng_inv %*% SG_Y_Ng_inv %*% SG_Ng_inv
  
  
  # log_det(SG_Z)
  source("Fn_log_det.R")
  chol_SG_Z <- chol(SG_Z)
  log_SG_Z_det <- log_det(chol_SG_Z)
  #str(log_SG_Z_det) # num 987
  
  
  # # construct joint Z, stack each Zi in df
  Z <- c()
  for (i in 1:p) {
    Z <- c(Z, df[[paste0("Z", i)]])
  }
  
  # str(Z)# num [1:1200] or use Fn_Stack_Z
  
  # neg_logL
  L <- length(Z) # different from n1 = length(Z1)
  neg_logL <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det - 
                   1/2 * t(Z) %*% SG_Z_inv %*% Z) # a 1 by 1 matrix
  
  neg_logL <- as.numeric(neg_logL) # a scalar
  
  
  # return scalar
  return(neg_logL)
  
}



















