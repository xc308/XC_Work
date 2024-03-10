#============
# 8 Mar. 2024
#============

# Aim:
  # 1D inf for 6 fields using Non cross-MRF, i.e., Tst9c Matern, 

# Method: 
  # SG_inv: using Tst9d, matern 
  source("Fn_TST9d_SpNReg_Thres_Matern_SG_SGInv.R")
  # df: df_TW, df_WL


#=========
# settings
#=========

#---------
# fit indx
#---------
# Ref: source("042_Tst_Fit_Obs_indx.R")
source("Fn_Tst_Fit_Obs_indx.R") # generate indx

Tst_Fit_Obs_indx <- Fn_Tst_Fit_Obs_indx(df = df_TW, num_folds = 4)

Fit_indx <- Tst_Fit_Obs_indx$Fit_indx # lst
fit_indx <- Fit_indx[[1]]


#-----
# pars
#-----
p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_Matern_6 <- All_paras(p = p, data = data_str)
all_pars_lst <- all_pars_lst_Matern_6 


# set ini
ini <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}

Vals

# together with tau2's 
theta <- c(Vals, rep(1, p)) # to get the neg_logL run


#---------------------------------------
# spatial grids, displacement, distance
#---------------------------------------
ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s)


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:200, 1:200]


# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:40000]


df <- df_TW
b = "Tri-Wave"

#==========
# neg_logL
#==========

neg_logL_Matern <- function(theta, ..., p, data_str, all_pars_lst, df, b = b,
                     fit_indx){
  
  source("Fn_TST9d_SpNReg_Thres_Matern_SG_SGInv.R")
  
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
  
  
  # select rows of df using Fit_indices to construct neg_logL for fit
  df_ft <- df[fit_indx, ]
  H_ft <- t(outer(df_ft$s, df_ft$s, FUN = "-"))
  #str(H_ft) # num [1:150, 1:150]
  Dvec_ft <- c(abs(H_ft))
  #str(Dvec_ft) # num [1:22500]
  
  # construct all obs Z under df_ft
  Z_ft <- c()
  for (i in 1:p) {
    Z_ft <- c(Z_ft, df_ft[[paste0("Z", i)]])
  }
  #str(Z_ft) # num [1:900]
  
  
  
  # construct SIGMA_Y, SIGMA_Y_inv for process Y
  SG_SG_inv_Y_ft <- TST9d_SpNormReg_SG_SGInv(p = p, data = data_str, 
                           A_mat = all_pars_lst[[1]],
                           dlt_mat = all_pars_lst[[2]],
                           sig2_mat = all_pars_lst[[3]],
                           kappa_mat = all_pars_lst[[4]],
                           d_vec = Dvec_ft, h = H_ft, b = b,
                           reg_ini = 1e-9,thres_ini = 1e-3)
  
  
  SG_Y_ft <- SG_SG_inv_Y_ft$SIGMA
  SG_Y_ft_inv <- SG_SG_inv_Y_ft$SIGMA_inv
  #str(SG_Y_ft) #  num [1:900, 1:900]
  
  # calculate SG_Ng
  ## 1st calcuate the # of parameters accumulated so far, 
  # so can connect theta components on top of current index
  # with measurement error tau2
  
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras(p = p, data = data_str)
  # for assign NA 
  
  SUM <- 0 
  for (i in 1:length(all_pars_lst)){
    s <- sum(is.na(all_pars_lst[[i]]))
    SUM <- SUM + s
  }
  
  
  # tau2 diag matrix
  #tau2_mat <- diag(theta[SUM+1], theta[SUM+2], ..., theta[SUM+p] )
  THETA <- c()
  for(i in 1:p){
    THETA <- c(THETA, theta[SUM + i])
  }
  tau2_mat <- diag(THETA)
  
  #n1 <- length(Z1) # total # of locations of univariate process
  n1 <- nrow(df_ft)
  I_sp_mat_ft <- I_sparse(size = n1, value = 1)
  SG_Ng_ft <- kronecker(tau2_mat, I_sp_mat_ft)
  SG_Ng_inv_ft <- solve(SG_Ng_ft)
  #str(SG_Ng_ft) # num [1:450]
  
  
  # SIGMA, SIGMA_inv for observation Z
  SG_Z_ft = SG_Y_ft + SG_Ng_ft 
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng_ft <- SG_Y_ft_inv + SG_Ng_inv_ft
  SG_Y_Ng_inv_ft <- chol2inv(chol(SG_Y_Ng_ft))
  
  SG_Z_inv_ft <- SG_Ng_inv_ft - SG_Ng_inv_ft %*% SG_Y_Ng_inv_ft %*% SG_Ng_inv_ft
  
  
  # log_det(SG_Z)
  source("Fn_log_det.R")
  chol_SG_Z_ft <- chol(SG_Z_ft)
  log_SG_Z_det_ft <- log_det(chol_SG_Z_ft)
  #str(log_SG_Z_det_ft) #  num 253
  
  # neg_logL
  L <- length(Z_ft) # different from n1 = length(Z1)
  neg_logL <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det_ft - 
                   1/2 * t(Z_ft) %*% SG_Z_inv_ft %*% Z_ft) # a 1 by 1 matrix
  
  neg_logL <- as.numeric(neg_logL) # a scalar [1] 1804.386
  
  
  # return scalar
  return(neg_logL)
  
}




