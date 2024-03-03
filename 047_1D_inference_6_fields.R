#=============
# 28 Jan. 2024
#=============

# Aim:
  # To go over the 1D inference using TST CAR for 6 fields
  
# Questions involved:
  # 1. since uni-CAR and matern has different parameters 
      # parameter list generate function All_pars() and 
      # the generated all_pars_lst still can be used in 
      # TST CAR? 

  # 2. if still can, then will the following steps involving
      # all_pars_lst be affected?

  # 3. more further questions: how shall we coordinate the 
      # the parameters in the list, the subscription of the 
      # position indicators within the list and the loop


# Method:
  # TST10_SpNormPert_SG_SGInv in 036_1D_simu_CAR_Chain_6F.R


#==========
# Settings
#==========
p = 6
data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)
all_pars_lst <- all_pars_lst_CAR_6 


#-----
# ini vals for run the neg_logL
#-----

ini <- c(1, 0.1, 1) # A, dlt, sig2, 
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}
Vals

# ini theta for assigning in the 1st part of neg_logL
theta <- c(Vals, rep(1, p))


#-----
# df and fit_indx
#-----
head(df) # in 046 generation of true process and noisy data
str(df)
Fit_indx
fit_indx <- Fit_indx[[1]]

chain = F



#=========
# neg_logL
#=========

neg_logL_CAR <- function(theta, ..., p, data_str, all_pars_lst, chain = F, 
                      df, fit_indx){
  
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
  #str(df_ft) # 150 obs 13 varialbes
  
  H_ft <- t(outer(df_ft$s, df_ft$s, FUN = "-"))
  #str(H_ft) # num [1:150, 1:150]
  
  # H_adj_ft
  H_adj_ft <- matrix(as.numeric(abs(H_ft) < 0.4), nrow(H_ft), nrow(H_ft))
  #str(H_adj_ft) # num [1:150, 1:150]
  diag(H_adj_ft) <- 0
  
  
  # phi and p.d. of (I - phi*H_adj_ft)
  eign_H_adj_ft <- eigen(H_adj_ft, symmetric = T, only.values = T)$val
  spec <- 1/ max(abs(eign_H_adj_ft)) # 0.1412107
  phi <- trunc(spec * 100) / 100     # 0.14
  
  
  Z_ft <- c(df_ft$Z1, df_ft$Z2, df_ft$Z3, df_ft$Z4, df_ft$Z5, df_ft$Z6)
  str(Z_ft) # num [1:900] 150 * 6
  
  
  # construct SIGMA_Y, SIGMA_Y_inv for process Y
  SG_SG_inv_Y_ft <- TST10_SpNormPert_SG_SGInv(p = p, data = data_str, 
                            chain = chain, A_mat = all_pars_lst[[1]], 
                            dlt_mat = all_pars_lst[[2]], 
                            sig2_mat = all_pars_lst[[3]], 
                            phi = phi, H_adj = H_adj_ft, h = H_ft)
  
  

  
  SG_Y_ft <- SG_SG_inv_Y_ft$SIGMA
  SG_Y_ft_inv <- SG_SG_inv_Y_ft$SIGMA_inv
  #str(SG_Y_ft) # num [1:900, 1:900] 150 * 6
  #str(SG_Y_ft_inv) # num [1:900, 1:900]
  
  # calculate SG_Ng
  ## 1st calcuate the # of parameters accumulated so far, 
  # so can connect theta components on top of current index
  # with measurement error tau2
  
  source("Fn_para_mat_construct.R")
  all_pars_lst <- All_paras_CAR(p = p, data = data_str)
  # for re-assign NA 
  
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
  n1 <- nrow(df_ft) # 150
  I_sp_mat_ft <- I_sparse(size = n1, value = 1)
  SG_Ng_ft <- kronecker(tau2_mat, I_sp_mat_ft)
  SG_Ng_ft_inv <- solve(SG_Ng_ft)
  #str(SG_Ng_ft) # num [1:900]
  
  
  # SIGMA, SIGMA_inv for observation Z
  SG_Z_ft = SG_Y_ft + SG_Ng_ft 
  #str(SG_Z_ft) # @ Dim     : int [1:2] 900 900
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng_ft <- SG_Y_ft_inv + SG_Ng_ft_inv
  SG_Y_Ng_ft_inv <- chol2inv(chol(SG_Y_Ng_ft))
  
  SG_Z_ft_inv <- SG_Ng_ft_inv - SG_Ng_ft_inv %*% SG_Y_Ng_ft_inv %*% SG_Ng_ft_inv
  
  
  # log_det(SG_Z)
  source("Fn_log_det.R")
  chol_SG_Z_ft <- chol(SG_Z_ft)
  log_SG_Z_det_ft <- log_det(chol_SG_Z_ft)
  #str(log_SG_Z_det_ft) # a number 1056
  
  # neg_logL
  L <- length(Z_ft) # 900 = 150*6
  neg_logL <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det_ft - 
                   1/2 * t(Z_ft) %*% SG_Z_ft_inv %*% Z_ft) # a 1 by 1 matrix
  
  neg_logL <- as.numeric(neg_logL) # a scalar
  
  
  # return scalar
  return(neg_logL)
  
}


#========
# Optim
#========

ini <- c(1, 0.1, 1) # A, dlt, sig2
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


## lower bound for each parameters, 
# NA: no lower bound
lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
                 rep(0.05, sum(is.na(all_pars_lst[[2]]))),
                 rep(0.001, sum(is.na(all_pars_lst[[3]]))),
                 rep(0.001, p))




optm_pars <- optim(par = c(Vals, rep(1, p)), # ini guess
                   fn = neg_logL_CAR,
                   p = p, data_str = hierarchy_data6, 
                   all_pars_lst = all_pars_lst_CAR_6, df = df, 
                   fit_indx = fit_indx,
                   method = "L-BFGS-B",
                   lower = lower_bound,
                   control = list(trace = 1, 
                                  pgtol = 1e-5,
                                  maxit = 3000))$par


optm_pars
# [1] 1.0000449 1.0002315
# [3] 0.9997245 1.0001348
# [5] 0.9998604 0.9999977
# [7] 1.0000037 0.9999943
# [9] 0.9999991 0.2334973
# [11] 0.2310299 0.2793296
# [13] 0.3908999 0.3932519
# [15] 0.5052313 0.6556379
# [17] 0.6127397 0.5988702
# [19] 0.9670511 1.1074911
# [21] 1.0904058 0.8978010
# [23] 1.2607710 1.1038912
# [25] 0.1491173 0.1270924
# [27] 0.2365520 0.4429848
# [29] 0.2290260 0.6449654






