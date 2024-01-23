#============
# 24 Jan.2024
#============

# Aim:
  # modify 041_inference_functions_1D
  # to account for the Fit_indx selection to allow for
  # multi-tasks, within one function, 
  # denoising, denoising and filling the gap, C.V.


#==========
# neg_logL
#==========
# Arguments:  
# theta: a vector of all parameters to be inferenced
# p: the number of variates
# data_str: hierarchy data structure
# all_pars_lst: a list of different parameter matrices
# df: data.frame, s, Y1-Y3, Z1-Z3; see 043
# fit_indx: a vector of indices from list Fit_indx in 042

theta <- c(Vals, rep(1, p))
#rm(theta)

p = 3
data_str <- hierarchy_data_3

source("Fn_para_mat_construct.R")
all_pars_lst <- All_paras(p = p, data = data_str)

# df : see 043 generation true process and noisy data
fit_indx <- Fit_indx[[1]]



neg_logL <- function(theta, ..., p, data_str, all_pars_lst, df, fit_indx){

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
  
  Z_ft <- c(df_ft$Z1, df_ft$Z2, df_ft$Z3)
  #str(Z_ft) # num [1:450]
  
  
  # construct SIGMA_Y, SIGMA_Y_inv for process Y
  SG_SG_inv_Y_ft <- TST11_SpNormPert_SG_SGInv(p = p, data = data_str, 
                                           A_mat = all_pars_lst[[1]],
                                           dlt_mat = all_pars_lst[[2]],
                                           sig2_mat = all_pars_lst[[3]],
                                           kappa_mat = all_pars_lst[[4]], 
                                           d_vec = Dvec_ft, h = H_ft)
  
  
  SG_Y_ft <- SG_SG_inv_Y_ft$SIGMA
  SG_Y_inv_ft <- SG_SG_inv_Y_ft$SIGMA_inv
  str(SG_Y_ft) # num [1:450, 1:450]
  
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
  str(SG_Ng_ft) # num [1:450]
  
  
  # SIGMA, SIGMA_inv for observation Z
  SG_Z_ft = SG_Y_ft + SG_Ng_ft 
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng_ft <- SG_Y_inv_ft + SG_Ng_inv_ft
  SG_Y_Ng_inv_ft <- chol2inv(chol(SG_Y_Ng_ft))
  
  SG_Z_inv_ft <- SG_Ng_inv_ft - SG_Ng_inv_ft %*% SG_Y_Ng_inv_ft %*% SG_Ng_inv_ft
  
  
  # log_det(SG_Z)
  source("Fn_log_det.R")
  chol_SG_Z_ft <- chol(SG_Z_ft)
  log_SG_Z_det_ft <- log_det(chol_SG_Z_ft)
  str(log_SG_Z_det_ft) # a number 120
  
  # neg_logL
  L <- length(Z_ft) # different from n1 = length(Z1)
  neg_logL <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det_ft - 
                   1/2 * t(Z_ft) %*% SG_Z_inv_ft %*% Z_ft) # a 1 by 1 matrix
  
  neg_logL <- as.numeric(neg_logL) # a scalar
  
  
  # return scalar
  return(neg_logL)
  
}


#========
# Optim
#========

ini <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
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
                 rep(0.001, sum(is.na(all_pars_lst[[4]]))),
                 rep(0.001, p))

# [1]    NA    NA    NA 0.050
# [5] 0.050 0.050 0.001 0.001
# [9] 0.001 0.001 0.001 0.001
# [13] 0.001 0.001 0.001

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


# final  value 429.405998 
# converged

optm_pars
# [1] 1.0001286 0.9999801
# [3] 0.9998906 0.5484820
# [5] 0.3266648 0.6296211
# [7] 0.8184632 0.8053495
# [9] 0.9103205 1.5640553
# [11] 1.1578236 1.8021678
# [13] 0.2044648 0.2717899
# [15] 0.2973535

