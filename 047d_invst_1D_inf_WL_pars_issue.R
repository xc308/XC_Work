#=============
# 4 Mar. 2024
#=============
# Aim:
  # investigate the optimization parameters and numerical stability





# Define optimization function
optm_pars_WL <- function(all_ini_Vals, p, hierarchy_data6, all_pars_lst_CAR_6, df, fit_indx, lower_bd, upper_bd) {
  # Define a function to log parameter values at each iteration
  log_parameter_values <- function(theta) {
    cat("Parameter Values:", theta, "\n")
  }
  
  # Define the negative log-likelihood function
  neg_logL_CAR <- function(theta, ..., p, data_str, all_pars_lst, b = "Wendland", 
                   df, Nb_radius = 0.4, fit_indx) {
    #
    
    
    source("Fn_TST10c_SpNReg_Thres_SG_SGInv_b_choice.R") # for TST10c
    
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

    
    # Call the logging function to print parameter values
    log_parameter_values(theta)
    
    # Return the negative log-likelihood value
    return(neg_logL)
  }
  
  # Perform optimization using optim()
  optimization_result <- optim(par = all_ini_Vals,
                               fn = neg_logL_CAR,
                               p = p,
                               data_str = hierarchy_data6,
                               all_pars_lst = all_pars_lst_CAR_6,
                               df = df,
                               fit_indx = fit_indx,
                               method = "L-BFGS-B",
                               lower = lower_bd,
                               upper = upper_bd,
                               control = list(trace = T, 
                                              pgtol = 1e-5, 
                                              maxit = 500))
  
  # Return the optimization result
  return(optimization_result)
}

# Call your optimization function with initial parameters
result <- optm_pars_WL(all_ini_Vals, p, hierarchy_data6, all_pars_lst_CAR_6, df, fit_indx, lower_bd, upper_bd)
