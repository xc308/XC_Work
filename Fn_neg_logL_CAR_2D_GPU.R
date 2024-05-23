#============
# 23 May 2024
#============
# Aim:
  # Function of neg_logL_CAR_2D_GPU


# Arg:
  # theta: a vector of initial parameters
  # p: the number of variates
  # data_str: graph structure among variate fields
  # all_pars_lst: list of parameters that reflect the graph structure
  # dsp_lon_mat: displacement Lon matrix
  # dsp_lat_mat: displacement Lat matrix
  # b: Tri-Wave or Wendland
  # phi: spatial decay parameter
  # H_adj: Adjacent matrix
  # df : data from containing observations Z


# Value:
  # a neg_logL scalar



neg_logL_CAR_2D_GPU <- function(theta, p, data_str, all_pars_lst, 
                                dsp_lon_mat, dsp_lat_mat, b = "Tri-Wave", phi, H_adj, df){
  
  source("Fn_TST12_SG_SGInv_CAR_2D_GPU.R")
  source("Fn_I_sparse.R")
  source("Fn_chol_inv_gpu.R") # input a gpumatrix, return chol inv
  
  
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
  SG_SG_inv_Y_gpu <- TST12_SG_SGInv_CAR_2D_GPU(p = p, data = data_str, 
                                               A_mat = all_pars_lst[[1]],
                                               dsp_lon_mat = dsp_lon_mat,
                                               dsp_lat_mat = dsp_lat_mat,
                                               dlt_lon_mat = all_pars_lst[[2]],
                                               dlt_lat_mat = all_pars_lst[[3]],
                                               b = b, phi = phi, 
                                               H_adj = H_adj, 
                                               sig2_mat = all_pars_lst[[4]],
                                               reg_ini = 1e-9, thres_ini = 1e-3)
  
  
  SG_Y_gpu <- SG_SG_inv_Y_gpu$SIGMA_gpu
  SG_Y_inv_gpu <- SG_SG_inv_Y_gpu$SIGMA_inv_gpu
  
  
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
  tau2_mat_gpu <- as.gpu.matrix(tau2_mat, device = "cuda")
  
  
  # total # of locations of univariate process
  n1 <- nrow(df) # 200
  I_sp_mat <- I_sparse(size = n1, value = 1)
  I_sp_gpu <- as.gpu.matrix(as.matrix(I_sp_mat), device = "cuda")
  
  SG_Ng_gpu <- tau2_mat_gpu %x% I_sp_gpu
  #SG_Ng <- kronecker(tau2_mat, I_sp_mat)
  
  SG_Ng_inv_gpu <- solve(SG_Ng_gpu)
  
  #str(SG_Ng_ft) # num [1:450]
  
  
  #  observation Z
  SG_Z_gpu = SG_Y_gpu + SG_Ng_gpu 
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng_gpu <- SG_Y_inv_gpu + SG_Ng_inv_gpu
  
  #SG_Y_Ng_inv <- chol2inv(chol(SG_Y_Ng))
  SG_Y_Ng_inv_gpu <- chol_inv_gpu(SG_Y_Ng_gpu)
  
  
  SG_Z_inv_gpu <- SG_Ng_inv_gpu - SG_Ng_inv_gpu %*% SG_Y_Ng_inv_gpu %*% SG_Ng_inv_gpu
  
  
  # log_det(SG_Z)
  #source("Fn_log_det.R")
  #chol_SG_Z <- chol(SG_Z)
  #log_SG_Z_det <- log_det(chol_SG_Z)
  log_SG_Z_det_gpu <- log(det(SG_Z_gpu))
  
  
  # construct joint Z, stack each Zi in df
  Z <- c()
  for (i in 1:p) {
    Z <- c(Z, df[[paste0("Z", i)]])
  }
  # str(Z)# num [1:1200] or use Fn_Stack_Z
  Z_gpu <- as.gpu.matrix(Z, device = "cuda")
  
  
  
  # neg_logL
  L <- length(Z_gpu) # different from n1 = length(Z1)
  neg_logL_gpu <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det_gpu - 
                       1/2 * t(Z_gpu) %*% SG_Z_inv_gpu %*% Z_gpu) # a 1 by 1 matrix on gpu
  
  neg_logL <- as.numeric(neg_logL_gpu, device = "cpu") # a scalar on cpu
  
  
  # return scalar
  return(neg_logL)
  
}

