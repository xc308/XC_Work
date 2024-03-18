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


#theta <- c(Vals, rep(1, p)) # with tau2s


#dsp_lon_mat <- DSP[, , 1]
#dsp_lat_mat <- DSP[, , 2]


# H_adj and phi from 046c
#b <- "Tri-Wave"

#df <- df_2D_TW


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
  #str(SG_Y) # num [1:1200, 1:1200]
  
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



#========
# Optim
#========

ini <- c(0.2, 0.1, 0.1, 0.5) # A, dlt_lon, dlt_lat, sig2
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


all_ini_Vals <- c(Vals, rep(0.1, p)) # with tau2s


## lower bound for each parameters, 
# NA: no lower bound
lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),  # A
                 rep(0.05, sum(is.na(all_pars_lst[[2]]))), # dlt_lon
                 rep(0.05, sum(is.na(all_pars_lst[[3]]))), # dlt_lat
                 rep(0.001, sum(is.na(all_pars_lst[[4]]))), # sig2
                 rep(0.001, p)) # tau2



optm_pars_CAR_2D_TW <- optim(par = all_ini_Vals, # ini guess
                   fn = neg_logL_CAR_2D,
                   p = p, data_str = hierarchy_data6, 
                   all_pars_lst = all_pars_lst_CAR_6_2D, 
                   dsp_lon_mat = DSP[, , 1], 
                   dsp_lat_mat = DSP[, , 2], 
                   b = "Tri-Wave", 
                   phi = phi, H_adj = H_adj,
                   df = df_2D_TW,
                   method = "L-BFGS-B",
                   lower = lower_bound,
                   control = list(trace = 0, 
                                  maxit = 300,
                                  pgtol = 1e-4))



#str(DSP[, , 1])
# num [1:200, 1:200]



optm_pars_CAR_2D_WL <- optim(par = all_ini_Vals, # ini guess
                             fn = neg_logL_CAR_2D,
                             p = p, data_str = hierarchy_data6, 
                             all_pars_lst = all_pars_lst_CAR_6_2D, 
                             dsp_lon_mat = DSP[, , 1], 
                             dsp_lat_mat = DSP[, , 2], 
                             b = "Wendland", 
                             phi = phi, H_adj = H_adj,
                             df = df_2D_WL,
                             method = "L-BFGS-B",
                             lower = lower_bound,
                             control = list(trace = 0, 
                                            maxit = 300,
                                            pgtol = 1e-4))


# 1st run results:
optm_pars_CAR_2D_WL$message
# [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"

optm_pars_CAR_2D_WL$convergence
# [1] 52

optm_pars_CAR_2D_WL$counts
# function gradient 
# 167      167 


optm_pars_CAR_2D_WL$value 
# [1] 1892.098

optm_pars_CAR_2D_WL$par

# [1] 0.200000000 0.200000000
#[3] 0.200000001 0.200000000
#[5] 0.200000000 0.200000000
#[7] 0.200000000 0.200000000
#[9] 0.200000000 0.050000000
#[11] 0.814168316 0.547179513
#[13] 0.050000000 0.064216076
#[15] 0.716281711 0.324964630
#[17] 0.840933681 0.999264695
#[19] 0.050000000 0.814168316
#[21] 0.547179513 0.050000000
#[23] 0.064216076 0.716281711
#[25] 0.324964630 0.840933681
#[27] 0.999264695 1.120650622
#[29] 0.526115555 0.641199971
#[31] 0.647073385 0.587038382
#[33] 1.024402358 0.002906765
#[35] 0.738870645 0.500614411
#[37] 0.709955483 0.519411667
#[39] 0.174092419


## 2nd run:
# adjust initial values based on the 1st run results
# A = 0.2
# dlt_lon = 0.8
# dlt_lat = 0.5
# sig2 = 1
# tau2 = 0.5


ini <- c(0.2, 0.8, 0.5, 1) # A, dlt_lon, dlt_lat, sig2
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


all_ini_Vals <- c(Vals, rep(0.5, p)) # with tau2s


optm_pars_CAR_2D_WL$message
#[1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

optm_pars_CAR_2D_WL$convergence
# [1] 0

optm_pars_CAR_2D_WL$counts
# function gradient 
# 112      112 

optm_pars_CAR_2D_WL$value
# [1] 1893.772


optm_pars_CAR_2D_WL$par
# [1] 0.20000000 0.20000000 0.20000000
#[4] 0.20000000 0.20000000 0.20000000
#[7] 0.20000000 0.20000000 0.20000000
#[10] 0.64915413 0.94172711 0.82823741
#[13] 0.70351314 0.66072452 0.81023674
#[16] 0.74358124 0.92863861 0.79851217
#[19] 0.36762732 0.63843495 0.52852900
#[22] 0.41271692 0.37802870 0.53312822
#[25] 0.46296913 0.66423517 0.55150740
#[28] 1.00800442 0.66374255 0.63249144
#[31] 0.77712224 0.65978567 0.97633003
#[34] 0.09494215 0.65355428 0.49134997
#[37] 0.58066433 0.47573630 0.19098267


