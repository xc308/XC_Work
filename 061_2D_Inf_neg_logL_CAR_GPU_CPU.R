#=============
# 18 Apr. 2024
#=============
# Aim:
  # GPU version of 055_2D_Inf_neg_logL_CAR


# Methods:
  # source("055_2D_Inf_neg_logL_CAR_2D.R")
  # GPUmatrix
  # df: 046d


#==============
# GPU settings
#==============
#------
# torch
#------
# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)


#-----------
# GPUmatrix
#------------

#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

system("nvidia-smi")

require(Matrix)


#============================
# multi-CPU core pre-settings
#============================
#install.packages("optimParallel")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(optimParallel)



#==============
# Pre-settings
#==============

#-----
# df
#-----

df_2D_TW_CAMS <- readRDS("df_2D_TW_CAMS.rds")
#head(df_2D_TW_CAMS)
#str(df_2D_TW_CAMS)
# data.frame':	200 obs. of  12 variables


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


#-------
# set some ini for theta to run the neg_logL
#-------
#Vals <- c(0.1, 0.2, 0.1, 0.1, 0.1,  
#          0.1, 0.1, 0.1, 0.1, 0.1,
#          rep(0.2, 5), 
#          rep(1, 5)) # w/o tau2s


#theta <- c(Vals, rep(0.2, p))


#----------
# 2D coords
#----------
ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
crds <- cbind(s, s)

# nrow(crds) [1] 200


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
# for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
#str(DSP[, , 1]) # num [1:200, 1:200]
#str(DSP[, , 2])


#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR

DIST <- as.matrix(dist(crds, diag = T, upper = T))

Nb_radius <- 0.6 # lag-4

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0

spec <- eigen(H_adj, symmetric = T, only.values = T)$val
phi <- 1/max(abs(spec)) # 0.1251127
phi <- trunc(phi * 100)/100 # [1] 0.12


#------
# try Nb_radius <- 8 50%
#------
# Expect: more data left in H_adj
        # more diverse the eigen vectors and eig value
        # smaller 1/max(abs(eig))
        # smaller phi

## see:Investigate_phi_Nb_radius.R




#==================
# neg_logL function
#==================


neg_logL_CAR_2D_GPU <- function(theta, p, data_str, all_pars_lst, 
                            dsp_lon_mat, dsp_lat_mat, b = "Tri-Wave", 
                            phi, H_adj, df){
  
  #source("Fn_TST12_SG_SGInv_CAR_2D_GPU.R")
  #source("Fn_I_sparse.R")
  #source("Fn_chol_inv_gpu.R") # input a gpumatrix, return chol inv
  
  
  
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
  
  #source("Fn_para_mat_construct.R")
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



#==================================
# Test neg_logL_CAR_2D_GPU function (success!)
#==================================
#neg_logL_CAR_2D_GPU(theta = theta, p = p, data_str = hierarchy_data_CAMS,
#                    all_pars_lst = all_pars_lst_CAR_2D_CMS, 
#                    dsp_lon_mat = DSP[, , 1], 
#                    dsp_lat_mat = DSP[, , 2], 
#                    b = "Tri-Wave", 
#                    phi = phi, H_adj = H_adj, 
#                    df = df_2D_TW_CAMS)

# [1] 1460.098



#========
# Optim
#========
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)
all_pars_lst <- all_pars_lst_CAR_2D_CMS

#-----------
# ini values
#-----------
ini <- c(0.2, 0.1, 0.1, 0.5) # A, dlt_lon, dlt_lat, sig2
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}


all_ini_Vals <- c(Vals, rep(0.1, p)) # with tau2s


#---------------------------------
## lower bound for each parameters
#---------------------------------
# NA: no lower bound
lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),  # A
                 rep(0.05, sum(is.na(all_pars_lst[[2]]))), # dlt_lon
                 rep(0.05, sum(is.na(all_pars_lst[[3]]))), # dlt_lat
                 rep(0.001, sum(is.na(all_pars_lst[[4]]))), # sig2
                 rep(0.001, p)) # tau2



#=============================
# multi-core parallel settings
#=============================

# Create a cluster
cl <- makeCluster(4)  

# Set the created cluster as the default for parallel computation
setDefaultCluster(cl = cl)  

# Specify the library path on each worker
clusterEvalQ(cl, {
  .libPaths("/bask/projects/v/vjgo8416-xchen")
})

# Load the packages on each worker
clusterEvalQ(cl, {
  library(torch)
  library(GPUmatrix)
  library(optimParallel)
})


# Load the core R package (no need to specify the path)
clusterEvalQ(cl, {
  library(Matrix)
})


# Source functions here to the current env
source("Fn_Check_par_node.R")
source("Fn_I_sparse.R")
source("Fn_chol_inv_gpu.R") # input a gpumatrix, return chol inv
source("Fn_shft_dist_mat.R") # construct shifted distance matrix using shft displacement for b function
source("Fn_Waves.R")
source("Fn_Wendland_32.R") # R = 0.5
source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
source("Fn_forceSym_GPU.R") # forceSym on GPU
source("Fn_Tst_sym_pd_GPU.R") # GPU version
source("Fn_Thres_tune_cov_GPU.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv on GPU

source("Fn_para_mat_construct.R")
source("Fn_TST12_SG_SGInv_CAR_2D_GPU.R")



# Save each function as a variable in the current Env
Check_par_node <- get("Check_par_node")
I_sparse <- get("I_sparse")
chol_inv_gpu <- get("chol_inv_gpu")
Shft_dst_mat <- get("Shft_dst_mat")
TriWave_2D <- get("TriWave_2D")
WendLd32_2D <- get("WendLd32_2D")
check_set_SpNorm_Reg <- get("check_set_SpNorm_Reg")
forceSym_gpu <- get("forceSym_gpu")
Tst_sym_pd_gpu <- get("Tst_sym_pd_gpu")
spress_cov <- get("spress_cov")
check_pd_gpu <- get("check_pd_gpu")
Thres_tune_cov_gpu <- get("Thres_tune_cov_gpu")

TST12_SG_SGInv_CAR_2D_GPU <- get("TST12_SG_SGInv_CAR_2D_GPU")

Para_A_mat <- get("Para_A_mat")
Para_Dlt_lon_mat <- get("Para_Dlt_lon_mat")
Para_Dlt_lat_mat <- get("Para_Dlt_lat_mat")
Para_sig2_mat <- get("Para_sig2_mat")
All_paras_CAR_2D <- get("All_paras_CAR_2D")



# Export the variables to the cluster workers
clusterExport(cl, c("Check_par_node", "I_sparse", "chol_inv_gpu",
                    "TriWave_2D", "WendLd32_2D", "Tst_sym_pd_gpu",
                    "check_set_SpNorm_Reg", "spress_cov", "check_pd_gpu",
                    "Thres_tune_cov_gpu", "Shft_dst_mat", "forceSym_gpu",
                    "TST12_SG_SGInv_CAR_2D_GPU", "All_paras_CAR_2D", "Para_A_mat",
                    "Para_Dlt_lon_mat", "Para_Dlt_lat_mat", "Para_sig2_mat"))




# export each variable name to each worker
clusterExport(cl, c("all_ini_Vals", "p", 
                    "hierarchy_data_CAMS", "all_pars_lst_CAR_2D_CMS", 
                    "DSP", "phi", "H_adj", "df_2D_TW_CAMS", "lower_bound"))  


#---------
# Tri-Wave 
#---------

optm_pars_CAR_2D_TW_GPU_CPU <- optimParallel(par = all_ini_Vals, # ini guess
                             fn = neg_logL_CAR_2D_GPU,
                             p = p, data_str = hierarchy_data_CAMS, 
                             all_pars_lst = all_pars_lst_CAR_2D_CMS, 
                             dsp_lon_mat = DSP[, , 1], 
                             dsp_lat_mat = DSP[, , 2], 
                             b = "Tri-Wave", 
                             phi = phi, H_adj = H_adj,
                             df = df_2D_TW_CAMS,
                             method = "L-BFGS-B",
                             lower = lower_bound,
                             control = list(maxit = 1,
                                            factr=.05/.Machine$double.eps))


stopCluster(cl)

optm_pars_CAR_2D_TW_GPU_CPU






