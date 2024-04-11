#============
# 8 Apr. 2024
#============

# Aim:
  # GPU version of TST12 function


# 

TST12_SG_SGInv_CAR_2D_GPU <- function(p, data, A_mat, dsp_lon_mat, dsp_lat_mat, 
                                  dlt_lon_mat, dlt_lat_mat, b = "Wendland",
                                  phi, H_adj, sig2_mat,
                                  reg_ini = 1e-9, thres_ini = 1e-3) {
  
  #source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Tst_sym_pd_GPU.R")
  source("Fn_check_set_SpNorm_Reg_GPU.R") # SpN + tune regularize number
  source("Fn_I_sparse.R")
  source("Fn_Thres_tune_cov_GPU.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  source("Fn_shft_dist_mat.R") # construct shifted distance matrix using shft displacement for b function
  source("Fn_chol_inv_gpu.R") # chol inv for gpu matrix
  
  # phi: controls the propostion of conditional spatial dependence 
  # sig2_mat: diag matrix containing sigma2 in CAR
  # H: adjacency matrix
  
  #I_sps <- I_sparse(size = nrow(H_adj), value = 1)
  I_sps <- diag(1, nrow(H_adj), nrow(H_adj))
  I_sps_gpu <- as.gpu.matrix(I_sps)
  c_inv <- I_sps_gpu - phi * H_adj  # gpu
  #C11_inv <- c_inv %*% as.gpu.matrix(I_sparse(size = nrow(H_adj), value = 1/sig2_mat[1, 1]))
  C11_inv <- c_inv %*% as.gpu.matrix(diag(1/sig2_mat[1, 1], nrow(H_adj), nrow(H_adj)))
  
  
  #C11 <- chol2inv(chol(C11_inv)) 
  
  C11 <- chol_inv_gpu(C11_inv) # gpu
  
  n <- nrow(C11)
  SIGMA <- C11 # gpu
  
  reg_num <- reg_ini
  
  restart <- T
  while(restart){
    restart <- F
    SIGMA_inv_pd <- T # initializa label
    
    # reset the SIGMA each while iteration
    SIGMA <- C11
    
    for(r in seq(2, p, by = 1)){
      
      PN <- Check_par_node(Node = r, data = data)
      R <- C <- NULL
      
      cat("r:", r, "\n")
      for(c in seq(1, (r-1), by = 1)){
        
        BT <- NULL
        C_rc <- 0
        for(t in c(PN)){
          
          if (b == "Tri-Wave") {
            
            shft_dst_mat <- Shft_dst_mat(dsp_lon_mat = dsp_lon_mat, dsp_lat_mat = dsp_lat_mat, 
                                         dlt1 = dlt_lon_mat[r, t], dlt2 =  dlt_lat_mat[r, t])
            
            B_rt <- TriWave_2D(shft_dst_mat = shft_dst_mat, A = A_mat[r, t])
            B_rt <- as.gpu.matrix(B_rt)
            
            #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
          }
          
          if (b == "Wendland") {
            
            shft_dst_mat <- Shft_dst_mat(dsp_lon_mat = dsp_lon_mat, dsp_lat_mat = dsp_lat_mat, 
                                         dlt1 = dlt_lon_mat[r, t], dlt2 =  dlt_lat_mat[r, t])
            
            B_rt <- WendLd32_2D(shft_dst_mat = shft_dst_mat, A = A_mat[r, t])
            B_rt <- as.gpu.matrix(B_rt)
            #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
          }
          
          ## spectral normalization of B_rt
          B_rt <- check_set_SpNorm_Reg_gpu(B_rt, reg_num = reg_num) # gpu
          #cat("B cond numb:", kappa(B_rt), "\n")
          
          B_rt <- as.matrix(B_rt) # matrix for rbind with NULL
          
          BT <- rbind(BT, t(B_rt)) # matrix
          C_rc <- C_rc + B_rt %*% SIGMA[((t-1)*n+1) : (t*n), ((c-1)*n+1): (c*n)] # gpu
          
          C_rc <- as.matrix(C_rc)
        }
        
        R <- cbind(R, C_rc)
        
        C_cr <- t(C_rc)
        C <- rbind(C, C_cr)
      }
      
      #Drr_inv <- c_inv %*% as.gpu.matrix(I_sparse(size = nrow(H_adj), value = 1/sig2_mat[r, r])) # gpu
      Drr_inv <- c_inv %*% as.gpu.matrix(diag(1/sig2_mat[r, r], nrow(H_adj), nrow(H_adj)))
      
      #D_rr <- chol2inv(chol(Drr_inv))
      D_rr <- chol_inv_gpu(Drr_inv) # gpu
      
      
      t <- c(PN)
      Subset_cols <- function(t) {
        start_col <- (t - 1) * n + 1
        end_col <- t * n
        
        result <- R[, start_col:end_col]
      }
      
      result_lst <- lapply(t, FUN = Subset_cols)
      R_subset <- do.call(cbind, result_lst)
      
      C_rr <- R_subset %*% BT + D_rr # gpu
      
      
      SG <- SIGMA # p-1 for inverse, gpu
      SG <- forceSymmetric(as.matrix(SG)) # matrix
      
      Col <- rbind(C, C_rr) # gpu, C_rr is gpu
      Row <- rbind(SG, R)
      SIGMA <- cbind(Row, Col) # gpu
      SIGMA <- forceSymmetric(as.matrix(SIGMA)) # matrix
      
      
      ## SIGMA_Inv
      #Drr_inv <- chol2inv(chol(D_rr)) # Schur complement inverse
      
      if (r == 2){
        SG_inv <- C11_inv # gpu
      }
      
      CDrr_in <- C %*% Drr_inv     # gpu
      #CDR_sym <- forceSymmetric(CDrr_in %*% R)
      
      #CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
      cat("condition number of C", kappa(C), "\n")
      #cat("condition number of CDinv", kappa(CDrr_in), "\n")
      #cat("condition number of CDinvR", kappa(CDR_sym), "\n")
      
      SGCD <- SG_inv %*% CDrr_in # gpu
      RSG <- R %*% SG_inv # gpu
      BK1 <- SG_inv + SGCD %*% RSG   # gpu
      BK1 <- forceSymmetric(as.matrix(BK1)) # dsy matrix
      
      #BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
      #BK1 <- forceSymmetric(BK1)
      BK2 <- - SG_inv %*% (C %*% Drr_inv)   # gpu
      BK2 <- as.matrix(BK2)
      BK3 <- - (Drr_inv %*% R) %*% SG_inv   # gpu
      BK3 <- as.matrix(BK3)
      BK4 <- forceSymmetric(as.matrix(Drr_inv)) # dsy matrix
      
      #cat("BK4", "\n")
      #Tst_sym_pd(BK4)
      
      
      SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4)) # matrix & dsymatrix
      #SG_inv <- SIGMA_inv
      SG_inv <- forceSymmetric(SIGMA_inv) # dsymatrix
      SG_inv <- as.gpu.matrix(as.matrix(SG_inv)) # gpu 
      
      
      # to condition to break the for loop and start all over again
      if (!as.logical(check_pd_gpu(SG_inv))) {
        restart <- T
        SIGMA_inv_pd <- F
        cat("New reg_num needed.", "\n")
        
        if (reg_num < 1) {
          reg_num <- reg_num * 10
          cat("Reg_num updated to:", reg_num, "\n")
        } else {
          reg_num <- reg_num + 0.1
          cat("Reg_num updated to:", reg_num, "\n")
        }
        
        break # break the following of for loop and start from begining again
      }
      
      # early perturb SG_inv if not p.d.
      cat("r", r, "\n")
      cat("SG_inv", "\n")
      Tst_sym_pd_gpu(SG_inv) # gpu
    }
    
    cat("Final reg_num:", reg_num, "\n")
    
    if (SIGMA_inv_pd){
      
      # Compute SIGMA_inv_ini, SG_inv constructed with the smallest possible reg_num
      SIGMA_inv_ini <- SG_inv * (abs(SG_inv) > thres_ini) # gpu
      
      # 1. tune threshold if SIGMA_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIGMA_inv and SIGMA
      SG_SG_inv_thres <- Thres_tune_cov_gpu(thres_ini = thres_ini, 
                                        cov_mat_thres = SIGMA_inv_ini,
                                        cov_mat = SG_inv) 
      
      return(list(SIGMA = as.matrix(SIGMA),
                  SIGMA_inv = SG_SG_inv_thres$SIGMA_inv))  
      
    }
  }
}


#===========================
# Test on 2D simulated data
#===========================

#----------
# 2D coords
#----------
ds <- 0.1
#s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)


crds <- cbind(s, s)
head(crds)
nrow(crds) # [1] 20. [1] 200
str(crds) # num [1:20, 1:2]# num [1:200, 1:2] 


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
  # for TST12_gpu function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
str(DSP[, , 1]) # all Lon displacement
# num [1:20, 1:20]; num [1:200, 1:200]
str(DSP[, , 2]) # all Lat displacement
# num [1:20, 1:20] ; num [1:200, 1:200]



#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR


#crds_gpu <- as.gpu.matrix(crds)
#crds_gpu 
# [ CPUDoubleType{20,2} ]



DIST <- as.matrix(dist(crds, diag = T, upper = T))
str(DIST)
# num [1:20, 1:20]; num [1:200, 1:200]

quantile(DIST)
#       0%       25%       50%       75%      100% 
# 0.0000000 0.4242641 0.8485281 1.4142136 2.6870058 

#       0%       25%       50%       75%      100% 
# 0.000000  3.818377  8.343860 14.142136 28.142850 



Nb_radius <- 0.8 # 50% of DIST matrix will be set to zero
Nb_radius <- 8

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0
str(H_adj) #  num [1:200, 1:200]

eig_val <- eigen(as.gpu.matrix(H_adj), symmetric = T, only.values = T)$val
spec <- eig_val@gm$real


phi <- 1/max(abs(spec)) # [1] 0.1344431; 0.0098534

phi <- trunc(phi * 100)/100  # [1] 0.13

phi <- trunc(phi * 1000)/1000  # [1] 0.009


#--------
# pars
#--------

hierarchy_data_CAMS <- data.frame(
  node_id = c(1, 2, 3, 4,  5, 5),
  par_id = c(NA, 1, 2, 3, c(4, 1))
)

p = 5

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)

source("Fn_set_ini_vals.R")
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[1]], ini_vals = 1)
dlt_lon_02 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[2]], ini_vals = 0.2)
dlt_lat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[3]], ini_vals = 0.4)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[4]], ini_vals = 1)



## Tri-Wave
SG_SGinv_CAR_6_2D_TW_GPU <- TST12_SG_SGInv_CAR_2D_GPU(p = 5, data = hierarchy_data_CAMS, 
                                              A_mat = A_1, 
                                              dsp_lon_mat = DSP[, , 1], 
                                              dsp_lat_mat = DSP[, , 2],
                                              dlt_lon_mat = dlt_lon_02, 
                                              dlt_lat_mat = dlt_lat_04, 
                                              b = "Tri-Wave", 
                                              phi =  phi, H_adj = H_adj,
                                              sig2_mat = sig2_mat_1, 
                                              reg_ini = 1e-9, thres_ini = 1e-3)


#---------------
# 200 locations
#---------------

# r: 2 
#condition number of C 2.312059e+18 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 2.395258e+16 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 2.443129e+16 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 2.5076e+16 
#r 5 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


#--------------
# 20 locations
#--------------
# r: 2 
#condition number of C 108750.4 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 179037.2 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 226578.7 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 195287.4 
#r 5 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 





SG_SGinv_CAR_6_2D_WL_GPU <- TST12_SG_SGInv_CAR_2D_GPU(p = 5, data = hierarchy_data_CAMS, 
                                                      A_mat = A_1, 
                                                      dsp_lon_mat = DSP[, , 1], 
                                                      dsp_lat_mat = DSP[, , 2],
                                                      dlt_lon_mat = dlt_lon_02, 
                                                      dlt_lat_mat = dlt_lat_04, 
                                                      b = "Wendland", 
                                                      phi =  phi, H_adj = H_adj,
                                                      sig2_mat = sig2_mat_1, 
                                                      reg_ini = 1e-9, thres_ini = 1e-3)

#---------------
# 200 locations
#---------------

#r: 2 
#condition number of C 8.680415e+17 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 1.198395e+17 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 1.452322e+17 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 1.192316e+17 
#r 5 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


#--------------
# 20 locations
#--------------
# r: 2 
#condition number of C 9656854 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 17817090 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 24990765 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 27624647 
#r 5 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 





#======
# Note
#======
as.gpu.matrix(matrix(rnorm(9), 3, 3)) # as.gpu.matrix(matrix obj)

try <- matrix(c(4, 1, 1, 2), 2, 2)
str(forceSymmetric(try))
# Formal class 'dsyMatrix' [package "Matrix"] with 5 slots
#..@ Dim     : int [1:2] 2 2
#..@ Dimnames:List of 2
#.. ..$ : NULL
#.. ..$ : NULL
#..@ x       : num [1:4] 4 1 1 2
#..@ uplo    : chr "U"
#..@ factors : list()
 
### so need to as.matrix of dsymatrix first before as.gpu.matrix


## cbind, rbind cannot be performed between dsymatrix and gpumatrix
# but only between matrix and dsymatrix, or matrix and gpumatrix



#------
# draft
#------
#I_sps <- diag(1, nrow(H_adj), nrow(H_adj))
#I_sps_gpu <- as.gpu.matrix(I_sps)
#chol_inv_gpu(I_sps_gpu)

#as.gpu.matrix(diag(1/2, nrow(H_adj), nrow(H_adj)))



