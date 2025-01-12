#==============
# 11 Dec. 2024
#==============
# Aim: 
# 1. Generate joint Sigma and Sigma_inv for 1st-stage plots, CI among p only
# 2. Calculate elapsed system wall time


# Settings:
# p = 10, n = 400, 600, 800
# pars: A = 0.1, dlt = 0.5, sig2 = 1, kappa = 2


# Method:
# TST9c in 032c



TST9c_SpNormPert_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat,
                                      d_vec, h, reg_ini = 1e-9, thres_ini = 1e-3) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Pert_Mat.R")
  source("Fn_I_sparse.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
  source("Fn_Thres_tune_cov.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  
  C11 <- Matern_32(Var = sig2_mat[1, 1], Kappa = kappa_mat[1, 1], d_vec = d_vec)
  n <- nrow(C11)
  SIGMA <- C11
  
  reg_num <- reg_ini # Initialize reg_num
  
  restart <- T
  while (restart) {
    restart <- F
    SIGMA_inv_pd <- TRUE  # Ini flag for p.d. of SIGMA_inv
    
    # Reset SIGMA at each iteration
    SIGMA <- C11
    
    for(r in seq(2, p, by = 1)){
      
      PN <- Check_par_node(Node = r, data = data)
      R <- C <- NULL
      
      cat("r:", r, "\n")
      for(c in seq(1, (r-1), by = 1)){
        
        BT <- NULL
        C_rc <- 0
        for(t in c(PN)){
          B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
          #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
          
          ## spectral normalization of B_rt
          B_rt <- check_set_SpNorm_Reg(B_rt, reg_num = reg_num)
          cat("B cond numb:", kappa(B_rt), "\n")
          
          BT <- rbind(BT, t(B_rt))
          C_rc <- C_rc + B_rt %*% SIGMA[((t-1)*n+1) : (t*n), ((c-1)*n+1): (c*n)]
        }
        
        R <- cbind(R, C_rc)
        
        C_cr <- t(C_rc)
        C <- rbind(C, C_cr)
      }
      
      D_rr <- Matern_32(Var = sig2_mat[r, r], Kappa = kappa_mat[r, r], d_vec = d_vec)
      
      t <- c(PN)
      Subset_cols <- function(t) {
        start_col <- (t - 1) * n + 1
        end_col <- t * n
        
        result <- R[, start_col:end_col]
      }
      
      result_lst <- lapply(t, FUN = Subset_cols)
      R_subset <- do.call(cbind, result_lst)
      
      C_rr <- R_subset %*% BT + D_rr
      
      
      SG <- SIGMA # p-1 for inverse
      SG <- forceSymmetric(SG)
      
      Col <- rbind(C, C_rr)
      Row <- rbind(SG, R)
      SIGMA <- cbind(Row, Col)
      SIGMA <- forceSymmetric(SIGMA)
      
      
      ## SIGMA_Inv
      Drr_inv <- chol2inv(chol(D_rr)) # Schur complement inverse
      
      if (r == 2){
        SG_inv <- chol2inv(chol(SG))
      }
      
      CDrr_in <- C %*% Drr_inv
      CDR_sym <- forceSymmetric(CDrr_in %*% R)
      
      cat("condition number of C", kappa(C), "\n")
      cat("condition number of CDinv", kappa(CDrr_in), "\n")
      cat("condition number of CDinvR", kappa(CDR_sym), "\n")
      
      BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
      BK1 <- forceSymmetric(BK1)
      BK2 <- - SG_inv %*% C %*% Drr_inv
      BK3 <- - Drr_inv %*% R %*% SG_inv
      BK4 <- forceSymmetric(Drr_inv)
      
      SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
      SG_inv <- forceSymmetric(SIGMA_inv)
      
      # Check if SIGMA_inv is positive definite
      if (!check_pd(SG_inv)) {
        restart <- T
        SIGMA_inv_pd <- F
        cat("New reg_num need.", "\n")
        
        if (reg_num < 1){
          reg_num <- reg_num * 10 # update regularize number
          cat("Reg_num updated to:", reg_num, "\n")
          
        } else {
          reg_num <- reg_num + 0.1
          cat("Reg_num updated to:", reg_num, "\n")
        }
        
        break
      }
      
      cat("r", r, "\n")
      cat("SG_inv", "\n")
      Tst_sym_pd(SG_inv)
      
    }
    
    cat("Final reg_num:", reg_num, "\n")
    
    # Check if SIGMA_inv is positive definite for all r
    if (SIGMA_inv_pd) {
      
      # Compute SIGMA_inv_ini, SG_inv constructed with the smallest possible reg_num
      SIGMA_inv_ini <- SG_inv * (abs(SG_inv) > thres_ini)
      
      # 1. tune threshold if SIGMA_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIGMA_inv and SIGMA
      SG_SG_inv_thres <- Thres_tune_cov(thres_ini = thres_ini, cov_mat_thres = SIGMA_inv_ini,
                                        cov_mat = SG_inv) 
      
      return(list(SIGMA = as.matrix(SIGMA),
                  SIGMA_inv = SG_SG_inv_thres$SIGMA_inv))
    }
    
    
    # check for condition to terminate while loop
    #if (check_pd(SG_SG_inv_thres$SIGMA_inv)){
    # break
    #}
    
  }
}


#=========
# Settings
#=========

#------------------------------------
# Location, displacements, distance
#------------------------------------

## Locations
#ds <- 0.1 # for with SpN plts
ds <- 0.05 # for w/o SpN plts, esp. SIGMA; also try for plts with SpN
# ds = 0.05 has better visualization effect for both with, w/o SpN

s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20]; num [1:40]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:400]


s <- seq(-15 + ds/2, 15 - ds/2, by = ds)
str(s) # num [1:600]

s <- seq(-20 + ds/2, 20 - ds/2, by = ds)
str(s) # num [1:800]


s <- seq(-25 + ds/2, 25 - ds/2, by = ds)
str(s) # num [1:1000]


## Displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) 

## Distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) 


#----------------
# data structure
#----------------

p = 10 
hierarchy_data10 <- data.frame(
  node_id = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10),
  par_id = c(NA, 1, 2, c(2, 3), c(1,4), c(1, 5), 3, 6, 6, 6)
)


#-----------
# Parameters
#-----------

source("Fn_para_mat_construct.R")
all_pars_lst_10 <- All_paras(p = 10, data = hierarchy_data10)


source("Fn_set_ini_vals.R")  
A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[1]], ini_vals = 0.1)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[4]], ini_vals = 2)



#===============================================
# Generate Joint Sigma and Sigma_inv in parallel
#===============================================

#------
# p = 10, n = 400, Tri-wave
#-----
SG_SG_inv_TriWave_p10_n400 <- TST9c_SpNormPert_SG_SGInv(p = 10, 
                                                        data = hierarchy_data10, 
                                                        A_mat = A_mat_0.1, 
                                                        dlt_mat = dlt_mat_0.5, 
                                                        sig2_mat = sig2_mat_1, 
                                                        kappa_mat = kappa_mat_2,
                                                        d_vec = D_vec, h = H, 
                                                        reg_ini = 1e-9, 
                                                        thres_ini = 1e-3)

#saveRDS(SG_SG_inv_TriWave_p10_n400, file = "SG_SG_inv_TriWave_p10_n400.rds",
#        compress = T)

#file_info <- file.info("SG_SG_inv_TriWave_p10_n400.rds")

#file_info$size
# [1] 123548975 = 123.5MB

SG_SG_inv_TriWave_p10_n400 <- readRDS("SG_SG_inv_TriWave_p10_n400.rds")



#------
# p = 10, n = 600, Tri-wave
#-----
system.time(
SG_SG_inv_TriWave_p10_n600 <- TST9c_SpNormPert_SG_SGInv(p = 10, 
                                                        data = hierarchy_data10, 
                                                        A_mat = A_mat_0.1, 
                                                        dlt_mat = dlt_mat_0.5, 
                                                        sig2_mat = sig2_mat_1, 
                                                        kappa_mat = kappa_mat_2,
                                                        d_vec = D_vec, h = H, 
                                                        reg_ini = 1e-9, 
                                                        thres_ini = 1e-3)
  )


#saveRDS(SG_SG_inv_TriWave_p10_n600, file = "SG_SG_inv_TriWave_p10_n600.rds")

#file_info <- file.info("SG_SG_inv_TriWave_p10_n600.rds")
#file_info$size
# [1] 268662749 = 268.7MB

# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 

##############################
#user       system  elapsed 
#2883.046   23.923 5231.029 
##############################


#------
# p = 10, n = 800, Tri-wave
#-----
SG_SG_inv_TriWave_p10_n800 <- TST9c_SpNormPert_SG_SGInv(p = 10, 
                                                        data = hierarchy_data10, 
                                                        A_mat = A_mat_0.1, 
                                                        dlt_mat = dlt_mat_0.5, 
                                                        sig2_mat = sig2_mat_1, 
                                                        kappa_mat = kappa_mat_2,
                                                        d_vec = D_vec, h = H, 
                                                        reg_ini = 1e-9, 
                                                        thres_ini = 1e-3)


# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 

#saveRDS(SG_SG_inv_TriWave_p10_n800, file = "SG_SG_inv_TriWave_p10_n800.rds",
#       compress = T)

#ff <- file.info("SG_SG_inv_TriWave_p10_n800.rds")
#ff$size
# [1] 458903966 = 458.9MB

SG_SG_inv_TriWave_p10_n800 <- readRDS("SG_SG_inv_TriWave_p10_n800.rds")



#==================
# SG, SG_inv plots
#==================






