#=============
# 11 Mar. 2024
#=============

# Aim: 
  # 1. Obtain the distance matrix in 2D
  # 2. Construct the H_adj matrix for UniCAR
  # 3. Modify the TST function for 2D


# Method:
  # Modify on top of TST10c in 034c




TST12_SG_SGInv_CAR_2D <- function(p, data, A_mat, dsp_lon_mat, dsp_lat_mat, 
                                  dlt_lon_mat, dlt_lat_mat, b = "Wendland",
                                  phi, H_adj, sig2_mat,
                                  reg_ini = 1e-9, thres_ini = 1e-3) {
  
  #source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
  source("Fn_I_sparse.R")
  source("Fn_Thres_tune_cov.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  source("Fn_shft_dist_mat.R") # construct shifted distance matrix using shft displacement for b function
  
  # phi: controls the propostion of conditional spatial dependence 
  # sig2_mat: diag matrix containing sigma2 in CAR
  # H: adjacency matrix
  
  I_sps <- I_sparse(size = nrow(H_adj), value = 1)
  c_inv <- I_sps - phi * H_adj
  C11_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = 1/sig2_mat[1, 1])
  
  C11 <- chol2inv(chol(C11_inv)) 
  n <- nrow(C11)
  SIGMA <- C11
  
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
            #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
          }
          
          if (b == "Wendland") {
            
            shft_dst_mat <- Shft_dst_mat(dsp_lon_mat = dsp_lon_mat, dsp_lat_mat = dsp_lat_mat, 
                                         dlt1 = dlt_lon_mat[r, t], dlt2 =  dlt_lat_mat[r, t])
            
            B_rt <- WendLd32_2D(shft_dst_mat = shft_dst_mat, A = A_mat[r, t])
            #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
          }
          
          ## spectral normalization of B_rt
          B_rt <- check_set_SpNorm_Reg(B_rt, reg_num = reg_num)
          #cat("B cond numb:", kappa(B_rt), "\n")
          
          BT <- rbind(BT, t(B_rt))
          C_rc <- C_rc + B_rt %*% SIGMA[((t-1)*n+1) : (t*n), ((c-1)*n+1): (c*n)]
        }
        
        R <- cbind(R, C_rc)
        
        C_cr <- t(C_rc)
        C <- rbind(C, C_cr)
      }
      
      Drr_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = 1/sig2_mat[r, r])
      D_rr <- chol2inv(chol(Drr_inv))
      
      
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
      #Drr_inv <- chol2inv(chol(D_rr)) # Schur complement inverse
      
      if (r == 2){
        SG_inv <- C11_inv
      }
      
      CDrr_in <- C %*% Drr_inv
      #CDR_sym <- forceSymmetric(CDrr_in %*% R)
      
      #CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
      cat("condition number of C", kappa(C), "\n")
      cat("condition number of CDinv", kappa(CDrr_in), "\n")
      #cat("condition number of CDinvR", kappa(CDR_sym), "\n")
      
      SGCD <- SG_inv %*% CDrr_in
      RSG <- R %*% SG_inv
      BK1 <- SG_inv + SGCD %*% RSG
      BK1 <- forceSymmetric(BK1)
      
      #BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
      #BK1 <- forceSymmetric(BK1)
      BK2 <- - SG_inv %*% (C %*% Drr_inv)
      BK3 <- - (Drr_inv %*% R) %*% SG_inv
      BK4 <- forceSymmetric(Drr_inv)
      
      #cat("BK4", "\n")
      #Tst_sym_pd(BK4)
      
      
      SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
      #SG_inv <- SIGMA_inv
      SG_inv <- forceSymmetric(SIGMA_inv)
      
      # to condition to break the for loop and start all over again
      if (!check_pd(SG_inv)) {
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
      Tst_sym_pd(SG_inv)
    }
    
    cat("Final reg_num:", reg_num, "\n")
    
    if (SIGMA_inv_pd){
      
      # Compute SIGMA_inv_ini, SG_inv constructed with the smallest possible reg_num
      SIGMA_inv_ini <- SG_inv * (abs(SG_inv) > thres_ini)
      
      # 1. tune threshold if SIGMA_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIGMA_inv and SIGMA
      SG_SG_inv_thres <- Thres_tune_cov(thres_ini = thres_ini, 
                                        cov_mat_thres = SIGMA_inv_ini,
                                        cov_mat = SG_inv) 
      
      return(list(SIGMA = as.matrix(SIGMA),
                  SIGMA_inv = SG_SG_inv_thres$SIGMA_inv))  
      
    }
  }
}


#============
# Test TST12
#============

#----------
# 2D coords
#----------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
crds <- cbind(s, s)
head(crds)
#         s     s
#[1,] -0.95 -0.95
#[2,] -0.85 -0.85
#[3,] -0.75 -0.75
#[4,] -0.65 -0.65
#[5,] -0.55 -0.55

nrow(crds) #[1] 20


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
    # for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
str(DSP[, , 1])
# num [1:20, 1:20]
str(DSP[, , 2])
# num [1:20, 1:20] 


#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
  # for H_adj and phi for UniCAR

DIST <- as.matrix(dist(crds, diag = T, upper = T))
quantile(DIST)
str(DIST) # num [1:20, 1:20]

## set Nb_radius = 0.6
Nb_radius <- 0.6


H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0
H_adj


spec <- eigen(H_adj, symmetric = T, only.values = T)$val
phi <- 1/max(abs(spec)) # [1] 0.1344431

phi <- trunc(phi * 100)/100  # [1] 0.13



#--------
# pars
#--------

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_6 <- All_paras_CAR_2D(p = 6, data =  hierarchy_data6)

source("Fn_set_ini_vals.R")
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[1]], ini_vals = 1)
dlt_lon_02 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[2]], ini_vals = 0.2)
dlt_lat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[3]], ini_vals = 0.4)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[4]], ini_vals = 1)



#------
# TST12
#------

## Tri-Wave
SG_SGinv_CAR_6_2D_TW <- TST12_SG_SGInv_CAR_2D(p = 6, data = hierarchy_data6, 
                      A_mat = A_1, dsp_lon_mat = DSP[, , 1], dsp_lat_mat = DSP[, , 2],
                      dlt_lon_mat = dlt_lon_02, dlt_lat_mat = dlt_lat_04, 
                      b = "Tri-Wave", phi =  phi, H_adj = H_adj,
                      sig2_mat = sig2_mat_1, reg_ini = 1e-9, thres_ini = 1e-3)



# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


## Wendland
SG_SGinv_CAR_6_2D_WL <- TST12_SG_SGInv_CAR_2D(p = 6, data = hierarchy_data6,
                      A_mat = A_1, 
                      dsp_lon_mat = DSP[, , 1], dsp_lat_mat = DSP[, , 2], 
                      dlt_lon_mat = dlt_lon_02, 
                      dlt_lat_mat = dlt_lat_04,
                      b = "Wendland", 
                      phi = phi, H_adj = H_adj, sig2_mat = sig2_mat_1,
                      reg_ini = 1e-9, thres_ini = 1e-3)

# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 





