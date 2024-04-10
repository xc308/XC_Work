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
  
  I_sps <- I_sparse(size = nrow(H_adj), value = 1) 
  I_sps_gpu <- as.gpu.matrix(I_sps, sparse = T)
  c_inv <- I_sps_gpu - phi * H_adj  # gpu
  C11_inv <- c_inv %*% as.gpu.matrix(I_sparse(size = nrow(H_adj), value = 1/sig2_mat[1, 1]))
  
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
      
      Drr_inv <- c_inv %*% as.gpu.matrix(I_sparse(size = nrow(H_adj), value = 1/sig2_mat[r, r])) # gpu
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
      BK1 <- forceSymmetric(as.matrix(BK1))
      
      #BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
      #BK1 <- forceSymmetric(BK1)
      BK2 <- - SG_inv %*% (C %*% Drr_inv)   # gpu
      BK3 <- - (Drr_inv %*% R) %*% SG_inv   # gpu
      BK4 <- forceSymmetric(as.matrix(Drr_inv)) # matrix
      
      #cat("BK4", "\n")
      #Tst_sym_pd(BK4)
      
      
      SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4)) # gpu
      #SG_inv <- SIGMA_inv
      SG_inv <- forceSymmetric(as.matrix(SIGMA_inv))
      SG_inv <- as.gpu.matrix(SG_inv) # gpu 
      
      
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
      SG_SG_inv_thres <- Thres_tune_cov(thres_ini = thres_ini, 
                                        cov_mat_thres = SIGMA_inv_ini,
                                        cov_mat = SG_inv) 
      
      return(list(SIGMA = as.matrix(SIGMA),
                  SIGMA_inv = SG_SG_inv_thres$SIGMA_inv))  
      
    }
  }
}





