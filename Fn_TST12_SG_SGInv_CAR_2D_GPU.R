#=============
# 18 Apr. 2024
#=============

# Aim:
  # Function TST12_SG_SGInv_CAR_2D_GPU for neg_logL source



TST12_SG_SGInv_CAR_2D_GPU <- function(p, data, A_mat, dsp_lon_mat, dsp_lat_mat, 
                                      dlt_lon_mat, dlt_lat_mat, b = "Tri-Wave",
                                      phi, H_adj, sig2_mat,
                                      reg_ini = 1e-9, thres_ini = 1e-3) {
  
  #source("Fn_Check_par_node.R")
  #source("Fn_Waves.R")
  #source("Fn_Wendland_32.R") # R = 0.5
  #source("Fn_Tst_sym_pd_GPU.R") # GPU version
  #source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
  #source("Fn_I_sparse.R")
  #source("Fn_Thres_tune_cov_GPU.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv on GPU
  #source("Fn_shft_dist_mat.R") # construct shifted distance matrix using shft displacement for b function
  #source("Fn_forceSym_GPU.R") # forceSym on GPU
  
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
        SG_inv_gpu <- as.gpu.matrix(as.matrix(SG_inv), device = "cuda")
      }
      
      # off-load relevant matrices on GPU
      SG_gpu <- as.gpu.matrix(as.matrix(SG), device = "cuda")
      C_gpu <- as.gpu.matrix(C, device = "cuda")
      Drr_inv_gpu <- as.gpu.matrix(as.matrix(Drr_inv), device = "cuda")
      R_gpu <- as.gpu.matrix(R, device = "cuda")
      
      
      # Construct each BK on GPU
      CDrr_in_gpu <- C_gpu %*% Drr_inv_gpu
      SGCD_gpu <- SG_inv_gpu %*% CDrr_in_gpu
      RSG_gpu <- R_gpu %*% SG_inv_gpu
      BK1 <- SG_inv_gpu + SGCD_gpu %*% RSG_gpu
      BK1 <- forceSym_gpu(BK1)
      
      #BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
      #BK1 <- forceSymmetric(BK1)
      BK2 <- - SG_inv_gpu %*% (C_gpu %*% Drr_inv_gpu)
      BK3 <- - (Drr_inv_gpu %*% R_gpu) %*% SG_inv_gpu
      BK4 <- forceSym_gpu(Drr_inv_gpu)
      
      #cat("BK4", "\n")
      #Tst_sym_pd(BK4)
      
      
      SIGMA_inv_gpu <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
      #SG_inv <- SIGMA_inv
      SG_inv_gpu <- forceSym_gpu(SIGMA_inv_gpu)
      
      # condition to break the for loop and start all over again
      if (!check_pd_gpu(SG_inv_gpu)) {
        restart <- T
        SIGMA_inv_pd <- F
        cat("New reg_num needed.", "\n")
        
        if (reg_num < 0.1) {
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
      cat("SG_inv_gpu", "\n")
      Tst_sym_pd_gpu(SG_inv_gpu)
    }
    
    cat("Final reg_num:", reg_num, "\n")
    
    if (SIGMA_inv_pd){
      
      # Compute SIGMA_inv_ini, SG_inv constructed with the smallest possible reg_num
      SIGMA_inv_ini <- SG_inv_gpu * (abs(SG_inv_gpu) > thres_ini)
      #class(SIGMA_inv_ini) 
      #Formal class 'gpu.matrix.torch' [package "GPUmatrix"] with 5 slots
      #@ gm      :Double [1:100, 1:100]
      
      
      # 1. tune threshold if SIGMA_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIGMA_inv and SIGMA
      SG_SG_inv_thres <- Thres_tune_cov_gpu(thres_ini = thres_ini, 
                                            cov_mat_thres = SIGMA_inv_ini,
                                            cov_mat = SG_inv_gpu) 
      
      return(list(SIGMA_gpu = as.gpu.matrix(as.matrix(SIGMA), device = "cuda"),
                  SIGMA_inv_gpu = SG_SG_inv_thres$SIGMA_inv_gpu))  
      
    }
  }
}












