#=============
# 8 Mar. 2024
#=============
# Aim:  
  # Function of SG, SG_inv generation using C.I. among p variate fields 
    # (non-cross-MRF, i.e. Matern + DAG)

# Aim:
# Function for generating SG, SG_inv using cross-MRF (CAR+DAG).

# Args:
    # p : variate field number
    # data: data structure reflecting DAG
    # A_mat: matrix of par A, amplitude in b function (TW, WL)
    # dlt_mat: matrix of par dlt, shift in b function
    # sig2_mat: matrix of par sig2 in Matern
    # kappa_mat: matrix of par kappa in Matern
    # h: matrix of displacement
    # d_vec: a vector of distance
    # b: b function, Wendland or Tri-Wave
    # reg_ini: regularization initial value to ensure matrix B numerically stable
    # thres_ini: the ini value to threshold the value of SG_inv to exact zero





TST9d_SpNormReg_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat,
                                      d_vec, h, b = "Wendland", reg_ini = 1e-9, 
                                      thres_ini = 1e-3) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  #source("Fn_Pert_Mat.R")
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
          
          if (b == "Tri-Wave") {
            B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
          }
          
          if (b == "Wendland") {
            B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
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
      
      #cat("condition number of C", kappa(C), "\n")
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
    
    
  }
}