#=============
# 26 Feb. 2024
#=============

# Aim:
  # Modify Tst10b allow for choice between Tri-Wave and Wendland function



TST10c_SpNReg_Thres_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat,
                                         phi, H_adj, h, b = "Wendland",
                                         reg_ini = 1e-9, 
                                         thres_ini = 1e-3) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  #source("Fn_Pert_Mat.R")
  #source("Fn_Pert_Mat_2.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
  source("Fn_I_sparse.R")
  source("Fn_Thres_tune_cov.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  
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


#==================
# Test the function
#==================
# Aim:
# to test the above constructed function workable


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6)

source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)


SG_SG_inv_wl_try <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                             A_mat = A_01, dlt_mat = dlt_05, 
                             sig2_mat = sig2_mat_1,
                             phi = phi, H_adj, h = H, 
                             b = "Wendland",
                             reg_ini = 1e-9, thres_ini = 1e-3) 

# r: 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


SG_SG_inv_TW_try <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                 A_mat = A_01, dlt_mat = dlt_05, 
                                                 sig2_mat = sig2_mat_1,
                                                 phi = phi, H_adj, h = H, 
                                                 b = "Tri-Wave",
                                                 reg_ini = 1e-9, thres_ini = 1e-3) 



#r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001  



#===================================
# More parameters combination check
#===================================


for (dlt in seq(0.5, 10, by = 1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = a)
    
    SG_SG_inv_tst_6 <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                                sig2_mat = sig2_mat_1,
                                                phi = phi, H_adj, h = H, 
                                                b = "Wendland",
                                                reg_ini = 1e-9, thres_ini = 1e-3) 
    
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA)
    
  }
}

## Results:
  # all p.d.


for (dlt in seq(10, 20, by = 1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = a)
    
    SG_SG_inv_tst_6 <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                    A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                                    sig2_mat = sig2_mat_1,
                                                    phi = phi, H_adj, h = H, 
                                                    b = "Wendland",
                                                    reg_ini = 1e-9, thres_ini = 1e-3) 
    
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA)
    
  }
}

## results: all p.d.


for (dlt in seq(20, 30, by = 1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = a)
    
    SG_SG_inv_tst_6 <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                    A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                                    sig2_mat = sig2_mat_1,
                                                    phi = phi, H_adj, h = H, 
                                                    b = "Wendland",
                                                    reg_ini = 1e-9, thres_ini = 1e-3) 
    
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_tst_6$SIGMA)
    
  }
}


# result:
# dlt: 21 
# A: 0.5 
# r: 2 
#Error in qr.default(if (d[1L] < d[2L]) t(z) else z) : 
#  NA/NaN/Inf in foreign function call (arg 1)

## conjecture:
  # dlt has upper bound: 
    # dlt <= max(abs(H))

  # i.e., dlt can only shift within the domain


## Investigate: see 034d

quantile(H)
#   0%   25%   50%   75%  100% 
# -19.9  -5.9   0.0   5.9  19.9












