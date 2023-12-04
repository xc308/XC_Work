#==============
# 24 Nov. 2023
#==============

# wave function versions:
# wave_v4: TST; phi = 1/2; /h-dlt/ < 2/dlt/
# wave_v5: TST2; phi = 2; /h-dlt/ < /dlt/
# wave_v6: TST3; phi = 1/2; /h-dlt/ < /dlt/
# wave_v7: TST4: phi = 2; /h-dlt/ < 2/dlt/


# Aim:
  # want to know the impact of domain of compact spprt
    # on the p.d. of SIGMA_inv; 
    # prefer larger domain or smaller domain?

    # when wave has the same decay speed parameter phi = 2
    # but wave_v7 vs wave_v5


# Conclusion:
    # wave_v7's value range [-7, 1], too large oscilation
    # too many negative values??
 
   
source("Fn_Waves.R")

TST4_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  
  C11 <- Matern_32(Var = sig2_mat[1, 1], Kappa = kappa_mat[1, 1], d_vec = d_vec)
  n <- nrow(C11)
  SIGMA <- C11
  
  for(r in seq(2, p, by = 1)){
    
    PN <- Check_par_node(Node = r, data = data)
    R <- C <- NULL
    
    for(c in seq(1, (r-1), by = 1)){
      
      BT <- NULL
      C_rc <- 0
      for(t in c(PN)){
        B_rt <- wave_v7(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        
        
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
    
    
    ## SIGMA_Inv
    Drr_inv <- chol2inv(chol(D_rr)) # Schur complement inverse
    
    if (r == 2){
      SG_inv <- chol2inv(chol(SG))
    }
    
    CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
    BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
    BK1 <- forceSymmetric(BK1)
    BK2 <- - SG_inv %*% C %*% Drr_inv
    BK3 <- - Drr_inv %*% R %*% SG_inv
    BK4 <- forceSymmetric(Drr_inv)
    
    
    SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
    #SG_inv <- SIGMA_inv
    SG_inv <- forceSymmetric(SIGMA_inv)
    
    
    if (r == p) return(
      list(SIGMA = as.matrix(forceSymmetric(SIGMA)), 
           SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv)))
    )
  }
}


#------
# data
#------
p = 5

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)

source("Fn_set_ini_vals.R")
A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================

# Method:
# use TST4_SG_SGInv with wave_v7 
  # (fast decay phi = 2; larger region supprt vs wave_v5)

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)



for (dlt in seq(0.1, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST4_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                       A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                       sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                       d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    
  }
}


## conslusion 
# all NO
# prefer smaller compact supprt region. 


#=================
# Peturb magnitude
#=================

Tst_Pert_w7 <- function(Dlt, A) {
  
  for (dlt in Dlt){
    cat("dlt", dlt, "\n")
    dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
    
    for (a in A){
      cat("A:", a, "\n")
      A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
      
      SG_SG_inv_5 <- TST4_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                         A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                         sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                         d_vec = D_vec, h = H)
      SG_Inv_5 <- SG_SG_inv_5$SIGMA_inv
      
      for (pert in 10^seq(-7, -1, by = 1)){
        cat("pert:", pert, "\n")
        
        I_pert <- I_spar(size = nrow(SG_Inv_5), value = pert)
        SG_inv_5_pert <- SG_Inv_5 + I_pert
        
        Tst_sym_pd(SG_inv_5_pert)
      }
    }
    
  }
}


Tst_Pert_w7(Dlt = seq(0.1, 0.9, by = 0.1), 
            A = seq(0.5, 1, by = 0.1))

# all NO




















