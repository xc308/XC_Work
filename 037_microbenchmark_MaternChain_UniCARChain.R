#============
# 5 Jan. 2024
#============
# Aim:
  # Want to compare the computational speed between
  # Matern + Chain and Uni-CAR + Chain
  

# Method:
  # microbenchmark the SIGMA, SIGMA_inv generation time
  # using TST10(uni-CAR + Chain) and TST11 (Matern + Chain)


# 
#=========
# Settings
#=========

install.packages("microbenchmark")
library(microbenchmark)


#---------------
# data structure
#---------------

p = 6
hierarchy_data6_chain <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 6),
  par_id = c(NA, 1, 2, 3, 4, 5)
)



p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1
#s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
#str(s) # num [1:20]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]; num [1:200, 1:200]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]; num [1:40000]


#--------------------------------------
# Sepration lag based adjacency matrix
#--------------------------------------

# separation lag
abs(H)

# radius for definition of neighbourhood
abs(H) < 0.4 # 3-order

as.numeric(abs(H) < 0.4) # vector
H_adj <- matrix(as.numeric(abs(H) < 0.4), nrow(H), nrow(H))
diag(H_adj) <- 0


#H_adj <- matrix(as.numeric(abs(H) < 0.3), nrow(H), nrow(H))
#diag(H_adj) <- 0 # 2nd-order


#---------------------
# phi and p.d.of H_adj
#---------------------

eigen_Hadj <- eigen(H_adj, symmetric = T, only.values = T)$val
1/ max(abs(eigen_Hadj)) # [1] 0.1577773; [1] 0.1412107

#phi = 0.15
phi = 0.14




#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6_chain)
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6)


source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)



#================
# Test functions
#================


TST10_UniCAR_Chain <- function(p, data, chain = F, A_mat, dlt_mat, 
                                      sig2_mat, phi, H_adj, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Pert_Mat.R")
  source("Fn_Pert_Mat_2.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpN_Pert.R")# lower kappa of B
  source("Fn_I_sparse.R")
  
  # phi: controls the propostion of conditional spatial dependence 
  # sig2_mat: diag matrix containing sigma2 in CAR
  # H: adjacency matrix
  
  I_sps <- I_sparse(size = nrow(H_adj), value = 1)
  c_inv <- I_sps - phi * H_adj
  C11_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = 1/sig2_mat[1, 1])
  
  C11 <- chol2inv(chol(C11_inv)) 
  n <- nrow(C11)
  SIGMA <- C11
  
  for(r in seq(2, p, by = 1)){
    
    PN <- Check_par_node(Node = r, data = data)
    R <- C <- NULL
    
    #cat("r:", r, "\n")
    for(c in seq(1, (r-1), by = 1)){
      
      
      BT <- NULL
      C_rc <- 0
      for(t in c(PN)){
        #B_rt <- wave_v9(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v7(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
        
        ## spectral normalization of B_rt
        B_rt <- check_set_SpNorm_pert(B_rt)
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
    #cat("condition number of C", kappa(C), "\n")
    #cat("condition number of CDinv", kappa(CDrr_in), "\n")
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
    
    
    # early perturb SG_inv if not p.d.
    #cat("r", r, "\n")
    #SG_inv <- Pert_Mat(SG_inv)
    SG_inv <- Pert_Mat_2(SG_inv)
    
    
    #cat("SG_inv", "\n")
    #Tst_sym_pd(SG_inv)
    
    
    
    if (r == p){
      if (chain){
        for(bdw in seq(0.3, 0.95, by = 0.05)){
          SG_inv_bd <- band(SG_inv, -bdw * nrow(SG_inv), bdw * ncol(SG_inv))
          
          if (all(eigen(SG_inv_bd, symmetric = T, only.values = T)$val > 0)){
            SG_inv <- as(SG_inv_bd, "sparseMatrix")
            
            cat("bdw：", bdw, "\n")
            break
          }
        }
        return(
          list(SIGMA = as.matrix(SIGMA), 
               #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
               SIGMA_inv = SG_inv # sparseMatrix formate, easy for inference
          )
        )
      }  
      
      return(
        list(SIGMA = as.matrix(SIGMA), 
             #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
             SIGMA_inv = as.matrix(SG_inv))
      )
    }
  }
}



TST11_Matern_Chain <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat,
                                      d_vec, h, chain = F) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Pert_Mat.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpN_Pert.R")# lower kappa of B
  
  C11 <- Matern_32(Var = sig2_mat[1, 1], Kappa = kappa_mat[1, 1], d_vec = d_vec)
  n <- nrow(C11)
  SIGMA <- C11
  
  for(r in seq(2, p, by = 1)){
    
    PN <- Check_par_node(Node = r, data = data)
    R <- C <- NULL
    
    #cat("r:", r, "\n")
    for(c in seq(1, (r-1), by = 1)){
      
      
      BT <- NULL
      C_rc <- 0
      for(t in c(PN)){
        #B_rt <- wave_v9(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v7(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
        
        ## spectral normalization of B_rt
        B_rt <- check_set_SpNorm_pert(B_rt)
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
    
    #CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
    #cat("condition number of C", kappa(C), "\n")
    #cat("condition number of CDinv", kappa(CDrr_in), "\n")
    #cat("condition number of CDinvR", kappa(CDR_sym), "\n")
    
    BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
    BK1 <- forceSymmetric(BK1)
    BK2 <- - SG_inv %*% C %*% Drr_inv
    BK3 <- - Drr_inv %*% R %*% SG_inv
    BK4 <- forceSymmetric(Drr_inv)
    
    #cat("BK4", "\n")
    #Tst_sym_pd(BK4)
    
    
    SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
    #SG_inv <- SIGMA_inv
    SG_inv <- forceSymmetric(SIGMA_inv)
    
    
    # early perturb SG_inv if not p.d.
    #SG_inv <- Pert_Mat(SG_inv)
    SG_inv <- Pert_Mat_2(SG_inv)
    
    
    #cat("SG_inv", "\n")
    #Tst_sym_pd(SG_inv)
    
    
    
    if (r == p){
      if (chain){
        for(bdw in seq(0.2, 0.95, by = 0.05)){
          SG_inv_bd <- band(SG_inv, -bdw * nrow(SG_inv), bdw * ncol(SG_inv))
          
          if (all(eigen(SG_inv_bd, symmetric = T, only.values = T)$val > 0)){
            SG_inv <- as(SG_inv_bd, "sparseMatrix")
            
            cat("bdw：", bdw, "\n")
            break
          }
        }
        return(
          list(SIGMA = as.matrix(SIGMA), 
               #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
               SIGMA_inv = SG_inv # sparseMatrix formate, easy for inference
          )
        )
      }  
      
      return(
        list(SIGMA = as.matrix(SIGMA), 
             #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
             SIGMA_inv = as.matrix(SG_inv))
      )
    }
  }
}




TST12_UniTap_Chain <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, 
                                  d_vec, h, chain = F) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Pert_Mat.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpN_Pert.R")# lower kappa of B
  
  C11 <- Matern_32(Var = sig2_mat[1, 1], Kappa = kappa_mat[1, 1], d_vec = d_vec)
  Tap <- WendLd_32(r = h, R = 0.5, dlt = 0, A = 1)
  C11 <- C11 * Tap
  
  
  n <- nrow(C11)
  SIGMA <- C11
  
  for(r in seq(2, p, by = 1)){
    
    PN <- Check_par_node(Node = r, data = data)
    R <- C <- NULL
    
    cat("r:", r, "\n")
    for(c in seq(1, (r-1), by = 1)){
      
      
      BT <- NULL
      C_rc <- 0
      for(t in c(PN)){
        #B_rt <- wave_v9(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v7(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
        
        ## spectral normalization of B_rt
        B_rt <- check_set_SpNorm_pert(B_rt)
        #cat("B cond numb:", kappa(B_rt), "\n")
        
        BT <- rbind(BT, t(B_rt))
        C_rc <- C_rc + B_rt %*% SIGMA[((t-1)*n+1) : (t*n), ((c-1)*n+1): (c*n)]
      }
      
      R <- cbind(R, C_rc)
      
      C_cr <- t(C_rc)
      C <- rbind(C, C_cr)
    }
    
    D_rr <- Matern_32(Var = sig2_mat[r, r], Kappa = kappa_mat[r, r], d_vec = d_vec)
    Tap <- WendLd_32(r = h, R = 0.5, dlt = 0, A = 1)
    D_rr <- D_rr * Tap
    
    
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
    
    #CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
    #cat("condition number of C", kappa(C), "\n")
    #cat("condition number of CDinv", kappa(CDrr_in), "\n")
    #cat("condition number of CDinvR", kappa(CDR_sym), "\n")
    
    BK1 <- SG_inv %*% (SG + CDR_sym) %*% SG_inv
    BK1 <- forceSymmetric(BK1)
    BK2 <- - SG_inv %*% C %*% Drr_inv
    BK3 <- - Drr_inv %*% R %*% SG_inv
    BK4 <- forceSymmetric(Drr_inv)
    
    #cat("BK4", "\n")
    #Tst_sym_pd(BK4)
    
    
    SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
    #SG_inv <- SIGMA_inv
    SG_inv <- forceSymmetric(SIGMA_inv)
    
    # early perturb SG_inv if not p.d.
    #cat("r", r, "\n")
    #SG_inv <- Pert_Mat(SG_inv)
    SG_inv <- Pert_Mat_2(SG_inv)
    
    
    #cat("SG_inv", "\n")
    #Tst_sym_pd(SG_inv)
    
    
    if (r == p){
      if (chain){
        for(bdw in seq(0.2, 0.95, by = 0.05)){
          SG_inv_bd <- band(SG_inv, -bdw * nrow(SG_inv), bdw * ncol(SG_inv))
          
          if (all(eigen(SG_inv_bd, symmetric = T, only.values = T)$val > 0)){
            SG_inv <- as(SG_inv_bd, "sparseMatrix")
            
            cat("bdw：", bdw, "\n")
            break
          }
        }
        return(
          list(SIGMA = as.matrix(SIGMA), 
               #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
               SIGMA_inv = SG_inv # sparseMatrix formate, easy for inference
          )
        )
      }  
      
      return(
        list(SIGMA = as.matrix(SIGMA), 
             #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
             SIGMA_inv = as.matrix(SG_inv))
      )
    }
  }
}





#================
# Microbenchmark
#================

SG_SG_inv_UniCAR_Chain <- TST10_UniCAR_Chain(p = 6, data = hierarchy_data6_chain, 
                   chain = T, A_mat = A_01, dlt_mat = dlt_05,
                   sig2_mat = sig2_mat_1, phi = phi, 
                   H_adj = H_adj, h = H)



SG_SG_inv_Matern_Chain <- TST11_Matern_Chain(p = 6, data = hierarchy_data6_chain, 
                   A_mat = A_01, dlt_mat = dlt_05,
                   sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                   d_vec = D_vec, h = H, chain = T)


SG_SG_inv_Tap_Chain <- TST12_UniTap_Chain(p = 6, data = hierarchy_data6_chain, 
                   A_mat = A_01, dlt_mat = dlt_05,
                   sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                   d_vec = D_vec, h = H, chain = T)




microbenchmark(
  SG_SG_inv_UniCAR_Chain <- TST10_UniCAR_Chain(p = 6, data = hierarchy_data6_chain, 
                                               chain = T, A_mat = A_01, dlt_mat = dlt_05,
                                               sig2_mat = sig2_mat_1, phi = phi, 
                                               H_adj = H_adj, h = H),
  
  SG_SG_inv_Tap_Chain <- TST12_UniTap_Chain(p = 6, data = hierarchy_data6_chain, 
                                            A_mat = A_01, dlt_mat = dlt_05,
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H, chain = T),
  
  SG_SG_inv_Matern_Chain <- TST11_Matern_Chain(p = 6, data = hierarchy_data6_chain, 
                                               A_mat = A_01, dlt_mat = dlt_05,
                                               sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                               d_vec = D_vec, h = H, chain = T)
  
   

)


# Results:
#      min       lq     mean    median       uq      max    neval
# 2.839134 2.962440 2.978158  2.969599 2.987456 3.104558    100
# 4.013218 4.130731 4.131932  4.141981 4.149930 4.246715    100 
# 4.016198 4.140458 4.165273  4.153939 4.180296 4.464462    100


#============================================
# Mircrobenchmark UniCAR and Taper NOT chain
#============================================



SG_SG_inv_UniCAR_NC <- TST10_UniCAR_Chain(p = 6, data = hierarchy_data6, 
                   chain = F, A_mat = A_01, dlt_mat = dlt_05,
                   sig2_mat = sig2_mat_1, phi = phi, 
                   H_adj = H_adj, h = H)



SG_SG_inv_Tap_NC <- TST12_UniTap_Chain(p = 6, data = hierarchy_data6, 
                                          A_mat = A_01, dlt_mat = dlt_05,
                                          sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                          d_vec = D_vec, h = H, chain = F)


microbenchmark(
  SG_SG_inv_UniCAR_NC <- TST10_UniCAR_Chain(p = 6, data = hierarchy_data6, 
                                            chain = F, A_mat = A_01, dlt_mat = dlt_05,
                                            sig2_mat = sig2_mat_1, phi = phi, 
                                            H_adj = H_adj, h = H),
  
  
  SG_SG_inv_Tap_NC <- TST12_UniTap_Chain(p = 6, data = hierarchy_data6, 
                                         A_mat = A_01, dlt_mat = dlt_05,
                                         sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                         d_vec = D_vec, h = H, chain = F)
  
  
  
)


# Results
#      min       lq     mean     median       uq      max neval
# 2.778585 2.897267 2.898103   2.903990 2.912605 3.043995   100
# 3.879586 3.902875 3.975860   4.003598 4.012794 4.154892   100








