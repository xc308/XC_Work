#============
# 9 Jan. 2024
#============

# Aim:
  # Link cross-MRF and mixed sp graph
  # uni-variate using Matern Taper

# Method:
  # source("038_Uni_Wendland_Tapper.R")
  # TST9_SpNormPert_SG_SGInv



TST12_UniTap_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, 
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
        cat("B cond numb:", kappa(B_rt), "\n")
        
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
    cat("condition number of C", kappa(C), "\n")
    cat("condition number of CDinv", kappa(CDrr_in), "\n")
    cat("condition number of CDinvR", kappa(CDR_sym), "\n")
    
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
    cat("r", r, "\n")
    #SG_inv <- Pert_Mat(SG_inv)
    SG_inv <- Pert_Mat_2(SG_inv)
    
    
    cat("SG_inv", "\n")
    Tst_sym_pd(SG_inv)
    
    
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


#=========
# Settings
#=========

#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:40]

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:40, 1:40]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:1600]


#----------------
# data structure
#----------------

p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


hierarchy_data6_chain <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 6),
  par_id = c(NA, 1, 2, 3, 4, 5)
)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6)
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6_chain)


source("Fn_set_ini_vals.R")
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)


for (dlt in seq(0.1, 1, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = a)
    
    SG_SG_inv_6 <- TST12_UniTap_SG_SGInv(p = 6, data = hierarchy_data6_chain, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_6$SIGMA_inv)
    
  }
}


#================
# Visualisation
#================

A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)

SG_SG_inv_6_Tap_a01d05 <- TST12_UniTap_SG_SGInv(p = 6, data = hierarchy_data6, 
                                               A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                                               sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                               d_vec = D_vec, h = H)




SG_SG_inv_6_Tap_a01d05_chain <- TST12_UniTap_SG_SGInv(p = 6, data = hierarchy_data6_chain, 
                                                A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                                                sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                                d_vec = D_vec, h = H, chain = T)

#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#bdw： 0.25 


#--------------
# Plot Wendland
#--------------

plt_Sig(Sigma = SG_SG_inv_6_Tap_a01d05$SIGMA, p = 6)
plt_Sig(Sigma = log(SG_SG_inv_6_Tap_a01d05$SIGMA), p = 6)
plt_Sig(Sigma = log(SG_SG_inv_6_Tap_a01d05$SIGMA_inv), p = 6)

plt_Sig(Sigma = log(abs(SG_SG_inv_6_Tap_a01d05$SIGMA_inv)), p = 6)


#-------
# Chain
#-------
plt_Sig(Sigma = SG_SG_inv_6_Tap_a01d05_chain$SIGMA, p = 6)
plt_Sig(Sigma = log(SG_SG_inv_6_Tap_a01d05_chain$SIGMA), p = 6)
chain_bd_025 <- as.matrix(SG_SG_inv_6_Tap_a01d05_chain$SIGMA_inv)
plt_Sig(Sigma = log(abs(chain_bd_025)), p = 6)








