#==============
# 22 Nov. 2023
#==============

# Aim:
  # intend to add some forcesymmetric part in 009_1D_simu_NEW
  # but first test here. 

# Modify:
  # TST2_build_SG_SGInv 
    # 1. change the source file to "Fn_Waves.R" contain wave_v5
    # 2. B_rt use wave_v5 (see visualisation of wave_v5.R)
    # 3. forcesymmetric(SG) in Line 67
    # 4. CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R) in Line82
    # 5. BK1, BK4, symmetric, SIGMA_inv symmetric
    
# Note: 1. and 2. are for 012_Investigate_pd_SGInv_TST2_Wave_v5_small_pert.R

# Further reading:
    # 012_Investigate_pd_SGInv_TST2_wave5_small_pert.R
  

TST2_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
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
        B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
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


## below are all based on wave_v4 ##

#======================================
# Test (below are all based on wave_v4)
#======================================
# Aim:
  # want to see for those non-p.d. SIGMA_inv
  # what can be the problems?
    # 1. diagonally dominant? if not, perturb
    # 2. all the diagonal elements > 0?
        # or will they be > 0 after perturb?


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


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1 
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) #[1:400]

    
#-----------------
# Test on the algo
#-----------------
SG_SG_inv_5_TST2 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                    A_mat = A_mat_0.5, dlt_mat = dlt_mat_0.5,
                    sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                    d_vec = D_vec, h = H)


Test_sym_pd(SG_SG_inv_5_TST2$SIGMA)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

Test_sym_pd(SG_SG_inv_5_TST2$SIGMA_inv)
# [1] "Symmetric: Yes"
# [1] "p.d.: No"


#=========================
# try perturb (success!!)
#=========================

I_01 <- I_spar(size = nrow(SG_SG_inv_5_TST2$SIGMA_inv), value = 0.1)
SIGMA_inv_5_pert <- SG_SG_inv_5_TST2$SIGMA_inv + I_01

Test_sym_pd(SIGMA_inv_5_pert)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


#-----------------------------------------
# Now want to know how small I can perturb
#-----------------------------------------

for (pert in 10^seq(-9, -1, by = 1)){
  I_prt <- I_spar(size = nrow(SG_SG_inv_5_TST2$SIGMA_inv), value = pert)
  
  SIGMA_inv_5_pert <- SG_SG_inv_5_TST2$SIGMA_inv + I_prt
  
  cat(pert,"\n")
  Test_sym_pd(SIGMA_inv_5_pert)
}


# 1e-09 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#1e-08 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#1e-07 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#1e-06 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#0.001 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#0.01 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#0.1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

 
## conclusion:
    # For A = 0.5, dlt = 0.5, ds = 0.1, kappa = 2
    # the smallest perturb is 1e-3. 
  
  
#================
# More A values
#================
  
# Aim:
  # want to know different A values form 0.5 on wards
  # the smallest pertub 

for(a in seq(0.5, 1, by = 0.1)){
  
  A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
  
  SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                      A_mat = A_mat_a, dlt_mat = dlt_mat_0.5,
                      sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                      d_vec = D_vec, h = H)
  cat("A:", a, "\n")
  
  for(pert in 10^seq(-9, 1, by = 1)){
    
    I_pert <- I_spar(size = nrow(SG_SG_inv_5$SIGMA_inv), value = pert)
    SG_inv_5_pert <- SG_SG_inv_5$SIGMA_inv + I_pert
    
    if (all(eigen(SG_inv_5_pert, symmetric = T, only.values = T)$value > 0)){
      cat("pert:", pert, "\n")
    }
  }  
}


# A: 0.5 
#pert: 0.001 
#pert: 0.01 
#pert: 0.1 
#pert: 1 
#pert: 10 
#A: 0.6 
#pert: 0.01 
#pert: 0.1 
#pert: 1 
#pert: 10 
#A: 0.7 
#pert: 1 
#pert: 10 
#A: 0.8 
#pert: 1 
#pert: 10 
#A: 0.9 
#A: 1
  









