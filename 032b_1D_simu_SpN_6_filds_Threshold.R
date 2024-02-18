#============
# 16 Feb.2024
#============

# Aim:
  # produce the 1D-simulation plots for part I of the paper
  # in particular, thresholding and set those value in the precision matrix
    # < threshold to exact zero. 


# Method:
  # source("032_1D_simu_plt_SpN_6_fileds.R")
  # TST9_SpNormPert_SG_SGInv

  # modify to thresholding the final run of SG_inv to exact zero for 
    # values below given threshold

  # threshold ini is set to 1e-3 (quite aggressive)
  # the function sourced "Fn_Thres_tune_cov.R" will automatically check if 
      # the resulting SG_inv after thresholding is p.d. or not
      # if not p.d., it will automatically tune the thresholding value
      # to the one that is 0.1 magnitude smaller
      # and will check the resulting SG_inv after new threshold is p.d. or not
      # until find a largest possible thresholding value that could set to values
      # below this threshold to exact zero while maitain the p.d. of SG_inv

# Values:
  # the function TST9b automatically return the constructed SIGMA and SIGMA_inv
    # that has values below the largest possible threshold being exact zero
    # both SIGMA and SIGMA_inv are p.d.



TST9b_SpNormPert_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat,
                                     d_vec, h, thres_ini = 1e-3) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  source("Fn_Pert_Mat.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpN_Pert.R")# lower kappa of B
  source("Fn_Thres_tune_cov.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  
  C11 <- Matern_32(Var = sig2_mat[1, 1], Kappa = kappa_mat[1, 1], d_vec = d_vec)
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
    SG_inv <- Pert_Mat_2(SG_inv) # larger pertub range
    
    
    cat("SG_inv", "\n")
    Tst_sym_pd(SG_inv)
    
    
    
    if (r == p) {
      
      # threshold SG_inv with thres_ini
      SG_inv_ini <- SG_inv * (abs(SG_inv) > thres_ini)
      
  
      # 1. tune threshold if SG_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIMGA_inv and SIGMA
      Thres_tune_cov(thres_ini = thres_ini, cov_mat_thres = SG_inv_ini,
                     cov_mat = SG_inv) 
      
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
str(s) # num [1:20]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]


# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]


#----------------
# data structure
#----------------

p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6)

A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)


#========================
# Test the TST9b function
#========================
# Aim: 
  # to test the TST9b function, to see if it automatically 
    # tune the thresholding value for SG_inv to ensure both
    # its sparsity and p.d.


SG_SG_inv_6_a01d05_Wend_Thres <- TST9b_SpNormPert_SG_SGInv(p = 6, data = hierarchy_data6, 
                          A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                          sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                          d_vec = D_vec, h = H, thres_ini = 1e-3)



# r 6 
# No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#============================================
# Experiement largest possible ini threshold
#============================================


SG_SG_inv_6_a01d05_Wend <- TST9_SpNormPert_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                    A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                                                    sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                                    d_vec = D_vec, h = H)


sg_inv <- SG_SG_inv_6_a01d05_Wend$SIGMA_inv

min(sg_inv) # [1] -1349.841
min(abs(sg_inv)) # [1] 9.580329e-09

sg_inv_ini <- sg_inv * (abs(sg_inv) > 1e-4)
length(which(sg_inv_ini == 0)) # 822

Tst_sym_pd(sg_inv_ini) # p.d. YES

length(sg_inv_ini) # 120*120 = 14400

length(which(sg_inv_ini == 0)) / length(sg_inv_ini) *100
# 5.708% sparse 






