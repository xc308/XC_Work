#=============
# 4 Jan. 2024
#=============

# Aim:
  # use the Matern uni-variate to see if chain-structure
  # is preserved and bandspase is possible even
  # after SpN + Reg (1)

  # if this is the case, then problem lies in the univariate model

# Method:
  # TST11_SpNormPert_SG_SGInv: from TST9 add chain = F, and Bandsparse
  # for chain-structure


# Conclusion:
  # the problem lies NOT in the SpN + Reg
  # but in the bandsparse position,
  # shall NOT in the generation process, 
  # but in the final SIGMA_inv

# base on the conclusion, we may modify below TST11
  # to set the bandsparse only in the final return


TST11_SpNormPert_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat,
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

#=========================
# Test on Chain structure
#=========================

#---------------
# data structure
#---------------

p = 6
hierarchy_data6_chain <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 6),
  par_id = c(NA, 1, 2, 3, 4, 5)
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]; num [1:200, 1:200]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]; num [1:40000]


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6_chain)

source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)


#------
# Test
#------
SG_SG_inv_matern_chain <- TST11_SpNormPert_SG_SGInv(p = 6, data = hierarchy_data6_chain,
                          A_mat = A_01, dlt_mat = dlt_05, 
                          sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2, 
                          d_vec = D_vec, h = H, 
                          chain = T)

# New Results:
  # r: 6 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #condition number of C 316088.5 
    #condition number of CDinv 37883.65 
    #condition number of CDinvR 4.906088e+18 
    #No need to perturb. 
    #SG_inv 
    #[1] "Symmetric: Yes"
    #[1] "p.d.: Yes"
    #bdw： 0.3 



# Old Results:
  # r 2 
    #bdw： 0.7 
    #No need to perturb. 
    #SG_inv 
    #[1] "Symmetric: Yes"
    #[1] "p.d.: Yes"
     
    #r 3 
    #bdw： 0.85 
    #No need to perturb. 

  # r 6 
    #No suitable pert found. 
    #Min & Max singular value: 1.059322 1.188263e+17 
    #Condition number is: 2.142002e+21 
    #SG_inv 
    #[1] "Symmetric: Yes"
    #[1] "p.d.: No"
        
 
# each iteration bandspase will induce the non-p.d.
# for final SIGMA_inv


# now remove each iteration bandsparse
SG_SG_inv_matern <- TST11_SpNormPert_SG_SGInv(p = 6, data = hierarchy_data6_chain,
                                                    A_mat = A_01, dlt_mat = dlt_05, 
                                                    sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2, 
                                                    d_vec = D_vec, h = H, 
                                                    chain = F)

# Result:
  # without bandsparse in the middel field
  # all p.d. and the speed is very fast

# Base on this result, want to know:
  # 1. whether there are chain structure in SIGMA_inv plot?
  # 2. whether final SG_inv will have bandsparse for easy inference
  # 3. think: your theory indicates final bandspare or 
      # in the process ?
  

#-------------------------------
# Check if band wd is applicable to final SIGMA_inv
#-------------------------------

SIGMA_inv_6_Chain <- SG_SG_inv_matern$SIGMA_inv

SIGMA_inv_Bd_Sps <- NULL
for(bdw in seq(0.2, 0.95, by = 0.05)){
  SG_inv_bd <- band(SIGMA_inv_6_Chain, -bdw * nrow(SIGMA_inv_6_Chain), 
                    bdw * ncol(SIGMA_inv_6_Chain))
  
  if (all(eigen(SG_inv_bd, symmetric = T, only.values = T)$val > 0)){
    SG_inv_6_chain_sps <- as(SG_inv_bd, "sparseMatrix")
    
    cat("bdw：", bdw, "\n")
    SIGMA_inv_Bd_Sps <- SG_inv_6_chain_sps
    break
  }
}
# bdw： 0.3
# bdw： 0.5


#--------------------
# Visualise SIGMA_inv
#--------------------
# Now want to visualise the bandsparse SIGMA_inv_Bd_Sps
plt_Sig(Sigma = as.matrix(SIGMA_inv_Bd_Sps), p = 6)
plt_Sig(Sigma = log(abs(as.matrix(SIGMA_inv_Bd_Sps))), p = 6)
plt_Sig(Sigma = log(as.matrix(SIGMA_inv_Bd_Sps)), p = 6)
# for asymmetry off-diag block

## result:
  # very good! chain-structured, with bdw = 0.3
  # beyond were all set to exact zero by bandsparse
  # meanwhile p.d. is still maintained. 


# what I learned: 
  # 1. the key problem lies NOT in the Regularize number
    # but in the bandsparse order;

  # 2. bandsparse is for the final SIGMA_inv
    # and cannot be used in the middle of generation process
    # otherwise will sabotash the p.d. of the final SIGMA_inv

  

#--------------------
# Visualise SIGMA
#--------------------

plt_Sig(Sigma = as.matrix(SG_SG_inv_matern$SIGMA), p = 6)
plt_Sig(Sigma = log(abs(as.matrix(SIGMA_inv_Bd_Sps))), p = 6)
plt_Sig(Sigma = log(as.matrix(SG_SG_inv_matern$SIGMA)), p = 6)

# result:
  # no sparse in SIGMA as expected






