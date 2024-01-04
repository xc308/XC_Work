#=============
# 4 Jan. 2024
#=============

# Aim:
  # Link Cross-MRF with mixed graphical spatial model framework
  # via model the inverse of uni-variate covariance using CAR, 
  # and cross-variate using chain structure to
  # Simulate the SIGMA and SIGMA_inv


# Method:
  source("034_1D_simu_SG_SG_inv_UniCAR.R")
  # TST10_SpNormPert_SG_SGInv


TST10_SpNormPert_SG_SGInv <- function(p, data, chain = F, A_mat, dlt_mat, 
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
  c_inv <- solve(I_sps - phi * H_adj)
  C11_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = sig2_mat[1, 1])
  
  C11 <- chol2inv(chol(C11_inv)) 
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
    
    Drr_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = sig2_mat[r, r])
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
    
    
    # early perturb SG_inv if not p.d.
    #cat("r", r, "\n")
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
# Simulate
#=========

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


H_adj <- matrix(as.numeric(abs(H) < 0.3), nrow(H), nrow(H))
diag(H_adj) <- 0 # 2nd-order


#---------------------
# phi and p.d.of H_adj
#---------------------

eigen_Hadj <- eigen(H_adj, symmetric = T, only.values = T)$val
1/ max(abs(eigen_Hadj)) # [1] 0.1577773; [1] 0.2249378

phi = 0.15
phi = 0.22


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6_chain)

source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)


SG_SG_inv_6_Wend_bd_chain <- TST10_SpNormPert_SG_SGInv(p = 6, data = hierarchy_data6_chain, 
                                                       chain = T, A_mat = A_01, 
                                                       dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                       phi = phi, H_adj = H_adj, h = H)


# NEW Result:
  # r: 6 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #B cond numb: 4.050954 
    #condition number of C 121.8894 
    #condition number of CDinv 1093.8 
    #No need to perturb. 
    #SG_inv 
    #[1] "Symmetric: Yes"
    #[1] "p.d.: Yes"
    #bdw： 0.35 




# OLD Result: 
  # 3rd-order neighbour adjacency matrix
      # r: 2 
      # bdw： 0.8 
      #No need to perturb.

  # 2nd-order neighbour adjancency matrix
    # r: 2 
    #bdw： 0.85 
    #No need to perturb. 



#==============
# Visualisaton
#==============

#----------
# Wendland
#---------

## SIGMA
plt_Sig(Sigma = SG_SG_inv_6_Wend_bd_chain$SIGMA, p = 6)
plt_Sig(Sigma = log(SG_SG_inv_6_Wend_bd_chain$SIGMA), p = 6)
plt_Sig(Sigma = log(abs(SG_SG_inv_6_Wend_bd_chain$SIGMA)), p = 6)



## SIGMA_inv
plt_Sig(Sigma = as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA_inv), p = 6)
plt_Sig(Sigma = log(as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA_inv)), p = 6)
plt_Sig(Sigma = log(abs(as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA_inv))), p = 6)



range(as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA_inv))
quantile(as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA_inv))
#        0%        25% 
# -4.9092586 -0.2761167 
# 50%        75% 
#  0.0000000  0.0000000 
# 100% 
# 11.7073933 

quantile(as.matrix(SG_SG_inv_6_Wend_bd_chain$SIGMA))
#         0%           25%           50% 
#-1.891355e+00 -1.762919e-01 -9.124712e-17 
#     75%          100% 
#  3.281188e-01  1.604297e+01 









#============
## problems:
#============
  # 1.in SIGMA: we observe each off-diag block
    # all has very high diag element 
    # this is due to Brt use both spN and Reg
    # the Reg use 1 as diag element, 
    # maybe too large

  # 2. in SIGMA_inv, origninal without SpN+Reg Brt
    # will resulting in standard chain structure
    # but now the chain structure is violated, 
    # the graph has no specific structure,


#-------------------
# Potential solution
#-------------------

# 1. check if it's because only 20 locations
    # increase location number

# 2. check if lower Reg number to 1e-5 or 1e-7
    # is the p.d. of different structure still robust?
  
# 3. check if lower Reg number will help
    # preserve the chain-structure
    # and allow bandsparse for easy inference

# 4. check if it's due to uni-varaite CAR inverse back
    # the marginal does not reflect the decreasing marginal
    # relationship

#---------
# Solution
#---------
# see 035_1D_simu_Metern_Chain_6F.R
# key problem:
  # bandsparse shall NOT in the generation process
  # but in the final SIGMA_inv

# TST10 now has been updated to move the bandsparse
  # into the final step;













