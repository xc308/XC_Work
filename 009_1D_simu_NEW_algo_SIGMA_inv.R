#==============
# 27 Oct. 2023
#==============

# Aim:
  # To extend the current 008_1D_simu_NEW_algo
    # to allow the algo generate the SIGMA_inv
    # simultaneously
    # on the computational order of p * O(n^3)
    # linear order in p, but still cubic order in n


build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Wave_V4.R")
  
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
        B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        
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
    
    Col <- rbind(C, C_rr)
    Row <- rbind(SG, R)
    SIGMA <- cbind(Row, Col)
    
    
    ## SIGMA_Inv
    Drr_inv <- chol2inv(chol(D_rr)) # Schur complement inverse
    
    if (r == 2){
      SG_inv <- chol2inv(chol(SG))
    }
    
    BK1 <- SG_inv %*% (SG + C %*% Drr_inv %*% R)%*% SG_inv
    BK2 <- - SG_inv %*% C %*% Drr_inv
    BK3 <- - Drr_inv %*% R %*% SG_inv
    BK4 <- Drr_inv
    
    SIGMA_inv <- rbind(cbind(BK1, BK2), cbind(BK3, BK4))
    SG_inv <- SIGMA_inv
    
    
    if (r == p) return(
      list(SIGMA = as.matrix(forceSymmetric(SIGMA)), 
           SIGMA_inv = as.matrix(SIGMA_inv))
      )
  }
}


#=======
# Test 
#=======
# Method:
  # the generated SIGMA and SIGMA_inv, 
  # their product shall be I. 


p = 5

#------
# data
#------
hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)
str(all_pars_lst_5)
# List of 4
#$ A_mat    : num [1:5, 1:5] 0 NA NA 0 0 0 0 NA NA 0 ...
#$ dlt_mat  : num [1:5, 1:5] 0 NA NA 0 0 0 0 NA NA 0 ...
#$ sig2_mat : num [1:5, 1:5] NA 0 0 0 0 0 NA 0 0 0 ...
#$ kappa_mat: num [1:5, 1:5] NA 0 0 0 0 0 NA 0 0 0 ...


source("Fn_set_ini_vals.R")
A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.1)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


#------------------------------------
# Location, separation lag, distance
#------------------------------------
ds <- 0.1 
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

# displacements between pairs of points
H <- outer(df$s, df$s, FUN = "-")
H <- t(H)  

# distance
D_vec <- as.double(c(abs(H))) #[1:400]


#-----------------
# Test on the algo
#-----------------
SG_SGInv_5 <- build_SG_SGInv(p = 5, data = hierarchy_data, 
               A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5,
               sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
               d_vec = D_vec, h = H)

str(SG_SGInv_5$SIGMA) # num [1:100, 1:100]
str(SG_SGInv_5$SIGMA_inv) # num [1:100, 1:100]

res_5 <- SG_SGInv_5$SIGMA %*% SG_SGInv_5$SIGMA_inv


#--------------------------------------
# symm & pd test of SIMGA and SIGMA_inv
#-------------------------------------
# symmetric & p.d. are required for further likelihood calculation
  # - det(SIGMA): symmetric is necessary
  # - SIGMA_inv: p.d. is necessary

# so it's better to have them both symmetric and p.d.

Test_sym_pd(SG_SGInv_5$SIGMA_inv)
# [1] "Symmetric: No"
# [1] "p.d.: Yes"

Test_sym_pd(forceSymmetric(SG_SGInv_5$SIGMA_inv))
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


Test_sym_pd(SG_SGInv_5$SIGMA)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"






