#=============
# 29 Nov. 2023
#=============

# Aim:
  # to understand the connection btw condition number of C and CDinvR
  # and p.d.of SG_inv

  # condition number links with numerical stability
  # so serves as an empirical indication of using b_function
  # that promotes significant diagonal dominance, 
  # c_ii = f(0) >> c_ij = f(xi - xj)


# Method:
  # modify TST5 to add condition number of both C and CDinvR

# Result:
  #



TST6_Pert_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Pert_Mat.R")
  
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
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
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
    
    CDrr_in <- C %*% Drr_inv
    CDR_sym <- forceSymmetric(CDrr_in %*% R)
    
    #CDR_sym <- forceSymmetric(C %*% Drr_inv %*% R)
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
    SG_inv <- Pert_Mat(SG_inv)
    
    
    cat("SG_inv", "\n")
    Tst_sym_pd(SG_inv)
    
    
    
    if (r == p) return(
      list(SIGMA = as.matrix(forceSymmetric(SIGMA)), 
           #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
           SIGMA_inv = as.matrix(forceSymmetric(SG_inv))
      )
    )
  }
}


#=======
# Test
#=======

#------
# data
#------
p = 5

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)


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
# use TST3_SG_SGInv with wave_v6 (slow decay phi = 1/2; same region supprt)

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)



for (dlt in seq(0.5, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST6_Pert_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    
  }
}


#=========
# Discover
#=========

# 1. condition numeber for all runs are very large
  # meaning any tiny errors in the matrix due to flop
  # can be magnified 

# 2. when CDR condition number is larger than e+22
  # majority of SG_inv has NO suitable perturbation


#==========
# Next step
#==========
# Try to use diag-scaling method to lower the condition number
# of CDR


#=========
# Feedback
#=========
# 1. though diag_scaling indeed can lower condition number
  # it will also change the final result of CDR

# 2. still resort to compare differnt wave function's 
  # condition number 
  # then conjecture the wave that promote significant
    # diag dominance has lower condition number

# 3. seperate the 3 matrix multiplication to 3 steps
  # to speed up the process


#---------------------
# Tst output (wave_v4)
#---------------------

r: 4 
condition number of C 2.46381e+16 
condition number of CDinvR 1.919463e+21 
r 4 
No need to perturb. 
SG_inv 
[1] "Symmetric: Yes"
[1] "p.d.: Yes"
r: 5 
condition number of C 5.675521e+16 
condition number of CDinvR 1.91837e+23 
r 5 
No suitable pert found. 
SG_inv 
[1] "Symmetric: Yes"
[1] "p.d.: No"
[1] "Symmetric: Yes"
[1] "p.d.: No"


#-------------------
# Tst output wave_v5
#-------------------

# significantly better condition number
  # most of CN is within (e+15 , e^20)


#=========
# Curious
#=========
# 1. what if decay par in wave is even large, 
  # so decay speed is even faster than wave_v5
  # to promote significant diag dominance
  # will the condition number be even better?


# 2. what if just use an p.d. function
  # what will the condition number be?


#=========
# Measures
#=========

# 1. create a wave_v8: 
  # narrower CP support domain
  # speedy decay phi = 2.5
  # to promote the f(0) >> f(xi-xj) in magnitude

# 2. observe its plots
  
# 3. test the corresponding condition number 
  # and p.d. of final SIGMA_inv

# 4. try Wendland function
  # a pure p.d. function
 

























