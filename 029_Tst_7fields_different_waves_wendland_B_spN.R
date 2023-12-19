#=============
# 18 Dec. 2023
#=============

# Aim:  
  # Follow 027_Tst, want to know if B's functional form
  # is critical to p.d. under spectral normalization of B

# Method:
  # use other versions of wave function

source("Fn_Waves.R")
source("Fn_Wendland_32.R")
#source("021_Test_p.d_SIGMA_inv_Wendland.R")


TST9_SpNormPert_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
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
        B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
        
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
    SG_inv <- Pert_Mat_2(SG_inv)
    
    
    cat("SG_inv", "\n")
    Tst_sym_pd(SG_inv)
    
    
    
    if (r == p) return(
      list(SIGMA = as.matrix(SIGMA), 
           #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
           SIGMA_inv = as.matrix(SG_inv)
      )
    )
  }
}


#=========
# Settings
#=========

#------
# Data
#------
p = 5

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)



p = 7
hierarchy_data7 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5), c(5, 3))
)



#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1 
ds <- 0.05
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:40]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:40, 1:40]

str(as.double(abs(H)))
# num [1:400] ; num [1:1600]


# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) #[1:400]



#=========
# Wave_v4
#=========

#---------
# ds = 0.1
#---------
wv_v4_0505 <- wave_v4(h = H, delta = 0.5, A = 0.5)
kappa(wv_v4_0505) # 153.5588
all_0_row_chck(wv_v4_0505)# NULL


#----------
# ds = 0.05
#----------
wv_v4_0505 <- wave_v4(h = H, delta = 0.5, A = 0.5)
kappa(wv_v4_0505) # [1]400.4099
all_0_row_chck(wv_v4_0505)# NULL



#========
# Wave_v7
#========

#---------
# ds = 0.1
#---------
wv_v7_0505 <- wave_v7(h = H, delta = 0.5, A = 0.5)
kappa(wv_v7_0505) # [1] 1490.949  
all_0_row_chck(wv_v7_0505)# NULL


#----------
# ds = 0.05
#----------
wv_v7_0505 <- wave_v7(h = H, delta = 0.5, A = 0.5)
kappa(wv_v7_0505) # 109.1679
all_0_row_chck(wv_v7_0505) # NULL


#============
# Wendland_32
#============

#---------
# ds = 0.1
#---------
WdLd_0505 <- WendLd_32(r = H, R = 0.5, dlt = 0.5, A = 0.5)
kappa(WdLd_0505) # Inf
all_0_row_chck(WdLd_0505) # [1] 20


#-----------
# ds = 0.05
#-----------
WdLd_0505 <- WendLd_32(r = H, R = 0.5, dlt = 0.5, A = 0.5)
kappa(WdLd_0505) # Inf
all_0_row_chck(WdLd_0505) # [1] 40




#============
# Parameters
#============
source("Fn_para_mat_construct.R")
#all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)


#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================
source("Fn_Tst_sym_pd.R")

#-------
# p = 7
#-------

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


#-----------
# Original B
#-----------
# Method: TST6_Pert_build_SG_SGInv
         # different waves

for (dlt in seq(0.5, 1, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = a)
    
    SG_SG_inv_7 <- TST6_Pert_build_SG_SGInv(p = 7, data = hierarchy_data7, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_7$SIGMA)
    
  }
}

## Results:
# Wave_v4
  # ds = 0.1
  # 1. many combinations not p.d.
  # 2. there are middle fields e.g., r = 5, 6
      # requires perturbation

  # ds = 0.05
  # 1. many NOT p.d. combinations for both SIGMA and SIGMA_inv
      # in which dlt = 1, A = 0.8~1; dlt = 0.8, A = 0.8~1; dlt = 0.9, A=0.8~1
      # not p.d. for both SIGMA and SIGMA_inv
  # 2. middle fields starting from 4 = 4
      # can be NOT p.d. or require large perturbation


# Wave_v7
  # ds = 0.1
  # 1. non of the combination ensure SIGMA and SIGMA_inv p.d. 
  

  # ds = 0.05
  # 1. non of the combination ensure SIGMA and SIGMA_inv p.d.
  # 2. only the first 3 fields are p.d.
    # the rest of fields are not p.d.
    # even the 3rd field still require perturb
  


#-----------------------
# Original B Wendland_32
#-----------------------
# Method: TST7_Pert_build_SG_SGInv
  
for (dlt in seq(0.5, 1, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = a)
    
    SG_SG_inv_7 <- TST7_Pert_build_SG_SGInv(p = 7, data = hierarchy_data7, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_7$SIGMA)
    
  }
}

## Results:
  # Wendland:
    # ds = 0.1:
      # all p.d.

    # ds = 0.05
      # many combinations p.d.



#------------------------
# Spectral Normalized B
#------------------------

for (dlt in seq(0.5, 1, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = a)
    
    SG_SG_inv_7 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_7$SIGMA)
    
  }
}

## Results:
# Wave_v4:
  # ds = 0.1
  # 1. all p.d. no need to perturb
  
  # ds = 0.05
  # 1. all p.d. speed very fast

# Wave_v7:
  # ds = 0.1
  # 1. all p.d. no need to pertub

# Wendland:
  # ds = 0.1
  # all p.d. but with CDrr_invR having much lower 

  # ds = 0.05
  # all p.d., few combinations require perturb 1e-4


#----------------------------
# Individual test for wave_v4
#----------------------------

# dlt = 0.9, A = 0.9

dlt_mat_09 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = 0.9)
A_mat_09 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.9)
A_mat_08 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.8)
A_mat_07 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.7)


SG_SG_inv_7_0709 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                                        A_mat = A_mat_07, dlt_mat = dlt_mat_09, 
                                        sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                        d_vec = D_vec, h = H)



# dlt = 09; A = 09
#r 7 
#No suitable pert found. 
#Min & Max singular value: 1.316385 1.139585e+17 
#Condition number is: 2.302826e+21 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"

# dlt = 09, A =08
# r 7 
#No suitable pert found. 
#Min & Max singular value: 0.271938 1.194123e+16 
#Condition number is: 2.377235e+20 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"

# dlt = 09, A =07
# r 7 
#No suitable pert found. 
#Min & Max singular value: 0.02515938 4.263511e+14 
#Condition number is: 7.805232e+18 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"


# ds <- 0.05
wv_v4_0909 <- wave_v4(h = H, delta = 0.9, A = 0.9)
kappa(wv_v4_0909) # [1]4.082295e+17
all_0_row_chck(wv_v4_0909) # NULL

wv_v4_0908 <- wave_v4(h = H, delta = 0.9, A = 0.8)
kappa(wv_v4_0908) #  4.277516e+18
all_0_row_chck(wv_v4_0908) # NULL

wv_v4_0907 <- wave_v4(h = H, delta = 0.9, A = 0.7)
kappa(wv_v4_0907) #  3.762503e+18
all_0_row_chck(wv_v4_0907) # NULL


##===========
# Conclusion
##===========

# 1. ds = 0.1, 
  # using spectral normalization can ensure
  # all versions of wave function p.d. SIGMA_inv
  # which is significantly improved than using
  # original B

# 2. ds = 0.05
  # wave_v4 and wave_v6: 3/60 combinations SIGMA_inv not p.d.
  # even spectral Norm has lower the condition number of B
  # to very small around 2, but still not p.d.

  # the reason could be due to the very small grid size
  # induce very close numbers, and wave function has both
  # +/- , so close numbers substraction induce cancellation error
  # thereby numerical instabilty

# 3. wendland_32 all p.d for ds = 0.1 and 0.05
  # this might thanks to the function has non-negative values
  # when the condition number of B has been lowered to ideal number
  # it won't suffer from addition cancelation error

  # this also could serve as a suggestion using
  # non-negative function as B can maintian 
  # better numerical stability











