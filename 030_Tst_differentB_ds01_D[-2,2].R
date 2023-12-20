#==============
# 20 Dec. 2023
#==============

# Aim:
  # want to know if retain grid size as ds = 0.1
  # but enlarge the domain to [-2, 2] 
  # the robustness of p.d. of SIGMA SIGMA_inv under
  # differnt B functions

# Method:
  # TST9_SpNormPert_SG_SGInv with different b functions

# conclusion:
  # 1. wave_v4
    # p = 5, all p.d.
    # p = 7, dlt = 0.4, A = 0.7, 1 SIGMA_inv NOT p.d. 
      # possible reason: similar quatity substraction cancelation error
  
  # 2. wave_v5
    # p = 5, all p.d.
    # p = 7, all p.d. require pertub 1e-6~ 0.01

  # 3. wave_v7
    # p = 5, all p.d.
    # p = 7, all p.d. no perturb required
  
  # 4. wendland
    # p = 5, all p.d. all combinations
    # p = 7, two combinations require perturb 1e-5, 1e-4


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
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
        
        ## spectral normalization of B_rt
        #B_rt <- check_set_SpNorm_pert(B_rt)
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
# settings
#=========

#---------------------------------
# Location, displacement, distance
#---------------------------------

ds <- 0.1
s <- seq(-2 + ds/2, 2 - ds/2, by = ds)
str(s) # num [1:40]

## displacement
H <- outer(s, s, FUN = "-")
H <- t(H)

## distance
D_vec <- as.double(abs(H))
str(D_vec) # num [1:1600]


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



#============
# B functions
#============

#--------
# Wave_v4
#--------
wv_v4_0505 <- wave_v4(h = H, delta = 0.5, A = 0.5)
kappa(wv_v4_0505)  # [1] 431.0894
all_0_row_chck(wv_v4_0505) # NULL


#--------
# Wave_v5
#--------
wv_v5_0505 <- wave_v5(h = H, delta = 0.5, A = 0.5)
kappa(wv_v5_0505)  # [1]  227789.9
all_0_row_chck(wv_v5_0505) # NULL


#--------
# Wave_v7
#--------
wv_v7_0505 <- wave_v7(h = H, delta = 0.5, A = 0.5)
kappa(wv_v7_0505) # [1] 662.7428
all_0_row_chck(wv_v7_0505) # NULL


#----------
# Wendland
#----------
WdLd_0505 <- WendLd_32(r = H, R = 0.5, dlt = 0.5, A = 0.5)
kappa(WdLd_0505) # [1] Inf
all_0_row_chck(WdLd_0505) # [1] 40




#============
# Parameters
#============
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)




#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================

#-----
# p = 5
#-----

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


#------------------------
# Spectral Normalized B
#------------------------

for (dlt in seq(0.5, 1, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST9_SpNormPert_SG_SGInv(p = 5, data = hierarchy_data, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_5$SIGMA)
    
  }
}


## result:
  # 1. wave_v4
    # p = 5, all p.d. for all combinations

  # 2. wave_v5
    # p = 5, all p.d. for all combiations

  # 3. wave_v7
    # p = 5, all p.d. for all combinations

  # 4. wendland
    # p = 5, all p.d. for all combinations


for (dlt in seq(0.1, 0.5, by = 0.1)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.1, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST9_SpNormPert_SG_SGInv(p = 5, data = hierarchy_data, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    cat("SIGMA:", "\n")
    Tst_sym_pd(SG_SG_inv_5$SIGMA)
    
  }
}

# Results:
  # 1. wave_v4
      # p = 5: all p.d. for all combinations 
  
  # 2. wave_v5
      # p = 5: all p.d. for all combinations

  # 3. wave_v7
      # p = 5: all p.d. for all combinations

  # 4. wendland
      # p = 5: all p.d. for all combinations


#-------
# p = 7
#-------

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


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

# results:
  # 1. wave_v4
    # p = 7, all p.d. for all combinations

  # 2. wave_v5
    # p = 7, all p.d. for all combinations

  # 3. wave_v7
    # p = 7, all p.d. for all combinations

  # 4. wendland
    # p = 7, all p.d. for all combinations


for (dlt in seq(0.1, 0.5, by = 0.1)){
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

## results:
  # 1. wave_v4:
    # dlt = 0.4, A = 1, A = 0.7 Not p.d.
    # the rest are all p.d.

  # 2. wave_v5:
    # all p.d. perturb magnitude 0.001 ~ 1e-6

  # 3. wave_v7:
    # all p.d. for all combinations

  # 4. wendland:
    # all p.d. for all combinations, 
    # two combinations require perturb magnitude 1e-5 ~1e-4
      


#----------------
# Individual test
#----------------

wv4_d04a1 <- wave_v4(h = H, delta = 0.4, A = 1)
kappa(wv4_d04a1) # [1] 133687.7
range(wv4_d04a1) # [1] -1  1

wv4_d04a1_spNpert <- check_set_SpNorm_pert(wv4_d04a1)
kappa(wv4_d04a1_spNpert) # [1] 2.727472
range(wv4_d04a1_spNpert)
# inspect wv4_d04a1_spNpert, observe elements with 
# same magnitude but oppostite sign, 
# this could contribute to cancellation error


wv4_d04a03 <- wave_v4(h = H, delta = 0.4, A = 0.3)
kappa(wv4_d04a03) # [1] 133687.7
range(wv4_d04a03) # [1] -0.3  0.3


dlt_mat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = 0.4)
A_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 1)

SG_SG_inv_7 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                                        A_mat = A_mat_1, dlt_mat = dlt_mat_04, 
                                        sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                        d_vec = D_vec, h = H)



# r: 7 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#B cond numb: 2.727472 
#condition number of C 469167 
#condition number of CDinv 6194920 
#condition number of CDinvR 3.423132e+21 
#r 7 
#No suitable pert found. 
#Min & Max singular value: 0.000242941 3.008697e+12 
#Condition number is: 2.868459e+16 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"


##===========================================
# Test the orignal B for different functioins
##===========================================
# Method:
  # TST9_SpNormPert_SG_SGInv but caption out
  # B spectral normalization

  
#----------
# Wave_v4
#----------

#-------
# p = 7
#-------

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


#------------------------
# Original B
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

# results:
# 1. wave_v4
# p = 7, only the first 3 fields are p.d., 
    # from the 4th fields till the 7th field,
    # none of the combinations is p.d.

# 2. wave_v5
# p = 7, only the first 4 fields are p.d.
    # the 5th, 6th 7th filed 
    # none is p.d.

# 3. wave_v7
# p = 7, B condition number for A = 1, dlt = 0.9
  # is 7.3277e+17
  # none of the combination SIGMA, SIGMA_inv
  # is p.d.

# 4. wendland
# p = 7, dlt = 1, A =1, perturb 1
        # dlt = 0.9, A = 1, NOT p.d.
        # dlt = 0.5, A = 0.9, perturb 1.5
        # dlt = 0.5, A = 0.8, perturb 0.1
        # dlt = 0.4, A = 1, NOT p.d.



for (dlt in seq(0.1, 0.5, by = 0.1)){
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



