#==============
# 20 Dec. 2023
#==============

# # Aim:
  # want to know if retain grid size as ds = 0.1
  # but enlarge the domain to [-10, 10] 
  # the robustness of p.d. of SIGMA SIGMA_inv under
  # differnt B functions


# Method:
# TST9_SpNormPert_SG_SGInv with different b functions


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
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

## displacement
H <- outer(s, s, FUN = "-")
H <- t(H)

## distance
D_vec <- as.double(abs(H))
str(D_vec) # num [1:40000]


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
## p = 5
wv_v4_0505 <- wave_v4(h = H, delta = 0.5, A = 0.5)
kappa(wv_v4_0505)  # [1] 32537.6
all_0_row_chck(wv_v4_0505) # NULL



#--------
# Wave_v5
#--------
wv_v5_0505 <- wave_v5(h = H, delta = 0.5, A = 0.5)
kappa(wv_v5_0505)  # [1]  7.81566e+24
all_0_row_chck(wv_v5_0505) # NULL


#--------
# Wave_v7
#--------
wv_v7_0505 <- wave_v7(h = H, delta = 0.5, A = 0.5)
kappa(wv_v7_0505) # [1]  6121.99
all_0_row_chck(wv_v7_0505) # NULL


#----------
# Wendland
#----------
WdLd_0505 <- WendLd_32(r = H, R = 0.5, dlt = 0.5, A = 0.5)
kappa(WdLd_0505) # [1] Inf
all_0_row_chck(WdLd_0505) # [1] 200




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
    # p = 5: all p.d. for all combinations, no perturb needed 

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
    # p = 7, 5 combinations NOT p.d.
    # the rest are all p.d.
    # pert magnitude 0.001~0.2
    

# 2. wave_v5
  # p = 7, 23 combinations NOT p.d.
    # the rest p.d.
    # pert magnitude 1e-4~0.01

# 3. wave_v7
    # p = 7, all p.d.

# 4. wendland
  # p = 7, all p.d. for all combinations
  # perturb magnitude: 1e-7 ~ 1e-4


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
    # p = 7, 20 combinations NOT p.d.
      # the rest p.d. perturb magnitude 1e-5 ~ 1e-4


# 2. wave_v5:
  # p = 7, 23 combinations NOT p.d.
      # the rest p.d. peturb magnitude 1e-4 ~ 0.01

# 3. wave_v7:
  # 6 combinations NOT p.d.
    # pert magnitude 0.1~0.4

# 4. wendland:
    # all p.d. for all combinations, 
    # peturb magnitude 1e-5~1e-4


#=========
# Original 
#=========


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

## Results:
  # p = 5
  # wave_v4: none of SIGMA_inv is p.d.
  # wave_v5: 12 combinations SIGMA_inv NOT p.d.



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

# Result:
  # p = 7
  # wave_4: none of SIGMA_inv is p.d.; few SIGMA is p.d.
  # wave_5: none of SIGMA_inv p.d; 
  # wendland: 9 SIGMA_inv not p.d; 
            # perturb magnitude 1e-5 ~ 0.3




