#============
# 17 Dec.2023
#============

# Aim:
  # apply spetral normalization onto matrix B
  
# Purpose:
  # Since C and R contains repeated multiplication of B 
  # across multiple iterations, 
  # such matrix B exponention would induce numerical
  # instability if there's no constrain on the norm of 
  # matrix B.

  # the changes in the further calculations of BK1, BK2, BK3 
  # that invovling C and R (or B exponentiation) will
  # grow at an ever-increasing rate or exponentially
  # resulting very high condition number of BK1, BK2, BK3
  # and the result of these three blocks would be 
  # unstable and diverge rapidly over iterations

  # such numerical instabilities will induce significant 
  # changes in the calculation of BK1, BK2, BK3, even
  # just in the presence of tiny errors, e.g., 
  # unavoidable rounding error in the flop architecture

  # consequently, the robustness of the positive definiteness 
  # of the SIGMA_inv = BK1, BK2 / BK2, BK4 is violated
  # over multiple iterations in a HM setting. 


# Method:
  # apply spectral_norm function onto B after its construction
  # the eigen value of B will then have its largest magnitude 
  # capped at 1

# Exepected impacts:
  # 1. the condition number of the repeated mulitiplication
    # of B and therefore that of the C, R will be much lower
  # 2. the changes in the calculation of BK1, BK2, BK3 
    # will be more stable in the presence of tiny errors,
    # e.g., rounding error
  # 3. lead to the whole iterative process numerically stable
  # 4. since the normalization is on regression coefficient matrix
    # such normalization would prevent overfitting
    # and promote a more parsimounoious model that 
    # generalize better to new data. 


TST9_SpNormPert_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
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
        #B_rt <- wave_v6(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        #B_rt <- wave_v4(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        
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


#---------
# wave_v5
#---------
# adjust wave_v5 to allow /h-dlt/ <= dlt
# such that /0 - 0.5/ won't be suppress to 0

wv_v5_0505 <- wave_v5(h = H, delta = 0.5, A = 0.5)
kappa(wv_v5_0505)

all_0_row_chck(wv_v5_0505)
# no row has all 0 elements

wv_v5_0505[1, 1] # [1] -0.5
wv_v5_0505[40, 40] # [1] -0.5

View(wave_v5)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)


source("Fn_set_ini_vals.R")
A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================
source("Fn_Tst_sym_pd.R")

#------
# p = 5
#------
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)

for (dlt in seq(0.5, 1, by = 0.2)){
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


# Results:
  # ds = 0.1, p = 5, all p.d. no perturb needed
  # ds = 0.05, p = 5, all p.d. no perturb needed


#-------
# p = 7
#-------

sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


for (dlt in seq(0.5, 1, by = 0.2)){
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


# Results:
# ds = 0.1, p = 7, all p.d. no perturb needed
# ds = 0.05, p = 7, dlt = 0.8, A = 0.7~1 r = 7, NOT p.d. 


#-----------------
# A = 1, dlt = 0.8
#-----------------

all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)


source("Fn_set_ini_vals.R")
A_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 1)
dlt_mat_0.8 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = 0.8)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


SG_SGINV_A1dlt08 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                         A_mat = A_mat_1, dlt_mat = dlt_mat_0.8, 
                         sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                         d_vec = D_vec, h = H)


# r 7 
#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#-----------------
# A = 0.8, dlt = 0.8
#-----------------
A_mat_08 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.8)

SG_SGINV_A08dlt08 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                                             A_mat = A_mat_08, dlt_mat = dlt_mat_0.8, 
                                             sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                             d_vec = D_vec, h = H)

# r 7 
#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#-------------------
# A = 0.6, dlt = 0.8
#-------------------
A_mat_06 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.6)

SG_SGINV_A06dlt08 <- TST9_SpNormPert_SG_SGInv(p = 7, data = hierarchy_data7, 
                                              A_mat = A_mat_06, dlt_mat = dlt_mat_0.8, 
                                              sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                              d_vec = D_vec, h = H)

# r 7 
#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"




for (dlt in seq(0.1, 0.5, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = dlt)
  
  for (a in seq(0.1, 0.5, by = 0.1)){
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






