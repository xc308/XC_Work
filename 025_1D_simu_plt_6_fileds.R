#============
# 6 Dec. 2023
#============

# Aim:
    # generate SIGMA SIGMA_inv for 7 fields
    # plot the corresponding graph


# Method:
    # TST7_Pert_build_SG_SGInv (Wendland 32)
    # TST8_Pert_build_SG_SGInv(wave_v5)


TST8_Pert_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Pert_Mat.R")
  source("Fn_Wendland_32.R")
  source("Fn_Tst_sym_pd.R")
  
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
        #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
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
    #cat("condition number of C", kappa(C), "\n")
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


#===========================
# Test a new data structure
#===========================

#----------------
# data structure
#----------------

p <- 7
TST8_Pert_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Pert_Mat.R")
  source("Fn_Wendland_32.R")
  source("Fn_Tst_sym_pd.R")
  
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
        #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
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
    #cat("condition number of C", kappa(C), "\n")
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


#===========================
# Test a new data structure
#===========================

#----------------
# data structure
#----------------

TST8_Pert_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Pert_Mat.R")
  source("Fn_Wendland_32.R")
  source("Fn_Tst_sym_pd.R")
  
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
        #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
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
    #cat("condition number of C", kappa(C), "\n")
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


#===========================
# Test a new data structure
#===========================

#----------------
# data structure
#----------------

TST8_Pert_build_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Pert_Mat.R")
  source("Fn_Wendland_32.R")
  source("Fn_Tst_sym_pd.R")
  
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
        #B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
        B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
        
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
    #cat("condition number of C", kappa(C), "\n")
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


#===========================
# Test a new data structure
#===========================

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

source("Fn_set_ini_vals.R")
A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
#A_mat_0.2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.2)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[4]], ini_vals = 2)


SG_SG_inv_6_wave_v5 <- TST8_Pert_build_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                                                sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                                d_vec = D_vec, h = H)

# r 6 
#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


SG_SG_inv_6_wave_v5_02 <- TST8_Pert_build_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                A_mat = A_mat_0.2, dlt_mat = dlt_mat_0.5, 
                                                sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                                d_vec = D_vec, h = H)




#-------------
# Plot wave_v5
#-------------
Tst_sym_pd(SG_SG_inv_6_wave_v5$SIGMA)
par(mfrow = c(1,1))

image.path
# [1] "./Results/"
jpeg(paste0(image.path, "SG_6_W5.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = SG_SG_inv_6_wave_v5$SIGMA, p = 6)
dev.off()


jpeg(paste0(image.path, "SG_6_inv_W5.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = SG_SG_inv_6_wave_v5$SIGMA_inv, p = 6)
dev.off()


jpeg(paste0(image.path, "SG_Inv_6_W5_log.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = log(abs(SG_SG_inv_6_wave_v5$SIGMA_inv)), p = 6)
dev.off()


#==========
# Wendland
#==========

SG_SG_inv_6_Wend <- TST8_Pert_build_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5, 
                                                sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                                d_vec = D_vec, h = H)

# r 6 
#No need to perturb. 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

jpeg(paste0(image.path, "SG_6_Wendland.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = SG_SG_inv_6_Wend$SIGMA, p = 6)
dev.off()


jpeg(paste0(image.path, "SG_6_inv_Wendland.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = SG_SG_inv_6_Wend$SIGMA_inv, p = 6)
dev.off()


jpeg(paste0(image.path, "SG_Inv_6_Wendland_log.jpeg"),
     width = 8, height = 7, units = "in",
     res = 300)
plt_Sig(Sigma = log(abs(SG_SG_inv_6_Wend$SIGMA_inv)), p = 6)
dev.off()


