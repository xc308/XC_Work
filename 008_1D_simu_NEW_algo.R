#=================
# 24 - 26 Oct.2023
#=================

# Aim:
  # 1. write the function to construct SIGMA
  # that have parameters corresponding to theta
  # in the log_like function and can be indexed from 
  # parameter matrix; 

  # 2. modify original make_SIGMA function 
    # to allow p = 2 being included in my own algo
    # does not need to separate

#======
# Algo
#======

build_SG <- function(p, data, A_mat, dlt_mat, sig2_mat, kappa_mat, d_vec, h) {
  
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
    
    if (r == p) return(as.matrix(SIGMA))
  }
}


#==========
# Test algo
#==========

p = 5 

#-----
# data
#-----
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
#List of 4
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
SG5 <- build_SG(p = 5, data = hierarchy_data, 
               A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5,
               sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
               d_vec = D_vec, h = H)

source("Fun_Tst_Sym_Pd.R")
Test_sym_pd(SG5) 
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


plt_Sig(Sigma = SG5, p = 5)


#----------------------
# Experiment on A vals
#----------------------

# Aim:
  # want to know how the change of value A in wave
  # change the corresponding SG

plt_Sig_ini <- function(Sigma, p, ini_vals, dlt) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  
  title <- paste("p =", p, "\nA =", ini_vals, ", delta = ", dlt)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = title)
  #paste("p = ", p, ",", " A = ", ini_vals, ", log")
  
}



for (ini_vals in seq(0.1, 1, by = 0.2)){
  
  A_mat <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = ini_vals)
  
  SG <- build_SG(p = 5, data = hierarchy_data, 
                 A_mat = A_mat, dlt_mat = dlt_mat_0.5,
                 sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                 d_vec = D_vec, h = H)
  
  plt_Sig_ini(Sigma = SG, p = 5, ini_vals = ini_vals, dlt = 0.5)
  #plt_Sig_ini(Sigma = log(SG), p = 5, ini_vals = ini_vals)
  
}




#=============================
# Test on Chain-structure data
#=============================

p = 7 

#-----
# data
#-----
hierarchy_data4 <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 6, 7),
  par_id = c(NA, 1, 2, 3, 4, 5, 6)
)


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data4)
str(all_pars_lst_7)

# List of 4
#$ A_mat    : num [1:7, 1:7] 0 NA 0 0 0 0 0 0 0 NA ...
#$ dlt_mat  : num [1:7, 1:7] 0 NA 0 0 0 0 0 0 0 NA ...
#$ sig2_mat : num [1:7, 1:7] NA 0 0 0 0 0 0 0 NA 0 ...
#$ kappa_mat: num [1:7, 1:7] NA 0 0 0 0 0 0 0 NA 0 ...

source("Fn_set_ini_vals.R")
A_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = 0.1)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


SG7 <- build_SG(p = 7, data = hierarchy_data4, 
                      A_mat = A_mat_0.1, dlt_mat = dlt_mat_0.5,
                      sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                      d_vec = D_vec, h = H)

source("Fun_Tst_Sym_Pd.R")
Test_sym_pd(SG7)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

plt_Sig(Sigma = SG7, p = 7)


plt_Sig_ini_chain <- function(Sigma, p, ini_vals, dlt) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  
  title <- paste("Chain:", "p =", p, "\nA =", ini_vals, ", delta = ", dlt)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = title)
  #paste("p = ", p, ",", " A = ", ini_vals, ", log")
  
}

plt_Sig_ini_chain_log <- function(Sigma, p, ini_vals, dlt) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  
  title <- paste("Chain_log:", "p =", p, "\nA =", ini_vals, ", delta = ", dlt)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = title)
  #paste("p = ", p, ",", " A = ", ini_vals, ", log")
  
}


plt_Sig_ini_chain_inv <- function(Sigma, p, ini_vals, dlt) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  
  title <- paste("Chain_inv:", "p =", p, "\nA =", ini_vals, ", delta = ", dlt)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = title)
  #paste("p = ", p, ",", " A = ", ini_vals, ", log")
  
}


plt_Sig_ini_chain_inv_log <- function(Sigma, p, ini_vals, dlt) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  
  title <- paste("Chain_inv_log:", "p =", p, "\nA =", ini_vals, ", delta = ", dlt)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = title)
  #paste("p = ", p, ",", " A = ", ini_vals, ", log")
  
}



for (ini_vals in seq(0.1, 1, by = 0.2)){
  
  A_mat <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = ini_vals)
  
  SG <- build_SG(p = 7, data = hierarchy_data4, 
                       A_mat = A_mat, dlt_mat = dlt_mat_0.5,
                       sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                       d_vec = D_vec, h = H)
  
  SG_inv <- chol2inv(chol(SG))
  #plt_Sig_ini_chain(Sigma = SG, p = 7, ini_vals = ini_vals, dlt = 0.5)
  plt_Sig_ini_chain_inv_log(Sigma = log(SG_inv), p = 7, ini_vals = ini_vals, dlt = 0.5)
  
}




## Conclusions:
  # Chain structure's SIGMA sparsity is not robust.
  # Chain: SIGMA inverse has robust sparse
  # So, when construct SIGMA, still need to 
  # involve O((np)^3)
  

## To do:
  # Save plts and piece them side by side
  # for discussion







