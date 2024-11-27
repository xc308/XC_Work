#============
# 20 Feb.2024
#============

# Aim:
  # On top of 034_1D_simu_SG_SG_inv_UniCAR that include a chain optional
  # add parts that will tune the Regularization number
  # and there will be thresholding part at the final SIGMA_inv

  # 6 Nov. 2024 percentage of exact zero for CI of cross-MRF


# Method:
  # source("034_1D_simu_SG_SG_inv_UniCAR.R")
  # TST10_SpNormPert_SG_SGInv

#=========
# Settings
#=========
install.packages("Matrix")
library(Matrix)



TST10b_SpNReg_Thres_SG_SGInv <- function(p, data, A_mat, dlt_mat, sig2_mat,
                                      phi, H_adj, h, reg_ini = 1e-9, 
                                      thres_ini = 1e-3) {
  
  source("Fn_Matern_32.R")
  source("Fn_Check_par_node.R")
  source("Fn_Waves.R")
  source("Fn_Wendland_32.R") # R = 0.5
  #source("Fn_Pert_Mat.R")
  #source("Fn_Pert_Mat_2.R")
  source("Fn_Tst_sym_pd.R")
  source("Fn_check_set_SpNorm_Reg.R") # SpN + tune regularize number
  source("Fn_I_sparse.R")
  source("Fn_Thres_tune_cov.R") # thresholding SIGMA_inv and return SIGMA and SIGMA_inv
  
  # phi: controls the propostion of conditional spatial dependence 
  # sig2_mat: diag matrix containing sigma2 in CAR
  # H: adjacency matrix
  
  I_sps <- I_sparse(size = nrow(H_adj), value = 1)
  c_inv <- I_sps - phi * H_adj
  C11_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = 1/sig2_mat[1, 1])
  
  C11 <- chol2inv(chol(C11_inv)) 
  n <- nrow(C11)
  SIGMA <- C11
  
  reg_num <- reg_ini
  
  restart <- T
  while(restart){
    restart <- F
    SIGMA_inv_pd <- T # initializa label
    
    # reset the SIGMA each while iteration
    SIGMA <- C11
    
    for(r in seq(2, p, by = 1)){
      
      PN <- Check_par_node(Node = r, data = data)
      R <- C <- NULL
      
      cat("r:", r, "\n")
      for(c in seq(1, (r-1), by = 1)){
        
        BT <- NULL
        C_rc <- 0
        for(t in c(PN)){
          B_rt <- wave_v5(h = h, delta = dlt_mat[r, t], A = A_mat[r, t])
          #B_rt <- WendLd_32(r = h, R = 0.5, dlt = dlt_mat[r, t], A = A_mat[r, t])
          
          
          ## spectral normalization of B_rt
          B_rt <- check_set_SpNorm_Reg(B_rt, reg_num = reg_num)
          cat("B cond numb:", kappa(B_rt), "\n")
          
          BT <- rbind(BT, t(B_rt))
          C_rc <- C_rc + B_rt %*% SIGMA[((t-1)*n+1) : (t*n), ((c-1)*n+1): (c*n)]
        }
        
        R <- cbind(R, C_rc)
        
        C_cr <- t(C_rc)
        C <- rbind(C, C_cr)
      }
      
      Drr_inv <- c_inv %*% I_sparse(size = nrow(H_adj), value = 1/sig2_mat[r, r])
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
      #cat("condition number of C", kappa(C), "\n")
      cat("condition number of CDinv", kappa(CDrr_in), "\n")
      cat("condition number of CDinvR", kappa(CDR_sym), "\n")
      
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
      
      # to condition to break the for loop and start all over again
      if (!check_pd(SG_inv)) {
        restart <- T
        SIGMA_inv_pd <- F
        cat("New reg_num needed.", "\n")
        
        if (reg_num < 1) {
          reg_num <- reg_num * 10
          cat("Reg_num updated to:", reg_num, "\n")
        } else {
          reg_num <- reg_num + 0.1
          cat("Reg_num updated to:", reg_num, "\n")
        }
        
        break # break the following of for loop and start from begining again
      }
      
      # early perturb SG_inv if not p.d.
      cat("r", r, "\n")
      cat("SG_inv", "\n")
      Tst_sym_pd(SG_inv)
    }
    
    cat("Final reg_num:", reg_num, "\n")
    
    if (SIGMA_inv_pd){
      
      # Compute SIGMA_inv_ini, SG_inv constructed with the smallest possible reg_num
      SIGMA_inv_ini <- SG_inv * (abs(SG_inv) > thres_ini)
      
      # 1. tune threshold if SIGMA_inv_ini is NOT p.d.,
      # 2. cov_mat construct with new thres
      # 3. check p.d. until cov_mat is p.d. with the updated largest possible thres
      # 4. return the thresholded and p.d. SIGMA_inv and SIGMA
      SG_SG_inv_thres <- Thres_tune_cov(thres_ini = thres_ini, 
                                        cov_mat_thres = SIGMA_inv_ini,
                                        cov_mat = SG_inv) 
      
      return(list(SIGMA = as.matrix(SIGMA),
                  SIGMA_inv = SG_SG_inv_thres$SIGMA_inv))  
        
    }
  }
}


#==================
# Test the function
#==================
# Aim:
  # to test the above constructed function workable
  # for further plts generation so as to formulate better conclusion
  # regarding the cross-MRF, whether it's only able to speed up the generation
  # or is able to further sparse the SIGMA_inv

#---------------
# data structure
#---------------

p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


p = 10 
hierarchy_data10 <- data.frame(
  node_id = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10),
  par_id = c(NA, 1, 2, c(2, 3), c(1,4), c(1, 5), 3, 6, 6, 6)
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.05
#s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
#str(s) # num [1:40]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:400]

s <- seq(-15 + ds/2, 15 - ds/2, by = ds)
str(s) # num [1:600]

s <- seq(-20 + ds/2, 20 - ds/2, by = ds)
str(s) # num [1:800]

s <- seq(-25 + ds/2, 25 - ds/2, by = ds)
str(s) # num [1:1000]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:40, 1:40]; num [1:400, 1:400]; num [1:600, 1:600]; num [1:1000, 1:1000]


# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:1600]; num [1:160000]


#--------------------------------------
# Sepration lag based adjacency matrix
#--------------------------------------

# separation lag
abs(H)

# radius for definition of neighbourhood
#abs(H) < 0.4 # 3-order
abs(H) < 0.2 # lag-3 for str(s) num [1:400]; num [1:1000, 1:1000]

#as.numeric(abs(H) < 0.4) # vector
#H_adj <- matrix(as.numeric(abs(H) < 0.4), nrow(H), nrow(H))
#diag(H_adj) <- 0


H_adj <- matrix(as.numeric(abs(H) < 0.2), nrow(H), nrow(H))
diag(H_adj) <- 0


#---------------------
# phi and p.d.of H_adj
#---------------------

eigen_Hadj <- eigen(H_adj, symmetric = T, only.values = T)$val
1/ max(abs(eigen_Hadj)) # [1] 0.07236703 ; [1] 0.1319365

#phi = 0.07236703
phi = 0.1319365
phi = trunc(phi * 100)/100
#phi = 0.15
#phi = 0.14
#[1] 0.13 


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data6)
all_pars_lst_10 <- All_paras(p = 10, data = hierarchy_data10)



source("Fn_set_ini_vals.R")
## p = 6
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)


## p = 10
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[1]], ini_vals = 1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[3]], ini_vals = 1)


#------------------
# SIGMA, SIGMA_inv
#------------------
SG_SGinv_CAR_SpNReg_thres_WL_a01d05 <- TST10b_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, A_mat = A_01,
                             dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                             phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                             thres_ini = 1e-3)

# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 

## New result for str(s): 400
# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"




SG_SGinv_CAR_SpNReg_thres_TW_a01d05 <- TST10b_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, A_mat = A_01,
                                                                    dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                    phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                    thres_ini = 1e-3)


# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001

## New result for str(s): 400
# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 

SG_SGinv_CAR_SpNReg_thres_TW_a1d05 <- TST10b_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, A_mat = A_1,
                                                                    dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                    phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                    thres_ini = 1e-3)




# p = 10, n = 400, Tri-wave

SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data10, A_mat = A_1,
                                                                   dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                   phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                   thres_ini = 1e-3)



# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 



# p =10, n = 600, Tri-wave
SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data10, A_mat = A_1,
                                                                       dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                       phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                       thres_ini = 1e-3)



# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


#p=10, n = 800
SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data10, A_mat = A_1,
                                                                       dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                       phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                       thres_ini = 1e-3)


# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 





# p = 10, n = 1000, Tri-wave

SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data10, A_mat = A_1,
                                                                       dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                       phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                       thres_ini = 1e-3)



# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 


#========
# Plots
#========
plt_Sig(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA, p = 6)
plt_Sig(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA_inv, p = 6)
plt_Sig(log(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA_inv), p = 6)
plt_Sig(log(abs(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA_inv)), p = 6)



#===================================
# exact zero percentage in SIGMA_inv
#===================================


## Wendland
length(which(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA_inv == 0))
# [1] 32914
# [1] 5446246 for str(s) 400

length(SG_SGinv_CAR_SpNReg_thres_WL_a01d05$SIGMA_inv)
# [1] 57600
# [1] 5760000


32914 / 57600 * 100
# [1] 57.14236%

5446246 / 5760000 * 100 # [1] 94.55288


## Tri-Wave
length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05$SIGMA_inv == 0))
# [1] 24572
# [1] 5355018


24572 / 57600 * 100
# 42.65972%

5355018/ 5760000 * 100 # 92.96906


# A = 1
length(which(SG_SGinv_CAR_SpNReg_thres_TW_a1d05$SIGMA_inv == 0))

#---------------------------------
# p = 10, n = 400, A = 1, Tri-wave
#---------------------------------
length(which(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv == 0))
# [1] 15551776

length(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv)
# [1] 16000000

15551776 / 16000000
# [1] 0.971986

#---------------------------------
# p = 10, n = 600, A = 1, Tri-wave
#---------------------------------

length(which(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv == 0))
# [1] 35311700

length(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv)
# [1] 36000000

35311700 / 36000000
#[1] 0.9808806


#---------------------------------
# p = 10, n = 800, A = 1, Tri-wave
#---------------------------------

length(which(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv == 0))
# [1] 63077944

length(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv)
# [1] 64000000

63077944 / 64000000

#---------------------------------
# p = 10, n = 1000, A = 1, Tri-wave
#---------------------------------
length(which(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv == 0))
# [1] 98845436

length(SG_SGinv_CAR_SpNReg_thres_TW_a1d05_p10$SIGMA_inv)
# [1] 100000000

98845436 / 100000000
# [1] 0.9884544


#--------------------------------
# compare with Uni Matern in 032c
#--------------------------------

## Wendland
length(which(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv == 0))
# [1] 30766

length(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv)
# [1] 57600

30766 / 57600 * 100
# [1] 53.41319


## Tri-Wave
length(which(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA_inv == 0))
# [1] 21250
21250 / 57600 * 100
# 36.89236%


## p = 6, n=400 see 032c


##-----------
# Conclusion
##-----------
# 1. there is indeed a different percentage of exact zero
  # in the SG_inv generated via uniCAR and uniMatern

# 2. we significantly feel the generation speed is super fast












