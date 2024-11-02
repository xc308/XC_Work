#==============
# 30 Oct. 2024
#==============

# Aim: 
  # Understand the exact-zero percentage in SIGMA_inv when 
    # CI among n only (Mardia 1988 case)


# Method:
  # 034b: TST10b


#---------------
# data structure
#---------------

p = 6
hierarchy_data_full <- data.frame(
  node_id = c(1, 2, rep(3, 2), rep(4,3), rep(5,4), rep(6, 5)),
  par_id = c(NA, 1, c(1, 2), c(1, 2, 3), 
             c(1, 2, 3, 4), c(1, 2, 3, 4, 5))
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.05
#s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
#str(s) # num [1:40]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:400]

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:400, 1:400]


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
abs(H) < 0.3 # lag-3 for str(s) num [1:400]


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
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data_full)

source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)


#------------------
# SIGMA, SIGMA_inv
#------------------
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full <- TST10b_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data_full, A_mat = A_01,
                                                                    dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                    phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                    thres_ini = 1e-4)





length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full$SIGMA_inv == 0))
# [1] 5183950

length(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full$SIGMA_inv)
# 5760000

#5183950 / 5760000 * 100

#89.09913
