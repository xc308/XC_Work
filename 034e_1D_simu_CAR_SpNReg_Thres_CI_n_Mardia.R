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


p = 10

hierarchy_data_full_10 <- data.frame(
  node_id = c(1, 2, rep(3, 2), rep(4,3), rep(5,4), rep(6, 5),
              rep(7, 6), rep(8, 7), rep(9, 8), rep(10, 9)),
  par_id = c(NA, 1, c(1, 2), c(1, 2, 3), 
             c(1, 2, 3, 4), c(1, 2, 3, 4, 5), seq(1, 6), seq(1, 7),
             seq(1, 8), seq(1, 9))
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
abs(H) < 0.2 # lag-3 for str(s) num [1:400]


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
# p = 6
source("Fn_para_mat_construct.R")
all_pars_lst_6 <- All_paras(p = 6, data = hierarchy_data_full)

source("Fn_set_ini_vals.R")
A_01 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 0.1)
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[1]], ini_vals = 1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_6[[3]], ini_vals = 1)


#p = 10
all_pars_lst_10 <- All_paras(p = 10, data = hierarchy_data_full_10)
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[1]], ini_vals = 1)
dlt_05 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_10[[3]], ini_vals = 1)



#------------------
# SIGMA, SIGMA_inv
#------------------
# p =6, n = 200
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full <- TST10b_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data_full, A_mat = A_01,
                                                                    dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                    phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                    thres_ini = 1e-3)





length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full$SIGMA_inv == 0))
# [1] 5183950

length(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full$SIGMA_inv)
# 5760000

#5183950 / 5760000 * 100

#89.09913


# p =10, n = 400
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data_full_10, A_mat = A_1,
                                                                         dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                         phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                         thres_ini = 1e-3)


# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv == 0))
# [1] 14297510


length(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv)
# [1] 16000000

14297510 / 16000000
# [1] 0.8935944


#p = 10, n = 600
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data_full_10, A_mat = A_1,
                                                                            dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                            phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                            thres_ini = 1e-3)





# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv == 0))
# [1] 33398422

length(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv)
#[1] 36000000

33398422 / 36000000
#[1] 0.9277339


# p =10, n = 800
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data_full_10, A_mat = A_1,
                                                                            dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                            phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                            thres_ini = 1e-3)



# r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

length(which(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv == 0))
# [1] 60514476

length(SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10$SIGMA_inv)
# [1] 64000000
60514476 / 64000000
# [1] 0.9455387



# p = 10, n = 1000
SG_SGinv_CAR_SpNReg_thres_TW_a01d05_full_10 <- TST10b_SpNReg_Thres_SG_SGInv(p = 10, data = hierarchy_data_full_10, A_mat = A_1,
                                                                            dlt_mat = dlt_05, sig2_mat = sig2_mat_1, 
                                                                            phi = phi, H_adj = H_adj, h = H, reg_ini = 1e-9,
                                                                            thres_ini = 1e-3)

#r 10 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09
# vector memory limit of 16.0 Gb reached


