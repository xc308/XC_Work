#======================
# 6 Dec. , 19 Dec. 2023
#======================

# Aim:
  # 1. when ds <- 0.01, location 200, H matrix 40000
  # 2. when ds <- 0.05, location 40, H matrix 1600
  # 3. 027_Tst reveal some combination test together 
      # will fail to be p.d. but when test individually
      # they will all be p.d., so would like to tst
      # those parameter combinations

  # Test robustness of SIGMA_inv




source("021_Test_p.d_SIGMA_inv_Wendland.R")



#=======
# Test
#=======

#------
# data
#------
p = 5

hierarchy_data3 <- data.frame(
  node_id = c(1, 2, 3, 4, 5),
  par_id = c(NA, 1, 2, 3, 4)
)



p = 7
hierarchy_data7 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5), c(5, 3))
)

p = 7
hierarchy_data7new <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5), c(6, 5, 3))
)



#------------------------------------
# Location, displacements, distance
#------------------------------------
#ds <- 0.1
ds <- 0.05
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:40]

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:40, 1:40]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:1600]

#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data3)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)
all_pars_lst_7_new <- All_paras(p = 7, data = hierarchy_data7new)



#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================
# Method:
# use TST6_SG_SGInv 
# but for ds <- 0.05
# str(D_vec)

#----------
# 5 Fields
#----------

source("Fn_set_ini_vals.R")
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)



for (dlt in seq(0.5, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST6_Pert_build_SG_SGInv(p = 5, data = hierarchy_data3, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    
  }
}


#--------
# Results
#--------

# most of r = 3 SIGMA_inv are p.d.
# but when r = 4, 5, p.d. not robust
# This is due to now location is 400 by 400,
# 5 fields, 
# robustness is not good. 


#----------------
# Individual test
#----------------

# A = 1, dlt = 0.9

dlt_mat_09 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.9)
A_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 1)


SG_SG_inv_5_d09a1 <- TST6_Pert_build_SG_SGInv(p = 5, data = hierarchy_data3, 
                                        A_mat = A_mat_1, dlt_mat = dlt_mat_09,
                                        sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                        d_vec = D_vec, h = H)



## Results:
  # r 5 
  # No suitable pert found. 
  # Min & Max singular value: 18719.19 1.90768e+24 
  # Condition number is: 2.448404e+32 
  # SG_inv 
  # [1] "Symmetric: Yes"
  # [1] "p.d.: No"


#---------
# 7 fields
#---------
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


for (dlt in seq(0.1, 0.5, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[2]], ini_vals = dlt)
  
  for (a in seq(0.1, 0.5, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST6_Pert_build_SG_SGInv(p = 5, data = hierarchy_data3, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
    
  }
}







#----------
# Results
#----------
# for small dlt and A, 
# all SIGMA_inv p.d.
# for 5 fields



#===========================
# Use Wendland function test
#===========================
# Method:
  # TST7_Pert_build_SG_SGInv (Wendland)
  # ds <- 0.05
  # H : 40
  # D_vec: 1600
  # total size of matrix for 5 fields
  # 40 * 5 by 40 *5 =. 200 by 200, total 40000


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data3)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)


source("Fn_set_ini_vals.R")
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7[[4]], ini_vals = 2)


for (dlt in seq(0.5, 1, by = 0.2)){
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
    
  }
}


#--------
# Results
#--------
# 5 fields, 40 * 40 locations, total matrix size 40000 
# Wendland all p.d.



# 7 fields, 40 * 40 locations, total 280 * 280
# for dlt 0.5-1, A 0.5 - 1
# dlt 0.9, A 0.9, 1 r = 7 NOT p.d.


#=======================================
# Test 7 fields with new graph structure
#=======================================
# 
#-----------------
# Data structure
#-----------------
p = 7
hierarchy_data7new <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7, 7),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5), c(6, 5, 3))
)



#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1
#ds <- 0.05
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:40]

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:40, 1:40]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:1600]

#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_7_new <- All_paras(p = 7, data = hierarchy_data7new)



#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================
# Method:
# use TST6_SG_SGInv (wave_v5)
# str(D_vec)

source("Fn_set_ini_vals.R")
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[4]], ini_vals = 2)


for (dlt in seq(0.5, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[1]], ini_vals = a)
    
    SG_SG_inv_7 <- TST6_Pert_build_SG_SGInv(p = 7, data = hierarchy_data7new, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    
  }
}


#--------
# Results
#--------
# wave_v5
# dlt = 0.5 -1, A = 0.5 - 1
# Not all p.d.for r = 7 , esp. dlt A approximate 1


#------------------
# Method: Wendland
#------------------
# Wendland
# New 7 fields structure
# ds <- 0.1

for (dlt in seq(0.5, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_7_new[[1]], ini_vals = a)
    
    SG_SG_inv_7 <- TST7_Pert_build_SG_SGInv(p = 7, data = hierarchy_data7new, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    
  }
}

#-------
# Result
#-------
# all p.d.






TST7_Pert_build_SG_SGInv









