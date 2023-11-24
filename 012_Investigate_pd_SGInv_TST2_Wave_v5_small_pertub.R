#==================
# 23 - 24 Nov. 2023 
#==================

#Aim:
    # Follow work done in 011_Investigate_pertub_on_wave_v4.R
    # Understand the ini value for dlt
    # and impact on wave function and p.d. of SIGMA_inv

# Discover:
    # need to narrow to compact support of a wave function
    # and faster the correlation decay

# Measures:
    # modified wave_v4 to wave_v5

# Investigate:
    # p.d. of resulting SIGMA_inv using
    # TST2_SG_SGInv function with wave_v5
  
# Conclusion: (24 Nov.)
    # much improved numerical stability of SIGMA_inv
    # dlt = 0.3, A = 1, p.d NO; pert: 1e-7
    # dlt = 0.4, A = 0.9, p.d. NO; pert: 1e-5
    # dlt = 0.7, A = 1, p.d. NO; pert: 1e-4


#------
# data
#------
p = 5

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1 
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) #[1:400]


# H
#      [,1] [,2] [,3] [,4]
# [1,]  0.0  0.1  0.2  0.3
#      [,5] [,6] [,7] [,8]
# [1,]  0.4  0.5  0.6  0.7
#      [,9] [,10] [,11] [,12]
# [1,]  0.8   0.9   1.0   1.1
#      [,13] [,14] [,15] [,16]
# [1,]   1.2   1.3   1.4   1.5
#      [,17] [,18] [,19] [,20]
# [1,]   1.6   1.7   1.8   1.9

## so when dlt = 0.5, 
# tri-wave function include 16 out of 20 values for s1
# include 17 out of 20 for s2
# include 18 out of 20 for s3
# include 19 out of 20 for s4
# include 20 out of 20 from s5 on wards

# this mean the correlation from s with other locations for
# a given row decays very slowly, 
# particularly when A is increasing from 0.5 onwards
# this implies the off-diag part C and R 
# in SIGMA_inv might not be contract,
# so p.d. is affected. 


#=================
# use smaller dlt
#=================

# if dlt set to realtively small value
# consider the compact support condition is 
# /h - dlt/ < 2/dlt/ 
# if dlt = 0.1, 2/dlt/ = 0.2, /h - 0.1/ < 0.2, very compact
# if dlt = 0.2, 2/dlt/ = 0.4, /h - 0.4/ < 0.8, max h = 1.2, 9 out of 20 preserved
# if dlt = 0.3, 2/dlt/ = 0.6, /h - 0.6/ < 1.2, max h = 1.2, 13 out of 20 preserved
# if dlt = 0.4, 2/dlt/ = 0.8, /h - 0.8/ < 1.6.


#--------------------------------------------------------
# Modified wave function (below are all based on wave_v5)
#--------------------------------------------------------
# with a narrower compact support range
source("Fn_Waves.R")
# wave_v5 
# read in conjunction with 010_Visualise_Wave_v5_3D_1D.R



#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)

source("Fn_set_ini_vals.R")
A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


#-----------------
# Test on the algo
#-----------------
SG_SG_inv_5_TST2 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                        A_mat = A_mat_0.5, dlt_mat = dlt_mat_0.5,
                                        sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                        d_vec = D_vec, h = H)


#-------------------
# Test p.d. and symm
#-------------------
Tst_sym_pd(SG_SG_inv_5_TST2$SIGMA)                                                                                                                      d_vec = D_vec, h = H)
Tst_sym_pd(SG_SG_inv_5_TST2$SIGMA_inv)

# [1] "Symmetric: Yes"
# [1] "p.d.: Yes" ! (compare "Investigate_pertub_on_wave_v4.R")


#-----------------
# More values of A
#-----------------

# Aim:
# want to know p.d. of SIGMA_inv for different A values form 0.5 on wards

for(a in seq(0.5, 1, by = 0.1)){
  
  A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
  
  SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                     A_mat = A_mat_a, dlt_mat = dlt_mat_0.5,
                                     sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                     d_vec = D_vec, h = H)
  cat("A:", a, "\n")
  Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
}

# A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

#============
# Conslusion
#============
# when the B function has narrow compact supprt
  # and the function value decays faster as /h - dlt/
  # gets larger
  # this ensures the B function and the off-diag block in
  # SIGMA_inv C and R are contract or convergent
  # theoretical support norm(C) < 1, 
  # manifestation faster correlation decay of B function
  
  # or C is highly likely column linearly indepdent, full column rank
  # BK1 and SIGMA_inv/BK1 p.d. are guaranteed. 


# Now for dlt = 0.5, even A value gets larger beyond 0.5
  # the resulting SIGMA_inv are p.d. for most of A values
  # except only one A = 0.9

# this is highly likely due to numerical issue
# so try pertub


#---------
# Perturb 
#---------
# A = 0.9, dlt = 0.5
  # SIGMA_inv is not p.d.

# try perturb and see what's the smallest value to perturb

for (pert in 10^seq(-9, -1, by = 1)){
  
  A_mat_09 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.9)
  SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                     A_mat = A_mat_09, dlt_mat = dlt_mat_0.5, 
                     sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                     d_vec = D_vec, h = H)
  
  I_pert <- I_spar(size = nrow(SG_SG_inv_5$SIGMA_inv), value = pert)
  SG_inv_pert <- SG_SG_inv_5$SIGMA_inv + I_pert
  
  cat("perturb:", pert, "\n")
  Tst_sym_pd(SG_inv_pert)
  
} 


# perturb: 1e-09 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#perturb: 1e-08 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#perturb: 1e-07 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#perturb: 1e-06 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#perturb: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#perturb: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#perturb: 0.001 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#perturb: 0.01 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#perturb: 0.1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

## Conclusion:
  # only need to perturb 1e-4, SIGMA_inv can be p.d.



#===================================
# More dlt and A values combination
#===================================

#A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
#dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)




for (dlt in seq(0.1, 1, by = 0.2)){
  cat("dlt:", dlt, "\n")
  dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
  
  for (a in seq(0.5, 1, by = 0.1)){
    cat("A:", a, "\n")
    A_mat_a <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = a)
    
    SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                       A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                       sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                       d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_5$SIGMA_inv)
  
  }
}


# dlt: 0.1 
#A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#dlt: 0.3 
#A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#dlt: 0.5 
#A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#dlt: 0.7 
#A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#dlt: 0.9 
#A: 0.5 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.6 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.7 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.8 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 0.9 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#A: 1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#-----------
# conclusion
#-----------
# discover: 
  # 1. dlt = 0.3, A = 1, p.d NO;
      # dlt = 0.4, A = 0.9, p.d. NO;
      # dlt = 0.7, A = 1, p.d. NO

  # 2. indeed due to numerical stability 
      # as 23 Nov. discover dlt = 0.5, A = 0.9, p.d. NO
      # but today(24 Nov.) such combination is p.d YES

  
#----------------------------------
# perturb 3 above p.d. NO situation
#----------------------------------

# dlt = 0.3, A = 1, p.d NO;
# dlt = 0.4, A = 0.9, p.d. NO;
# dlt = 0.7, A = 1, p.d. NO


#A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
#dlt_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.5)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


for (dlt in c(0.3, 0.7)){
  cat("dlt", dlt, "\n")
  
  for (pert in 10^seq(-7,-1, by = 1)){
    cat("pert:", pert, "\n")
    
    dlt_mat_d <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = dlt)
    A_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 1)
    
    SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                       A_mat = A_mat_1, dlt_mat = dlt_mat_d, 
                                       sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                       d_vec = D_vec, h = H)
    
    I_pert <- I_spar(size = nrow(SG_SG_inv_5$SIGMA_inv), value = pert)
    SG_inv_pert <- SG_SG_inv_5$SIGMA_inv + I_pert
    Tst_sym_pd(SG_inv_pert)
  }
}


# dlt 0.3 
#pert: 1e-07 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 1e-06 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.001 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.01 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#dlt 0.7 
#pert: 1e-07 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#pert: 1e-06 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#pert: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#pert: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.001 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.01 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

## conclusion:
  # dlt = 0.3, A = 1, pert = 1e-7
  # dlt = 0.7, A = 1, pert = 1e-4



# dlt = 0.4, A = 0.9, p.d. NO;
for (pert in 10^seq(-7,-1, by = 1)){
  cat("pert:", pert, "\n")
  
  dlt_mat_0.4 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.4)
  A_mat_0.9 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.9)
  
  SG_SG_inv_5 <- TST2_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                     A_mat = A_mat_0.9, dlt_mat = dlt_mat_0.4, 
                                     sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                     d_vec = D_vec, h = H)
  
  I_pert <- I_spar(size = nrow(SG_SG_inv_5$SIGMA_inv), value = pert)
  SG_inv_pert <- SG_SG_inv_5$SIGMA_inv + I_pert
  Tst_sym_pd(SG_inv_pert)
}

## conslusion:
  # pert 1e-5

# pert: 1e-07 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#pert: 1e-06 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#pert: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.001 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.01 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#pert: 0.1 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"



