#==============
# 27 Nov. 2023
#==============

# Aim:
  # test the function Pert_Mat(M) in Fn_Pert_Mat.R

source("Fn_Pert_Mat.R")
source("014_Investigate_p.d._condition_SIGMA_inv_Haville.R")

# Test Aim:
  # dlt = 0.9, A = 0.5, wave_v4
  # dlt = 0.1, A = 0.6, wave_v4
  # dlt = 0.9, A = 1, wave_v4

# Method:
  # use TST4_build_SG_SGInv


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


#-----------
# Parameters
#-----------
source("Fn_para_mat_construct.R")
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data)

source("Fn_set_ini_vals.R")
A_mat_0.5 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.5)
A_mat_0.6 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 0.6)
A_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[1]], ini_vals = 1)

dlt_mat_0.1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.1)
dlt_mat_0.9 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[2]], ini_vals = 0.9)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[3]], ini_vals = 1)
kappa_mat_2 <- Fn_set_ini_vals(pars_mat = all_pars_lst_5[[4]], ini_vals = 2)


SG_SG_inv_5_tst4 <- TST4_build_SG_SGInv(p = 5, data = hierarchy_data, 
                                   A_mat = A_mat_1, dlt_mat = dlt_mat_0.9, 
                                   sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                   d_vec = D_vec, h = H)


SG_inv_v4_dlt09_A05 <- SG_SG_inv_5_tst4$SIGMA_inv
SG_inv_v4_dlt01_A06 <- SG_SG_inv_5_tst4$SIGMA_inv 
SG_inv_v4_dlt09_A1 <- SG_SG_inv_5_tst4$SIGMA_inv



#---------------------------
# Test on Pert_Mat function
#---------------------------

SG_inv_v4_dlt09_A05_pert <- Pert_Mat(SG_inv_v4_dlt09_A05)
# smallest pert: 0.01 

SG_inv_v4_dlt01_A06_pert <- Pert_Mat(SG_inv_v4_dlt01_A06)
# No need to perturb. 

SG_inv_v4_dlt09_A1_pert <- Pert_Mat(SG_inv_v4_dlt09_A1)
# No suitable pert found.


#-----------
# Conclusion
#-----------
# Need to assign Pert_Mat result to a name to store pert result
# for wave_v4, dlt too close to boundary of interal [-1, 1],
# and A is 1, 
