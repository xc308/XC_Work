#=============
# 30 Nov. 2023
#=============

# Aim:
  # try different graphical structure to see if wave_v5's
  # p.d. is roubust. 


# Method:
  # use TST6_Pert_build_SG_SGInv
  # different graph structure


# conclusion: 
  # Wave_v5:
    # 1. under chain structure (p = 5)
      # the p.d. of SIGMA_inv is robust
      # without any perturbation.
    # 2. under randomly-selected graph structure (p = 7)
      # p.d. of SIGMA_inv is robust, without any perturb.
        




source("018_Compare_Cond_Numb_waves.R")


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
all_pars_lst_5 <- All_paras(p = 5, data = hierarchy_data3)
all_pars_lst_7 <- All_paras(p = 7, data = hierarchy_data7)


#=========================================================
# Test under all dlt and A combinations, p.d. of SIGMA_inv
#=========================================================
# Method:
# use TST6_SG_SGInv 
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
    
    SG_SG_inv_7 <- TST6_Pert_build_SG_SGInv(p = 5, data = hierarchy_data7, 
                                            A_mat = A_mat_a, dlt_mat = dlt_mat_d, 
                                            sig2_mat = sig2_mat_1, kappa_mat = kappa_mat_2,
                                            d_vec = D_vec, h = H)
    
    Tst_sym_pd(SG_SG_inv_7$SIGMA_inv)
    
  }
}


#============
# Conclusion
#============

# Wave_v5 
  # 1. under chain structure 
  # the p.d. of SIGMA_inv is robust
  # without any perturbation.

  # 2. under 7 field randomly-selected graph structure
  # robust p.d. 















