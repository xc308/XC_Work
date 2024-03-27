#==============
# 26 Mar. 2024
#==============

# Aim:
  # This is to test desired SG and SG_inv can be constructed by submitting 
    # to Bask. 

# Want to know:
  # where is the final result stored?


#=========
# Settings
#=========

source("Fn_TST12_SG_SGInv_CAR_2D.R")


#----------
# 2D coords
#----------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
crds <- cbind(s, s)


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
  # for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)

# DSP[, , 1] # all Lon DSP
# DSP[, , 2] # all Lat DSP





#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR

DIST <- as.matrix(dist(crds, diag = T, upper = T))

## set Nb_radius = 0.6
Nb_radius <- 0.6

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0


spec <- eigen(H_adj, symmetric = T, only.values = T)$val
phi <- 1/max(abs(spec)) # [1] 0.1344431

phi <- trunc(phi * 100)/100  # [1] 0.13


#--------
# pars
#--------

hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_6 <- All_paras_CAR_2D(p = 6, data =  hierarchy_data6)

source("Fn_set_ini_vals.R")
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[1]], ini_vals = 1)
dlt_lon_02 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[2]], ini_vals = 0.2)
dlt_lat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[3]], ini_vals = 0.4)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_6[[4]], ini_vals = 1)


## Tri-Wave
SG_SGinv_CAR_6_2D_TW <- TST12_SG_SGInv_CAR_2D(p = 6, data = hierarchy_data6, 
                                              A_mat = A_1, 
                                              dsp_lon_mat = DSP[, , 1], 
                                              dsp_lat_mat = DSP[, , 2],
                                              dlt_lon_mat = dlt_lon_02, 
                                              dlt_lat_mat = dlt_lat_04, 
                                              b = "Tri-Wave", phi =  phi, 
                                              H_adj = H_adj,
                                              sig2_mat = sig2_mat_1, 
                                              reg_ini = 1e-9, thres_ini = 1e-3)

#saveRDS(SG_SGinv_CAR_6_2D_TW, file = "SG_SGinv_CAR_6_2D_TW.rds")

#saveRDS(SG_SGinv_CAR_6_2D_TW$SG, file = "SG_SGinv_CAR_6_2D_TW.rds")
str(SG_SGinv_CAR_6_2D_TW$SG_inv)


