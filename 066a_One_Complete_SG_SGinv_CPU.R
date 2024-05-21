#============
# 21 May 2024
#============

# Aim:
  # want to know how long will it take to construct one complete SG, SG_inv
  # only on CPU, no offload to GPU


# Method:
  # 054: TST12_SG_SGInv_CAR_2D

# Data:
  # 


library(Matrix)



#-------------------------------
# check BLAS and OPENBLAS info: only 36 works
#-------------------------------
#install.packages("RhpcBLASctl")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(RhpcBLASctl)


#cat("Check Current BLAS Library", "\n")
#sessionInfo()

cat("Check the current number of BLAS threads", "\n")
blas_get_num_procs()


#==============
# Pre-settings
#==============

#-----
# df
#-----
df_Lon_Strp_1_Srt <- readRDS("df_Lon_Strip_1_Sort_new.rds")


#----------
# 2D coords
#----------
crds <- cbind(df_Lon_Strp_1_Srt$Lon, df_Lon_Strp_1_Srt$Lat)
#head(crds)


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
# for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
#str(DSP[, , 1]) # num [1:3793, 1:3793]
#str(DSP[, , 2]) # num [1:3793, 1:3793]

#DSP[, , 1][1, 1:20]  #for dlt_lon ini
# [1] 0.00 0.00 0.00 0.00 0.75
#[6] 0.75 0.75 0.75 0.75 1.50
#[11] 1.50 1.50 1.50 2.25 2.25
#[16] 2.25 3.00 3.00 3.00 3.75

#DSP[, , 2][1, 1:20] #for dlt_lat ini
# [1]  0.00  0.75  1.50  4.50
#[5] -0.75  0.00  0.75  1.50
#[9]  4.50 -0.75  0.00  0.75
#[13]  4.50 -0.75  0.00  0.75
#[17] -0.75  0.00  0.75 -1.50




#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR

DIST <- as.matrix(dist(crds, diag = T, upper = T))

Nb_radius <- 1.5 # degree 

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0

#dim(H_adj) # 3793 3793
#length(which(H_adj != 0)) # 27772

spec <- eigen(H_adj, symmetric = T, only.values = T)$val
#max(abs(spec)) # [1] 7.976106

phi <- 1/max(abs(spec)) # [1] 0.1253745
phi <- trunc(phi * 100)/100 # [1] 0.12


#--------
# pars
#--------

hierarchy_data_CAMS <- data.frame(
  node_id = c(1, 2, 3, 4,  5, 5),
  par_id = c(NA, 1, 2, 3, c(4, 1))
)

p = 5

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)

source("Fn_set_ini_vals.R")
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[1]], ini_vals = 1)
dlt_lon_02 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[2]], ini_vals = 0.2)
dlt_lat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[3]], ini_vals = 0.4)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CMS[[4]], ini_vals = 1)



#==========================
# TST12 function (non-GPU)
#==========================

source("Fn_TST12_SG_SGInv_CAR_2D.R")


system.time(
SG_SGinv_CAR_Lon_Strip_1_CPU <- TST12_SG_SGInv_CAR_2D(p = 5, data = hierarchy_data_CAMS, 
                                             A_mat = A_1, 
                                             dsp_lon_mat = DSP[, , 1], 
                                             dsp_lat_mat = DSP[, , 2],
                                             dlt_lon_mat = dlt_lon_02, 
                                             dlt_lat_mat = dlt_lat_04, 
                                             b = "Tri-Wave", 
                                             phi =phi, H_adj = H_adj,
                                             sig2_mat = sig2_mat_1, 
                                             reg_ini = 1e-9, thres_ini = 1e-3)
)



