#=============
# 25 May 2024
#=============
# Aim:
# Want to know one complete construction time for SG,and SG_inv 
# with part of the code offloaded to GPU

# Method:
# 059

# data:
# CAMS Lon_strip_31_new.rds (n = 5910, 6.5GB)

library(Matrix)

#==============
# GPU settings
#==============
#------
# torch
#------
# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)

cat("inter threads:", "\n")
torch_get_num_interop_threads()

cat("intra threads:", "\n")
torch_get_num_threads()


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

#blas_set_num_threads(36) will have error 139 exit

#cat("Updated BLAS threads:", "\n")
#blas_get_num_procs()



#-----------
# GPUmatrix
#------------

#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

system("nvidia-smi")


#===========================
# One complete construction use CAMS df_Lon_Strip_1_new.rds
#===========================

#-----
# df
#-----
df_Lon_Strp_31_Srt <- readRDS("df_Lon_Strip_31_Sort_new.rds")


#----------
# 2D coords
#----------
crds <- cbind(df_Lon_Strp_31_Srt$Lon, df_Lon_Strp_31_Srt$Lat)
#head(crds)


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
# for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
#str(DSP[, , 1]) # num [1:5910, 1:5910]
#str(DSP[, , 2]) # num [1:5910, 1:5910]

#DSP[, , 1][1, 1:160]  #for dlt_lon ini
# [49] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
# [57] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
# [65] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.75
# [73] 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75
#[81] 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75

#DSP[, , 2][1, 1:20] #for dlt_lat ini
#[1]  0.00  0.75  1.50  2.25
#[5]  3.00  3.75  4.50  5.25
#[9]  6.00  6.75  7.50 15.75
#[13] 17.25 18.00 18.75 19.50
#[17] 20.25 21.00 21.75 22.50


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

#spec <- eigen(H_adj, symmetric = T, only.values = T)$val
#max(abs(spec)) # [1] 7.979811

#phi <- 1/max(abs(spec)) # [1] 0.1253745
#phi <- trunc(phi * 100)/100 # [1] 0.12
phi <- 0.12


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



#===========================
# TST12_SG_SGInv_CAR_2D_GPU_RM
#===========================
source("Fn_TST12_SG_SGInv_CAR_2D_GPU.R")

system.time(
  SG_SGinv_CAMS_Lon_Strip_31_GPU <- TST12_SG_SGInv_CAR_2D_GPU(p = 5, data = hierarchy_data_CAMS, 
                                                             A_mat = A_1, 
                                                             dsp_lon_mat = DSP[, , 1], 
                                                             dsp_lat_mat = DSP[, , 2],
                                                             dlt_lon_mat = dlt_lon_02, 
                                                             dlt_lat_mat = dlt_lat_04, 
                                                             b = "Tri-Wave", 
                                                             phi =  phi, H_adj = H_adj,
                                                             sig2_mat = sig2_mat_1, 
                                                             reg_ini = 1e-9, thres_ini = 1e-3)
  
)



