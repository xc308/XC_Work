#==============
# 18 Mar. 2024
#==============

# Aim:
  # use the Lon/Lat-sorted df_Res_log_16_sorted to do 2D TST12, inference


# Methods

  # source("057_Data_process_df_Res_log_16_sorted.R")
  # Data: df_Res_log_16_sorted_DAG_V1

  # TST12:source("Fn_TST12_SG_SGInv_CAR_2D.R")

install.packages("dplyr")
library(dplyr)

#==========
# Settings
#==========
df_Res_log_16_sorted_DAG_V1 <- readRDS("df_Res_log_16_sorted_DAG_V1.rds")

df_Res_log_16_sorted_V1_Lon1 <- df_Res_log_16_sorted_DAG_V1 %>%
  filter(Lon_div_indx == 1)

#str(df_Res_log_16_sorted_V1_Lon1)
# 'data.frame':	3793 obs. of  11 variables:
# $ Lon         : num  -179 -179 -179 -179 -178 ...
#$ Lat         : num  71.2 68.2 67.5 66.8 71.2 ...
#$ Year        : num  2016 2016 2016 2016 2016 ...
#$ Z1          : num  0.477 1.115 1.186 0.956 0.338 ...
#$ Z2          : num  -0.3311 -0.0333 0.0737 0.037 -0.3935 ...
#$ Z3          : num  -0.7537 -0.1456 -0.0494 -0.144 -0.8353 ...
#$ Z4          : num  -0.5515 0.0583 0.0629 -0.065 -0.6494 ...
#$ Z5          : num  1.236 0.876 0.949 1.154 1.275 ...
#$ PM25_res    : num  0.455 0.467 0.48 0.46 0.442 ...
#$ Lat_div_indx: num  4 4 4 4 4 4 4 4 4 4 ...
#$ Lon_div_indx: num  1 1 1 1 1 1 1 1 1 1 ...




#----------
# 2D coords
#----------
coords_Lon_1 <- cbind(df_Res_log_16_sorted_V1_Lon1$Lon, df_Res_log_16_sorted_V1_Lon1$Lat)
#str(coords_Lon_1)
# num [1:3793, 1:2]



#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
# for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP_Lon1 <- make_DSP_mat(crds = coords_Lon_1)
#str(DSP_Lon1)
# num [1:3793, 1:3793, 1:2]
#DSP_Lon1[, , 1] # all DSP_Lon
#DSP_Lon1[, , 2] # all DSP_Lat

#str(DSP_Lon1[, , 1]) # num [1:3793, 1:3793]


#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR

DIST <- as.matrix(dist(coords_Lon_1, diag = T, upper = T))
#quantile(DIST)
#       0%       25%       50%       75%      100% 
# 0.00000  16.65270  27.17651  40.83886 126.01116 
#str(DIST) # num [1:3793, 1:3793]

## set Nb_radius = 5
Nb_radius <- 10


H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0

sum(H_adj == 1) / length(H_adj)
# [1] 0.1056757 [1] 0.3286574 [1] 0.07168505 [1] 0.02957242
# so, H_adj has only 30% of non-zero data, very sparse


spec <- eigen(H_adj, symmetric = T, only.values = T)$val
phi <- 1/max(abs(spec)) # [1] 0.002062418  0.00763647

phi <- trunc(phi * 1000)/1000  # [1] 0.002 0.007


#c_inv <- (I_sparse(size = nrow(H_adj), value = 1) - phi * H_adj)
#Tst_sym_pd(c_inv)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"



#--------
# pars
#--------

hierarchy_data_CAMS_V1 <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 5),
  par_id = c(NA, 1, 2, 3, c(1, 4))
)





source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CAMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS_V1)


source("Fn_set_ini_vals.R")
A_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CAMS[[1]], ini_vals = 1)
dlt_lon_02 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CAMS[[2]], ini_vals = 0.2)
dlt_lat_04 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CAMS[[3]], ini_vals = 0.4)
sig2_mat_1 <- Fn_set_ini_vals(pars_mat = all_pars_lst_CAR_2D_CAMS[[4]], ini_vals = 1)


## Tri-Wave
SG_SGinv_CAR_2D_CAMS_TW <- TST12_SG_SGInv_CAR_2D(p = 5, data = hierarchy_data_CAMS_V1, 
                                              A_mat = A_1, 
                                              dsp_lon_mat = DSP_Lon1[, , 1], 
                                              dsp_lat_mat = DSP_Lon1[, , 2],
                                              dlt_lon_mat = dlt_lon_02, 
                                              dlt_lat_mat = dlt_lat_04, 
                                              b = "Tri-Wave", phi = phi, 
                                              H_adj = H_adj,
                                              sig2_mat = sig2_mat_1, 
                                              reg_ini = 1e-9, thres_ini = 1e-3)


# Result on Laptop
# r: 2 
#condition number of C 4.774764e+13 
#condition number of CDinv 6.572813e+12 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 8.296721e+13 
#condition number of CDinv 1.420728e+13 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 1.148135e+14 
#condition number of CDinv 2.024515e+13 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 1.359929e+14 
#condition number of CDinv 2.204811e+13 
#Error: vector memory exhausted (limit reached?)
#In addition: Warning messages:
#  1: In asMethod(object) :
#  sparse->dense coercion: allocating vector of size 1.7 GiB
#2: In asMethod(object) :
#  sparse->dense coercion: allocating vector of size 1.7 GiB


# Result on GPU:
#r: 2 
#condition number of C 4.774772e+13 
#condition number of CDinv 6.572811e+12 
#r 2 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 3 
#condition number of C 8.296721e+13 
#condition number of CDinv 1.420727e+13 
#r 3 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 4 
#condition number of C 1.148135e+14 
#condition number of CDinv 2.024515e+13 
#r 4 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#r: 5 
#condition number of C 1.359934e+14 
#condition number of CDinv 2.204805e+13 
#r 5 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#> 






