#==============
# 14 Mar. 2024
#==============

# Aim:
  # Generate 2D noisy data and true processes for 2D inference

# Method:
  # ref: source("046b_generation_True_Y_noisy_data_Tst10c.R")


#=========
# Settings: grid locations, displacement, Adj matrix, phi range
#=========

#---------------
# 2D Grid locations
#---------------
# Aim: for df_2D

ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

coords <- cbind(s, s)
head(coords)
tail(coords)


#----------------
# 2D displacement
#----------------
# Aim: for TST12

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = coords)
str(DSP)
# num [1:200, 1:200, 1:2]

DSP[, , 1] # all Lon DSP
DSP[, , 2] # all Lat DSP


#--------------
# 2D H_adj, phi
#--------------
# Aim: for UniCAR in 2D

Dist <- as.matrix(dist(coords, diag = T, upper = T))
str(Dist) # num [1:200, 1:200]

Nb_radius <- 0.8 # lag-5
H_adj <- matrix(as.numeric(Dist < Nb_radius), nrow(Dist), nrow(Dist))
diag(H_adj) <- 0
H_adj


eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))

phi <- trunc(spec * 100)/100
# 0.1


#-------
# df_2D
#-------

df_2D_TW <- data.frame(lon = s, lat = s)
df_2D_WL <- data.frame(lon = s, lat = s)

n1 <- n2 <- n3 <- n4 <- n5 <- n6 <- nrow(df_2D_TW)

p <- 6
n <- p * n1  # [1] 1200


head(df_2D_TW)
#    lon   lat
# 1 -9.95 -9.95
# 2 -9.85 -9.85



#==========================
# graph structure: 6 fields
#=========================

hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


data_str <- hierarchy_data6


#-----
# pars
#-----
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6_2D <- All_paras_CAR_2D(p = p, data = data_str)


#====================
# setting parameters: for ture processes
#====================

# Aim:
# set the parameters for the ture processes, noisy data, C.V.

# A21 = 0.7
# A31 = A32 = 0.8
# A42 = A43 = 0.65
# A54 = 0.7
# A61 = A63 = A65 = 0.6

# dlt_lon_21 = 0.2
# dlt_lon_31 = dlt_lon_32 = 0.3
# dlt_lon_42 = dlt_lon_43 = 0.4
# dlt_lon_54 = 0.5
# dlt_lon_61 = dlt_lon_63 = dlt_lon_65 = 0.6

# dlt_lat_21 = 0.2
# dlt_lat_31 = dlt_lat_32 = 0.2
# dlt_lat_42 = dlt_lat_43 = 0.2
# dlt_lat_54 = 0.2
# dlt_lat_61 = dlt_lat_63 = dlt_lat_65 = 0.1


# sig2 = 1 all

# tau2_1~6 0.40 0.45 0.50 0.55 0.60 0.65

tau2s <- seq(0.4, 1, by = 0.05)[1:6]

pars_true <- c(0.7, 0.8, 0.8, 0.65, 0.65, 0.5, rep(0.6, 3), 
               0.2, 0.3, 0.3, 0.4, 0.4, 0.5, rep(0.6, 3),
               rep(0.2, 6), rep(0.1, 3), 
               rep(1, 6)) # w/o tau2s


source("Fn_vec2mat.R")
all_pars_lst_CAR_6_true_2D <- vec_2_mat(vector = pars_true, 
                                     all_pars_lst = all_pars_lst_CAR_6_2D)


all_pars_lst_CAR_6_true_2D[[1]]


#=========================
# Construct true processes
#=========================
source("Fn_TST12_SG_SGInv_CAR_2D.R") # for TST


## Tri-Wave
SG_SG_inv_Y_true_TW_2D <- TST12_SG_SGInv_CAR_2D(p = p, data = hierarchy_data6, 
                      A_mat = all_pars_lst_CAR_6_true_2D[[1]],
                      dsp_lon_mat = DSP[, , 1],
                      dsp_lat_mat = DSP[, , 2], 
                      dlt_lon_mat = all_pars_lst_CAR_6_true_2D[[2]],
                      dlt_lat_mat = all_pars_lst_CAR_6_true_2D[[3]],
                      b = "Tri-Wave", 
                      phi = phi, H_adj = H_adj, 
                      sig2_mat = all_pars_lst_CAR_6_true_2D[[4]],
                      reg_ini = 1e-9, thres_ini = 1e-3)



#r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#new thres: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


## Wendland
SG_SG_inv_Y_true_WL_2D <- TST12_SG_SGInv_CAR_2D(p = p, data = hierarchy_data6, 
                      A_mat = all_pars_lst_CAR_6_true_2D[[1]],
                      dsp_lon_mat = DSP[, , 1], 
                      dsp_lat_mat = DSP[, , 2], 
                      dlt_lon_mat = all_pars_lst_CAR_6_true_2D[[2]],
                      dlt_lat_mat = all_pars_lst_CAR_6_true_2D[[3]],
                      b = "Wendland", 
                      phi = phi, H_adj = H_adj, 
                      sig2_mat = all_pars_lst_CAR_6_true_2D[[4]],
                      reg_ini = 1e-9, thres_ini = 1e-3)


# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: No"
#new thres: 1e-05 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#----------------
# True processes
#----------------

## Tri-Wave
SG_Y_true <- SG_SG_inv_Y_true_TW_2D$SIGMA
str(SG_Y_true) # num [1:1200, 1:1200]



## Wendland
SG_Y_true <- SG_SG_inv_Y_true_WL_2D$SIGMA
str(SG_Y_true) # num [1:1200, 1:1200]



#===================================================
# simulate samples of true processes and noisy data
#===================================================
# Y ~ MVN(0, SG_Y_true)
# (Y - 0)(SG_Y_true)^{-0.5} = Z ~ MVN(0, I)
# Y = SG_Y_true^{0.5} * MVN(0, 1)

smp_Y_true <- t(chol(SG_Y_true)) %*% rnorm(nrow(SG_Y_true))
str(smp_Y_true)
# num [1:1200, 1]


sub_vec <- length(smp_Y_true) %/% p # 200

each_smp_Y_true_lst <- lapply(1:p, function(i) smp_Y_true[((i-1)*sub_vec + 1): (i*sub_vec)])
str(each_smp_Y_true_lst)
# List of 6
#$ : num [1:200] 1.04 1.91 2.89 3.51 2.31 ...
#$ : num [1:200] -0.2252 -0.2465 0.0438 2.0744 -0.1348 ...
#$ : num [1:200] 2.91 1.57 2.13 4.08 2.04 ...
#$ : num [1:200] 2.96 2.68 4.49 4.58 4.37 ...
#$ : num [1:200] 5.35 4.24 5.9 7.79 7.55 ...
#$ : num [1:200] 9.4 9.7 12.4 12.8 11.3 ...


#----------
# Tri-Wave
#----------
for (i in 1:length(each_smp_Y_true_lst)){
  #df_2D_TW[paste("smp_Y", i, sep = "")] <- each_smp_Y_true_lst[[i]]
  df_2D_TW[paste("Z", i, sep = "")] <- each_smp_Y_true_lst[[i]] + tau2s[i] * rnorm(n1)
}

head(df_2D_TW)



#----------
# Wendland
#----------
for (i in 1:length(each_smp_Y_true_lst)){
  #df_2D_WL[paste("smp_Y", i, sep = "")] <- each_smp_Y_true_lst[[i]]
  df_2D_WL[paste("Z", i, sep = "")] <- each_smp_Y_true_lst[[i]] + tau2s[i] * rnorm(n1)
}

head(df_2D_WL)
tail(df_2D_WL)

str(df_2D_TW)
# 'data.frame':	200 obs. of  14 variables:


#-----------------------------
# Save above as single obj RDS
#-----------------------------
saveRDS(df_2D_TW, file = "df_2D_TW.rds")
loaded_df <- readRDS("df_2D_TW.rds")
head(loaded_df)
saveRDS(df_2D_WL, file = "df_2D_WL.rds")
loaded_df_WL <- readRDS("df_2D_WL.rds")
head(loaded_df_WL)

file.exists("df_2D_TW.rds")
# [1] TRUE







