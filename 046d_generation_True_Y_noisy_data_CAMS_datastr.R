#=============
# 19 Apr. 2024
#=============

# Aim:
  # Generate true processes and noisy data from data with CAMS RA structure
    # for GPU inference testing

# Method:
  # ref:source("046c_generation_True_Y_noisy_data_TST12_2D.R")


#=========
# Settings: grid locations, displacement, Adj matrix, phi 
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

#Nb_radius <- 0.8 # lag-5
Nb_radius <- 0.6 # lag-4
H_adj <- matrix(as.numeric(Dist < Nb_radius), nrow(Dist), nrow(Dist))
diag(H_adj) <- 0
H_adj


eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))

phi <- trunc(spec * 100)/100
# [1] 0.12


#-------
# df_2D
#-------

df_2D_TW <- data.frame(lon = s, lat = s)
#df_2D_WL <- data.frame(lon = s, lat = s)

n1 <- n2 <- n3 <- n4 <- n5 <- nrow(df_2D_TW)

p <- 5
n <- p * n1  # [1] 1200


head(df_2D_TW)


#==========================
# graph structure: CAMS RA
#=========================

hierarchy_data_CAMS <- data.frame(
  node_id = c(1, 2, 3, 4,  5, 5),
  par_id = c(NA, 1, 2, 3, c(4, 1))
)

p = 5

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)


#====================
# setting parameters: for ture processes
#====================

# Aim:
# set the parameters for the ture processes, noisy data, C.V.

# A21 = 0.7
# A32 = 0.8
# A43 = 0.65
# A54 = 0.7
# A51 = 0.6

# dlt_lon_21 = 0.2
# dlt_lon_32 = 0.3
# dlt_lon_43 = 0.4
# dlt_lon_54 = 0.5
# dlt_lon_51 = 0.6

# dlt_lat_21 = 0.3
# dlt_lat_32 = 0.3
# dlt_lat_43 = 0.3
# dlt_lat_54 = 0.3
# dlt_lat_51 = 0.3


# sig2 = 0.5 all

# tau2_1~5 0.40 0.45 0.50 0.55 0.60 

tau2s <- seq(0.4, 1, by = 0.05)[1:5]

pars_true <- c(0.7, 0.8, 0.65, 0.5, 0.6, 
               0.2, 0.3, 0.4, 0.5, 0.6,
               rep(0.3, 5), 
               rep(0.5, 5)) # w/o tau2s


source("Fn_vec2mat.R")
all_pars_lst_CAR_CAMS_true_2D <- vec_2_mat(vector = pars_true, 
                                        all_pars_lst = all_pars_lst_CAR_2D_CMS)




#=========================
# Construct true processes
#=========================
source("Fn_TST12_SG_SGInv_CAR_2D.R")

## Tri-Wave
SG_SG_inv_Y_true_TW_2D <- TST12_SG_SGInv_CAR_2D(p = p, data = hierarchy_data_CAMS, 
                                                A_mat = all_pars_lst_CAR_CAMS_true_2D[[1]],
                                                dsp_lon_mat = DSP[, , 1],
                                                dsp_lat_mat = DSP[, , 2], 
                                                dlt_lon_mat = all_pars_lst_CAR_CAMS_true_2D[[2]],
                                                dlt_lat_mat = all_pars_lst_CAR_CAMS_true_2D[[3]],
                                                b = "Tri-Wave", 
                                                phi = phi, H_adj = H_adj, 
                                                sig2_mat = all_pars_lst_CAR_CAMS_true_2D[[4]],
                                                reg_ini = 1e-9, thres_ini = 1e-3)



# SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#----------------
# True processes
#----------------

## Tri-Wave
SG_Y_true <- SG_SG_inv_Y_true_TW_2D$SIGMA
str(SG_Y_true) # num [1:1000, 1:1000]
Tst_sym_pd(SG_Y_true)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


smp_Y_true <- t(chol(SG_Y_true)) %*% rnorm(nrow(SG_Y_true))
str(smp_Y_true)
# num [1:1000, 1]


sub_vec <- length(smp_Y_true) %/% p # 200

each_smp_Y_true_lst <- lapply(1:p, function(i) smp_Y_true[((i-1)*sub_vec + 1): (i*sub_vec)])
str(each_smp_Y_true_lst)
# List of 5


#----------
# Tri-Wave
#----------
for (i in 1:length(each_smp_Y_true_lst)){
  #df_2D_TW[paste("smp_Y", i, sep = "")] <- each_smp_Y_true_lst[[i]]
  df_2D_TW[paste("Z", i, sep = "")] <- each_smp_Y_true_lst[[i]] + tau2s[i] * rnorm(n1)
}

head(df_2D_TW)

df_2D_TW_CAMS <- df_2D_TW

saveRDS(df_2D_TW_CAMS, file = "df_2D_TW_CAMS.rds")
#readRDS("df_2D_TW_CAMS.rds")





