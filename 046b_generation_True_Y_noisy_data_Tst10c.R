#=============
# 26 Feb. 2024
#=============

# Aim:
  # generate true processes and noisy data for 1D-simu using graph in Fig12 in paper
  # using TST10b_SpNReg_Thres_SG_SGInv in 034b

  # with such a df, one can use for neg-logL, optm, cokrig


library(dplyr)

#=============================
# 6-field processes parameters 
#=============================

# 1D pars:
  # A_mat, dlt_mat, sig2_mat, tau2_1, ..., tau2_6

# if use hierarchy_data_6
  # total pars = 9 + 9 + 6 + 6 = 30

# for realiable inference results, each par require at least 20-30 obs
  # so total obs at least 30 * 20 = 600

# 6 fields, each field has at least 100 obs

# base on these we decide below grid size, domain


#=========
# Settings: grid locations, displacement, Adj matrix, phi range
#=========
ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds)
str(s) # num [1:200]

H <- t(outer(s, s, FUN = "-"))
str(H) # num [1:200, 1:200]

H_adj <- matrix(as.numeric(abs(H) < 0.4), nrow(H), nrow(H))
diag(H_adj) <- 0

eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1 / max(abs(eig)) # [1] 0.06666361
phi <- trunc(spec * 100)/100 # trunc only preserves integer part 6
# [1] 0.14

df <- data.frame(s)
n1 <- n2 <- n3 <- n4 <- n5 <- n6 <- nrow(df)

p <- 6
n <- p * n1  # [1] 2400

#================
# graph structure: 6 fields
#================

hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


data_str <- hierarchy_data6

source("Fn_para_mat_construct.R")
all_pars_lst_CAR_6 <- All_paras_CAR(p = p, data = data_str)


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

# dlt21 = 0.2
# dlt31 = dlt32 = 0.3
# dlt42 = dlt43 = 0.4
# dlt54 = 0.5
# dlt61 = dlt63 = dlt65 = 0.6

# sig2 = 1 all

# tau2_1~6 0.40 0.45 0.50 0.55 0.60 0.65

tau2s <- seq(0.4, 1, by = 0.05)[1:6]

pars_true <- c(0.7, 0.8, 0.8, 0.65, 0.65, 0.7, rep(0.6, 3), 
               0.2, 0.3, 0.3, 0.4, 0.4, 0.5, rep(0.6, 3),
               rep(1, 6))


source("Fn_vec2mat.R")
all_pars_lst_CAR_6_true <- vec_2_mat(vector = pars_true, all_pars_lst = all_pars_lst_CAR_6)


#=========================
# Construct true processes
#=========================
source("034c_1D_simu_CAR_SpNReg_Thres_b_choice.R")

SG_SG_inv_Y_true_WL <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                             A_mat = all_pars_lst_CAR_6_true[[1]],
                             dlt_mat = all_pars_lst_CAR_6_true[[2]], 
                             sig2_mat = all_pars_lst_CAR_6_true[[3]],
                              phi = phi, H_adj, h = H, 
                              b = "Wendland",
                              reg_ini = 1e-9, thres_ini = 1e-3) 

# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001 
#new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


SG_SG_inv_Y_true_TW <- TST10c_SpNReg_Thres_SG_SGInv(p = 6, data = hierarchy_data6, 
                                                    A_mat = all_pars_lst_CAR_6_true[[1]],
                                                    dlt_mat = all_pars_lst_CAR_6_true[[2]], 
                                                    sig2_mat = all_pars_lst_CAR_6_true[[3]],
                                                    phi = phi, H_adj, h = H, 
                                                    b = "Tri-Wave",
                                                    reg_ini = 1e-9, thres_ini = 1e-3) 

# r 6 
#SG_inv 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"
#Final reg_num: 1e-09 
#ini thres: 0.001

## SG_inv using Tri-Wave is sparser than those using Wendland


## True processes
SG_Y_true <- SG_SG_inv_Y_true_TW$SIGMA
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
#head(smp_Y_true)
#tail(smp_Y_true)


sub_vec <- length(smp_Y_true) %/% p # 200

each_smp_Y_true_lst <- lapply(1:p, function(i) smp_Y_true[((i-1)*sub_vec + 1): (i*sub_vec)])

str(each_smp_Y_true_lst)
# List of 6
#$ : num [1:200] 0.662 0.469 -1.249 -1.598 -1.655 ...
#$ : num [1:200] -0.707 -0.847 -1.171 0.85 -1.176 ...
#$ : num [1:200] 1.1 -0.122 -0.982 -1.452 -1.33 ...
#$ : num [1:200] 0.572 0.273 -2.267 -1.417 -0.764 ...
#$ : num [1:200] -0.652 2.278 0.97 0.832 2.844 ...
#$ : num [1:200] -1.069922 -3.293922 -0.02897 -0.000217 0.979236 ...


for (i in 1:length(each_smp_Y_true_lst)){
  #df[paste("smp_Y", i, sep = "")] <- each_smp_Y_true_lst[[i]]
  df[paste("Z", i, sep = "")] <- each_smp_Y_true_lst[[i]] + tau2s[i] * rnorm(n1)
}

head(df)


##--------------------------------------
# Alternative way to generate noisy data
##--------------------------------------

Zi <- rnorm(n1, mean = each_smp_Y_true_lst[[i]], sd = sqrt(tau2s[i]))












