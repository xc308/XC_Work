#============
# 29 Jan.2024
#============

# Aim:
  # Generate true process and corresponding noisy data for 1D simu 
    # for 6 fields with graph structure hierarchy_data6

  # with such df, we can use them to do neg_logL, optim and cokrig

library(dplyr)

#=============================
# 6-field processes parameters 
#=============================

# parameters in 1D
  # A_mat, dlt_mat, sig2_mat, tau2_1, ..., tau2_6
  
# if use hierarchy_data6 graph structure
  # total pars = 9 + 9 + 6 + 6 = 30

# to obtain a reliable inference, each pars require at least 15-20 obs
  # so obs at least 30 * 20 = 600

# now we have 6 field, each field at least 100 obs

## base on these anaylsis, we decide below grid size and total grids


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
str(H_adj) # num [1:200, 1:200]

eig <- eigen(H_adj, symmetric = T, only.values = T)$val
spec <- 1/max(abs(eig))  #  0.1412107
phi <- trunc(spec * 100) / 100 # trunc only preserve the integer part
# [1] 0.14


df <- data.frame(s)
n1 <- n2 <- n3 <- n4 <- n5 <- n6 <- nrow(df) # 200

p <- 6
n <- p * n1   #  1200 


#================
# graph structure: 6 fields
#================
  
p = 6
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

# A21 = 0.6
# A31 = A32 = 0.5
# A42 = A43 = 0.4
# A54 = 0.5
# A61 = A63 = A65 = 0.6

# dlt21 = 0.2
# dlt31 = dlt32 = 0.3
# dlt42 = dlt43 = 0.4
# dlt54 = 0.5
# dlt61 = dlt63 = dlt65 = 0.6

# sig2 = 1 all

# tau2_1~6 0.40 0.45 0.50 0.55 0.60 0.65
tau2s <- seq(0.4, 1, by = 0.05)[1:6]
# [1] 0.40 0.45 0.50 0.55 0.60 0.65


# A, delt, sig2
pars_true <- c(0.6, 0.5, 0.5, 0.4, 0.4, 0.5, rep(0.6, 3), 
  0.2, 0.3, 0.3, 0.4, 0.4, 0.5, rep(0.6, 3),
  rep(1, 6))
str(pars_true) # num [1:24]


all_pars_lst_CAR_6_true <- vec_2_mat(vector = pars_true, all_pars_lst = all_pars_lst_CAR_6)


## Wendland
SG_SG_inv_Y_true <- TST10_SpNormPert_SG_SGInv(p = p, data = data_str, 
                                            chain = F, A_mat = all_pars_lst_CAR_6_true[[1]], 
                                            dlt_mat = all_pars_lst_CAR_6_true[[2]], 
                                            sig2_mat = all_pars_lst_CAR_6_true[[3]], 
                                            phi = phi, H_adj = H_adj, h = H)


SG_Y_true <- SG_SG_inv_Y_true$SIGMA
str(SG_Y_true) #num [1:1200, 1:1200]


#===================================================
# simulate samples of true processes and noisy data
#===================================================
set.seed(0303)

# Y ~ MVN (0, SG_Y_true)
smp_Y_true <- t(chol(SG_Y_true)) %*% rnorm(nrow(SG_Y_true))
str(smp_Y_true) # num [1:1200, 1]

sub_vec <- length(smp_Y_true) %/% p

each_smp_Y_T_lst <- lapply(1:p, function(i) smp_Y_true[((i-1)*sub_vec + 1): (i*sub_vec)])
str(each_smp_Y_T_lst) # List of 6


for (i in 1:p){
  #df[paste("smp_Y", i, sep = "")] <- each_smp_Y_T_lst[[i]]
  df[paste("Z",i, sep = "")] <- each_smp_Y_T_lst[[i]] + tau2s[i]*rnorm(n1)
  
}

head(df)







