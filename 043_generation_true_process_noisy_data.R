#==============
# 22 Jan. 2024
#==============

# Aim:
  # 1. Generate ture processes and corresponding noisy data Z
      # using pre-specified parameters

  # 2. The noisy data is used in neg_logL, and optim for inference
      # such that may compare with the true pre-specified parameters;
      # source(041_inference_functions.R)

  # 3. the optim_pars then used as cokring parameters to predict the 
      # true processes or a collection of specified missings

  # 4. may further cross-validation, compare different scenarios: 
      # a. multivariate vs univariate (A = 0) CV metric, e.g. MSE, RMSE
      # b. UniMatern vs UniCAR, if result similar, then UniCAR is much efficiency

install.packages("dplyr")
library(dplyr)



#====================
# Tri-variate process
#====================

# parameter matrices: 
    # A_mat, dlt_mat, sig2_mat, kappa_mat, tau2_1, tau2_2, tau2_3;
    # total 15 parameters; 

# to obtain realiable inference parameters, each parameters require 20-30 observations
    # so at least require 15*30 = 450 observations

# for tri-variate process, each process at least 150 observations



#=========
# Settings: grid locations, displacement, distance
#=========

ds <- 0.1
s <- seq(-10 + ds/2, 10 - ds/2, by = ds) # num [1:200]
str(s)

H <- outer(s, s, FUN = "-")
H <- t(H)
H[1:10, 1:10]
str(H) # num [1:200, 1:200]

D_vec <- c(abs(H))
str(D_vec) # num [1:40000]
D_vec[1:10]


df <- data.frame(s)
n1 <- n2 <- n3 <- nrow(df)
n <- n1 + n2 + n3



#=======================================
# Hierarchy Data structure: 3 processes
#=======================================

hierarchy_data_3 <- data.frame(
  node_id = c(1, 2, 3, 3 ),
  par_id = c(NA, 1, c(2, 1))
)


#====================
# setting parameters
#====================

# Aim:
# For generating joint SIGMA, the true processes 
# on which desired noisy data are produced


# the True parameters used to generate the hidden true process

# A21 = 0.3; A31 = 0.5; A32 = 0.7
# dlt21 = 0.5; dlt31 = 0.3; dlt32 = 0.7
# sig2 = 1; 
# kappa1 = 1; kappa2 = 1.5; kappa3 = 2
# tau2_1 = 0.45; tau2_2 = 0.5; tau2_3 = 0.55


all_pars_lst <- All_paras(p = 3, data = hierarchy_data_3)
true_pars <- c(0.3, 0.5, 0.7, 0.5, 0.3, 0.7, 1, 1, 1, 1, 1.5, 2)

indx <- 1
for (lst in 1:length(all_pars_lst)){
  for(i in 1:nrow(all_pars_lst[[lst]])){
    for (j in 1:ncol(all_pars_lst[[lst]])){
      if (is.na(all_pars_lst[[lst]][i, j])){
        all_pars_lst[[lst]][i, j] <- true_pars[indx]
        indx <- indx + 1
      }
    }
  }
}

all_pars_lst


#=========================
# Construct true processes
#=========================
# Aim:
  # for generation of noisy data, need true processes and additive error
  # so first use above true parameters generate SIGMA of true processes

# Method:
  # Matern + non-chain
  # TST11

#----------------------------------------
# Generate SIGMA_Y of the true processes
#----------------------------------------
SG_SG_inv_Y <- TST11_SpNormPert_SG_SGInv(p = 3, data = hierarchy_data_3, 
                                         A_mat = all_pars_lst[[1]],
                                         dlt_mat = all_pars_lst[[2]],
                                         sig2_mat = all_pars_lst[[3]],
                                         kappa_mat = all_pars_lst[[4]], 
                                         d_vec = D_vec, h = H)

SG_Y_3 <- SG_SG_inv_Y$SIGMA
Test_sym_pd(SG_Y_3)


#---------------------------------------
# Generate the samples of true processes
#---------------------------------------

# sample_Y ~ MVN(0, SG_Y_3)
# sample_Y = 0 + t(chol(SG_Y_3)) %*% rnorm(n)
sample_Y <- t(chol(SG_Y_3)) %*% rnorm(n)
str(sample_Y) # num [1:600, 1]

p = 3
sample_Y1 <- sample_Y[1:((p - (p - 1))*n1)]
str(sample_Y1) # num [1:200]

sample_Y2 <- sample_Y[((p - (p - 1))*n1 + 1):((p - (p - 2))*n1)]
str(sample_Y2) # num [1:200]

sample_Y3 <- sample_Y[((p - (p - 2))*n1 + 1) : ((p - (p - 3))*n1)]
str(sample_Y3) # num [1:200]


df <- df %>% 
  mutate(smp_Y1 = sample_Y1, 
         smp_Y2 = sample_Y2, 
         smp_Y3 = sample_Y3)

head(df)
#      s      smp_Y1      smp_Y2     smp_Y3
#1 -9.95 0.182526845 -0.07904852  0.5509568
#2 -9.85 0.247444343 -0.08859714  0.1847341
#3 -9.75 0.229683300 -0.22369363 -0.1467838



#=========================
# Construct Noisy data
#=========================

# Noisy data Z = true process + tau2
# tau2_1 = 0.45; tau2_2 = 0.5; tau2_3 = 0.55

Z1 <- sample_Y1 + 0.45*rnorm(n1)
str(Z1) # num [1:200]

Z2 <- sample_Y2 + 0.5*rnorm(n2)
str(Z2) # num [1:200]

Z3 <- sample_Y3 + 0.55*rnorm(n3)
str(Z3) # num [1:200]


df <- df %>% 
  mutate(Z1 = Z1, Z2 = Z2, Z3 = Z3)

head(df)
#     s      smp_Y1      smp_Y2     smp_Y3         Z1          Z2          Z3
#1 -9.95 0.182526845 -0.07904852  0.5509568  0.1154027 -0.05305937  0.25250491
#2 -9.85 0.247444343 -0.08859714  0.1847341  0.4209296 -0.23281324  1.57004844


Z <- c(df$Z1, df$Z2, df$Z3)
str(Z) # num [1:600]


#================================
# Test the 3 types of indx in 042
#================================

str(df[Fit_indx[[1]], ])

#data.frame':	150 obs. of  7 variables:
# $ s     : num  -4.95 -4.85 -4.75 -4.65 -4.55 -4.45 -4.35 -4.25 -4.15 -4.05 ...
# $ smp_Y1: num  -1.106 -0.912 -0.724 -0.507 -0.337 ...
# $ smp_Y2: num  -2.83 -2.42 -1.98 -1.51 -1.19 ...
# $ smp_Y3: num  -4.5 -3.62 -2.88 -2.33 -1.67 ...
# $ Z1    : num  -1.256 -1.17 -0.738 -0.695 -0.263 ...
# $ Z2    : num  -3.18 -2.74 -1.62 -1 -1.4 ...
# $ Z3    : num  -3.73 -3.73 -3.13 -1.91 -1.8 ...


fit_indx <- Fit_indx[[1]]

str(df[fit_indx, ])




#==========
# Inference
#==========
# use the above simulated noisy data to do inference for parameters
  # via optimize the neg_logL function

# see 044














