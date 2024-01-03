#=============
# 1-3 Jan. 2024
#=============

# Aim:
  # Link Cross-MRF with mixed graphical spatial model framework
  # via model the inverse of uni-variate covariance using CAR, 
  # inverse them back for the construction of SIGMA, p*O(n)
  # 


# Methods:
  # 1. construct inverse of SIGMA_11, Drr by uni-variate CAR
      # a. ordered seperation lag in 1D
      # b. set radius as neighbourhood criteria
      # c. turn logic matrix to numerica
      # d. obtain H matrix
      # e. calculate eigen value of H
      # f. set phi to be close to 1/max abs(eigen(H))
      # g. check p.d. of univariate covariance matrix

  # 2. inverse back inverse of univariate cov obtain univariate cov


#===================
# simulation setting
#===================

## Data for Adjacency matrix: separation lag
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20]

DIST <- abs(t(outer(s, s,  FUN = "-")))
abs(outer(s, s, FUN = "-"))

quantile(DIST)
#  0%  25%  50%  75% 100% 
# 0.0  0.3  0.6  1.0  1.9 

## radius
# neighbour radius 0.4, 3-order Neighbour
try <- (DIST < 0.4)
str(try) # logi [1:20, 1:20] TRUE TRUE TRUE TRUE FALSE FALSE
Try <- as.matrix(DIST < 0.4)
str(Try)
all(Try == try) # [1] TRUE

H <- as.numeric(as.matrix(DIST < 0.4))
str(H) #num [1:400] a vector

H <- matrix(H, 20, 20)
diag(H) <- 0
# banded, diag H = 0


#=================
# parameter range
#=================

# eigen value of H
eign_H <- eigen(H, symmetric = T, only.values = T)$val

# $values
#[1]  6.3380470  5.2718922  3.7420321  2.1243407
#[5]  0.6252840  0.0609321 -0.2617542 -0.3501865
#[9] -0.4664053 -0.5227871 -1.0000000 -1.0000000
#[13] -1.0000000 -1.4641555 -1.5575934 -1.7012950
#[17] -1.8189404 -2.0323873 -2.2639353 -2.7230880


#-------------------------------------
# Properties of diagonally zero matrix
#-------------------------------------
# since trace of a matrix is sum of all diag elements
# trace(H) = 0
# meanwhile, trace(H) = sum of all eigen vals (H), 
# SO eigen values of H must be both positive and negative

1/max(abs(eign_H)) # [1] 0.1577773

phi = 0.15


#===============
# SIGMA_INV: CAR
#===============

I_sps <- I_sparse(size = nrow(H), value = 1)
sigma_inv <- solve(I_sps - phi * H)
str(sigma_inv)

Tst_sym_pd(sigma_inv)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

SIGMA_inv <- sigma_inv * I_sparse(size = nrow(H), value = 1)








 


