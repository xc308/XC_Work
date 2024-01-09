#============
# 9 Jan. 2024
#============

# Aim:
  # Matern + Taper to deal with univariate-large scale


#---------------
# data structure
#---------------

p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


#------------------------------------
# Location, displacements, distance
#------------------------------------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:40]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]; num [1:200, 1:200]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]; num [1:40000]


#===============
# Taper Matern
#===============

source("Fn_Matern_32.R")

C <- Matern_32(Var = 1, Kappa = 2, d_vec = D_vec)
B <- WendLd_32(r = H, R = 0.5, dlt = 0, A = 0.1)
str(B) # num [1:200, 1:200]

CB <- C * B
Tst_sym_pd(CB)

# What I Learned:
  # Wendland cannot be shift


















