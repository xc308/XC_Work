#==============
# 11 Mar. 2024
#==============

# Aim:
  # Wendland and Tri-Wave function in 2D 

# Method:
  # shifted distance: 051_2D_shifted_distance_matrix.R


#=====================
# 2D Wendland function
#=====================

# Aim:
# construct b function using Wendland in 2D space

# Arg:
# shft_dst: a matrix of distance caluculated based on shifted displacements in lon, lat
# R: support radius beyond which, function value drops to zero
# A: amplitude


WendLd_32_2D <- function(shft_dst, R, A){
  
  R <- round(quantile(shft_dst)[3], 1) # find the 50% quantile of shft_dst
  ifelse(shft_dst <= R, A*(1 - shft_dst/R)^4 * (1 + 4 * shft_dst/R), 0)
  
}

B_WL_2D <- WendLd_32_2D(shft_dst = Dst_shft_mat, R = 0.9, A = 1)

str(B_WL_2D) 
num [1:20, 1:20]


B_WL_2D_dlt0204 <- WendLd_32_2D(shft_dst = Shft_dst_mat_0204, R = 0.9, A = 1)
str(B_WL_2D_dlt0204)
quantile(B_WL_2D_dlt0204)
#          0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 1.838579e-05 1.913940e-01 8.219221e-01


#----------------------------------------
# Modify: R is decided by 50% of shft_dst
#----------------------------------------

WendLd_32_2D <- function(shft_dst, A){
  
  R <- round(quantile(shft_dst)[3], 1) # find the 50% quantile of shft_dst
  ifelse(shft_dst <= R, A*(1 - shft_dst/R)^4 * (1 + 4 * shft_dst/R), 0)
  
}


B_WL_2D_dlt0204 <- WendLd_32_2D(shft_dst = Shft_dst_mat_0204, A = 1)
str(B_WL_2D_dlt0204)
quantile(B_WL_2D_dlt0204)
#    0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 1.838579e-05 1.913940e-01 8.219221e-01




#=====================
# 2D Tri-Wave function
#======================

TriWave_2D <- function(shft_dst, R, A) {
  
  A * (1 - 2 * (shft_dst / R)^2) * (shft_dst <= R) 
}


TriWave_2D_dlt0204 <- TriWave_2D(shft_dst = Shft_dst_mat_0204, R = 0.9, A = 1)


## R set 0.9 see 051_2D_shfited_distance_matrix.R
q_Shft_dst_mat_0204 <- quantile(Shft_dst_mat_0204)
# 50%: 0.8602325 
str(q_Shft_dst_mat_0204)
# Named num [1:5] 0.141 0.447 0.86 1.421 3.114

round(q_Shft_dst_mat_0204[3], 1) # 0.9 

round(30.68,1) # [1] 30.7

str(round(quantile(Shft_dst_mat_0204)[3], 1)) 
# Named num 0.9
  # - attr(*, "names")= chr "50%"


str(TriWave_2D_dlt0204)
quantile(TriWave_2D_dlt0204)
#        0%        25%        50%        75%       100% 
#-0.8271605  0.0000000  0.0000000  0.5061728  0.9506173 


#-----------------------------------------
# Modify: R is decided by 50% of shft_dst
#-----------------------------------------
TriWave_2D <- function(shft_dst, A) {
  
  R <- round(quantile(shft_dst)[3], 1) # find the 50% quantile of shft_dst
  A * (1 - 2 * (shft_dst / R)^2) * (shft_dst <= R) 
}


TriWave_2D_dlt0204 <- TriWave_2D(shft_dst = Shft_dst_mat_0204, A = 1)

quantile(TriWave_2D_dlt0204)
#0%        25%        50%        75%       100% 
#-0.8271605  0.0000000  0.0000000  0.5061728  0.9506173 


#========
# SpN+Reg
#========

source("Fn_check_set_SpNorm_Reg.R")
B_TriWave_SpNReg <- check_set_SpNorm_Reg(TriWave_2D_dlt0204, reg_num = 1e-9)

str(B_TriWave_SpNReg)
# Formal class 'dgeMatrix' [package "Matrix"] with 4 slots

quantile(B_TriWave_SpNReg@x)
#        0%        25%        50%        75%       100% 
#-0.1666665  0.0000000  0.0000000  0.1019899  0.1915421 



B_WL_SpNReg <- check_set_SpNorm_Reg(B_WL_2D_dlt0204, reg_num = 1e-9)

str(B_WL_SpNReg)
quantile(B_WL_SpNReg@x)
#0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 5.352062e-06 5.571436e-02 2.392597e-01 






#-----
# V2: depreach
#-----

TriWave_2D_V2 <- function(shft_dst, dlt1, dlt2, A) {
  
  dlt_dst <- sqrt(dlt1^2 + dlt2^2)
  A * (1 - 2 * (shft_dst / dlt_dst)^2) * (shft_dst <= dlt_dst) 
}

TriWave_2D_dlt0204_V2 <- TriWave_2D_V2(shft_dst = Shft_dst_mat_0204, dlt1 = 0.2, dlt2 = 0.4, A = 1)

sqrt(dlt1^2 + dlt2^2) # [1] 0.5830952
quantile(Shft_dst_mat_0204)


quantile(TriWave_2D_dlt0204_V2)
#0%  25%  50%  75% 100% 
#-1.0  0.0  0.0  0.0  0.8  


# Conclusion:
# Tri-Wave version 1 the better. 

