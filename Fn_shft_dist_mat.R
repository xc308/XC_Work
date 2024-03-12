#=============
# 12 Mar. 2024
#=============

# Aim:  
  # function to generate the shifted distance matrix using shifted displacemnts
  # for 2D Tri-Wave / Wendland function as b function

# Args:
  # dsp_lon_mat: matrix of lon displacement
  # dsp_lat_mat: matrix of lat displacement
  # dlt1: shift on the lon displacement
  # dlt2: shift on the lat displacement

# Value:
  # a matrix of shifted distance 


#========================
# Fn_shifted_distance_mat
#========================

Shft_dst_mat <- function(dsp_lon_mat, dsp_lat_mat, dlt1, dlt2) {
  
  dsp_lon_shft_mat <- dsp_lon_mat - dlt1
  dsp_lat_shft_mat <- dsp_lat_mat - dlt2
  
  
  Dst_shft_mat <- matrix(0, nrow(dsp_lat_shft_mat), ncol(dsp_lat_shft_mat))
  for (i in 1:nrow(dsp_lat_shft_mat)) {
    for (j in 1:nrow(dsp_lat_shft_mat)) {
      
      Dst_shft_mat[i, j] <- sqrt((dsp_lon_shft_mat[i, j])^2 + (dsp_lat_shft_mat[i, j])^2)
      
    }
  }
  
  return(Dst_shft_mat)
  
}


#======
# Test
#======

# shifted distance when dlt1 = 0.2, dlt2 = 0.4
#Shft_dst_mat_0204 <- Shft_dst_mat(dsp_lon_mat = DSP[, , 1], dsp_lat_mat = DSP[, , 2], 
#                                  dlt1 = 0.2, dlt2 = 0.4)



#quantile(Shft_dst_mat_0204)
#       0%       25%       50%       75%      100% 
#0.1414214 0.4472136 0.8602325 1.4212670 3.1144823 

# for R set: can be R = 0.9, such that ~50% of distance will be zero


