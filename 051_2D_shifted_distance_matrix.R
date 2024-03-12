#=============
# 29 Feb.2024
#=============

# Aim:
  # construct shifted distance matrix to be put into Wendland / Tri-Wave function
  
# 2D Wendland function using 2D shifted displacements constructed shifted distance

# Method:
  # 2D displacement: h_lon, h_lat: source("050_2D_coords_displacement.R")
  # 2D shifted displacement: h_lon - dlt1, h_lat - dlt2
  # shifted distance (shd) based on shifted displacements
  # put the shd into Wendland/Tri-Wave function



#=============================
# 2D shifted distance matrix
#=============================
# steps:
  # 1. displacements in 2D, lon, lat, DSP
  # 2. shifts dlt1, dlt2 for lon, lat
  # 3. shifted displacement matrice for (lon - dlt1) and (lat - dlt2)
  # 4. distance using elements in shf displacement matrices lon-dlt1, lat-dlt2
  # 5. put the shifted distance into Wendland function

#m_try <- matrix(rnorm(9), 3, 3) 
#m_try - 1
  
dlt1 <- 0.3
dlt2 <- 0.5

## Displacement (DSP)
dsp_lon <- DSP[, , 1]
dsp_lat <- DSP[, , 2]

dsp_lon_shft <- dsp_lon - dlt1
dsp_lat_shft <- dsp_lat - dlt2
dim(dsp_lat_shft) #20 20


# shifted distance matrix
Dst_shft_mat <- matrix(0, nrow(dsp_lat_shft), ncol(dsp_lat_shft))
for (i in 1:nrow(dsp_lat_shft)) {
  for (j in 1:nrow(dsp_lat_shft)) {
    
    Dst_shft_mat[i, j] <- sqrt((dsp_lon_shft[i, j])^2 + (dsp_lat_shft[i, j])^2)
    
  }
}
str(Dst_shft_mat)
# num [1:20, 1:20]


#========================
# Fn_shifted_distance_mat
#========================

# Aim:
  # Function to generate the matrix of shifted distance using 
      # shifted displacement;
  # for WendLand/Tri-Wave function in 2D as b function

# Arg:
  # dsp_lon_mat: matrix of lon displacement
  # dsp_lat_mat: matrix of lat displacement
  # dlt1: shift on the lon displacement
  # dlt2: shift on the lat displacement

# Value:
  # a matrix of shifted distance 

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

# shifted distance when dlt1 = 0.2, dlt2 = 0.4
Shft_dst_mat_0204 <- Shft_dst_mat(dsp_lon_mat = DSP[, , 1], dsp_lat_mat = DSP[, , 2], 
             dlt1 = 0.2, dlt2 = 0.4)


quantile(Shft_dst_mat_0204)
#       0%       25%       50%       75%      100% 
#0.1414214 0.4472136 0.8602325 1.4212670 3.1144823 

# for R set: can be R = 0.9, such that ~50% of distance will be zero




