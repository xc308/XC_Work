#==============
# 29 Feb. 2024
#==============

# Aim:
  # 1. 2D coords 
  # 2. 2D displacement for 2D Wendland or Tri-Wave function


#===============================
# 2D coords, displacement matrix
#===============================

# Aim:  
# 1. generate coords consisting of Lon, Lat
# 2. generate displacement matrix for Lon, Lat of each pair of locations


#----------
# 2D coords
#----------
ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
crds <- cbind(s, s)
head(crds)
#         s     s
#[1,] -0.95 -0.95
#[2,] -0.85 -0.85
#[3,] -0.75 -0.75
#[4,] -0.65 -0.65
#[5,] -0.55 -0.55

nrow(crds) #[1] 20

#-----------------
# 2D Displacement
#-----------------

DSP <- array(0, dim = c(nrow(crds), nrow(crds), 2))
for (j in 1:nrow(crds)){
  for (i in 1:nrow(crds)){
    DSP[i, j, 1] <- crds[j, 1] - crds[i, 1] # lon displacement for all pairs
    DSP[i, j, 2] <- crds[j, 2] - crds[i, 2] # lat displacement
  }
}

DSP
# , , 1 (Lon displacement)

#[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#[1,]  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8   0.9   1.0   1.1   1.2
#[2,] -0.1  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7   0.8   0.9   1.0   1.1
#[3,] -0.2 -0.1  0.0  0.1  0.2  0.3  0.4  0.5  0.6   0.7   0.8   0.9   1.0
#[4,] -0.3 -0.2 -0.1  0.0  0.1  0.2  0.3  0.4  0.5   0.6   0.7   0.8   0.9


# , , 2 (Lat dispacement)

#[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#[1,]  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8   0.9   1.0   1.1   1.2
#[2,] -0.1  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7   0.8   0.9   1.0   1.1
#[3,] -0.2 -0.1  0.0  0.1  0.2  0.3  0.4  0.5  0.6   0.7   0.8   0.9   1.0
#[4,] -0.3 -0.2 -0.1  0.0  0.1  0.2  0.3  0.4  0.5   0.6   0.7   0.8   0.9


str(DSP)
# num [1:20, 1:20, 1:2]

DSP[, , 1] # ALL Lon displacement
DSP[, , 2] # all Lat displacement



#-------
# Verify: the magnitude of displacement  = Euclidean distance
#-------
sqrt(DSP[, , 1]^2 + DSP[, , 2]^2)


#=====================
# 2D Wendland function
#=====================
# see 049






