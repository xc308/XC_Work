#=============
# 12 Mar. 2024
#=============

# Aim:
  # Function to generate displacement matrix from 2D coords

# Args:
  # crds: coordinates in 2D space

# Value:
  # a array, with
  # [, , 1]: all Lon displacement matrix
  # [, , 2]: all Lat displacement matrix


#==========
# Function
#==========
make_DSP_mat <- function(crds){
  DSP <- array(0, dim = c(nrow(crds), nrow(crds), 2))
  for (j in 1:nrow(crds)){
    for (i in 1:nrow(crds)){
      DSP[i, j, 1] <- crds[j, 1] - crds[i, 1] # lon displacement for all pairs
      DSP[i, j, 2] <- crds[j, 2] - crds[i, 2] # lat displacement
    }
  }
  
  return(DSP)
}


#======
# Test
#======
#make_DSP_mat(crds = crds)





