#=============
# 16 Apr. 2024
#=============

# Aim:
  # GPUmatrix version of Fn_make_DSP_mat.R


# Arg:
  # crds_gpu: coordinates in 2D space stored on "cuda"

# Value:
  # a array stored on "cuda", with
  # [, , 1]: all Lon displacement matrix
  # [, , 2]: all Lat displacement matrix



#==========
# Function
#==========
make_DSP_mat <- function(crds_gpu){
  DSP_gpu <- as.array(0, dim = c(nrow(crds_gpu), nrow(crds_gpu), 2), device = "cuda")
  for (j in 1:nrow(crds_gpu)){
    for (i in 1:nrow(crds_gpu)){
      DSP_gpu[i, j, 1] <- crds_gpu[j, 1] - crds_gpu[i, 1] # lon displacement for all pairs
      DSP_gpu[i, j, 2] <- crds_gpu[j, 2] - crds_gpu[i, 2] # lat displacement
    }
  }
  
  return(DSP_gpu)
}


