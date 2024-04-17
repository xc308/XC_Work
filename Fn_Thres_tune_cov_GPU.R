#=============
# 8 Apr. 2024
#=============

# Aim:
  # To transform to Function Thres_tune_cov to be suitable for GPU matrix 


source("Fn_Tst_sym_pd_GPU.R")

spress_cov <- function(cov_mat, threshold) {
  
  cov_mat[abs(cov_mat) < threshold] <- 0
  return(cov_mat)
}

#a.gpu[abs(a.gpu) < 0.2] <- 0, gpu matrix can directly substract and assign
#a.gpu

# GPUmatrix
#torch_tensor
#-0.6265  1.5953  0.4874
#0.0000  0.3295  0.7383
#-0.8356 -0.8205  0.5758
#[ CPUDoubleType{3,3} ]


check_pd_gpu <- function(cov_mat) {
  
  eigen_val <- eigen(cov_mat, symmetric = T, only.values = T)$val
  eig_val_real <- Re(eigen_val) # gpumatrix on cuda
  min_eign <- min(eig_val_real)
  
  return(min_eign > 0) # a logical T, F 
  
}
#eig_val <- eigen(a.gpu, symmetric = T, only.values = T)$val
#min(eig_val@gm$real) > 0




#--------------------------------------------
# construct a cov_mat with initial threshold
#--------------------------------------------
  # cov_mat_thres: cov after ini thres 1e-3
  # cov_mat: original SG_inv in the final run  

Thres_tune_cov_gpu <- function(thres_ini, cov_mat_thres, cov_mat = SG_inv){
  
  thres <- thres_ini 
  cat("ini thres:", thres, "\n")
  
  while(!check_pd_gpu(cov_mat_thres)){     # not p.d.
    thres_new <- thres * 0.1
    cat("new thres:", thres_new, "\n")
    
    if (thres_new < 1e-15) {
      cat("no threshold can make SG_inv p.d.", "\n")
      thres_new <- 0 # no value is set to exact zero, no threshold
      #break  # Exit loop if threshold becomes 0
    }
    
    cov_mat_thres <- spress_cov(cov_mat, threshold = thres_new)
    Tst_sym_pd_gpu(cov_mat_thres) # cat new test result
    
    thres <- thres_new # this iteration new become next iter thres to further multiply 0.1
  }
  
  
  return(list(SIGMA_inv_gpu = cov_mat_thres)) # GPUmatrix on cuda
  
}


#===================
# Test on GPU matrix
#===================

#SG_SGinv_WL_6_GPU <- as.gpu.matrix(SG_SG_inv_6_A01dlt05_Wend$SIGMA_inv)
#class(SG_SGinv_WL_6_GPU)
#str(SG_SGinv_WL_6_GPU@gm) # Double [1:1200, 1:1200]


#try_GPU_thres_tune <- Thres_tune_cov_gpu(thres_ini = 1e-3, cov_mat_thres = SG_SGinv_WL_6_GPU, 
#                   cov_mat = SG_SGinv_WL_6_GPU)


#try_GPU_thres_tune$SIGMA_inv

#saveRDS(SG_SG_inv_6_A01dlt05_Wend, file = "SG_SG_inv_6_A01dlt05_Wend.rds")
#readRDS("SG_SG_inv_6_A01dlt05_Wend.rds")

