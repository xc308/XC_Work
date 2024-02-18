#==================================================
# Function for thresholding Cov mat 16-18 Feb.2024
#==================================================

# Aim:
  # initial thresholding = 1e-3
  # check p.d. for the SG_inv after thres
  # if not p.d., then thresholding = 1e-4
  # check p.d. again, 
  # if not 1e-5, and so on, until 1e-15
  
# Fun1: spress_cov: threshold the given cov_mat
# Fun2: check_pd: via the min eigen of the cov_mat after thres, return logic
# Fun3: Thres_tune_cov: 
            # thres_ini: 1e-3
            # cov_mat_thres: the cov after ini thres 1e-3,
            # cov_mat = SG_inv: orignially constructed w/o thres 1e-3


thres_ini <- 1e-3
#thres_lst <- 10^-(0:5)


spress_cov <- function(cov_mat, threshold) {
  
  cov_mat[abs(cov_mat) < threshold] <- 0
  return(cov_mat)
}


check_pd <- function(cov_mat) {
  
  eigen <- eigen(cov_mat, symmetric = T, only.values = T)$val
  min_eign <- min(eigen)
  
  return(min_eign > 0) # a logic 
  
}

# construct a cov_mat with initial threshold
  # cov_mat_thres: cov after ini thres 1e-3
  # cov_mat: original SG_inv in the final run  
            

Thres_tune_cov <- function(thres_ini, cov_mat_thres, cov_mat = SG_inv){
  
  thres <- thres_ini 
  cat("ini thres:", thres, "\n")
  
  while(!check_pd(cov_mat_thres)){
    thres_new <- thres * 0.1
    cat("new thres:", thres_new, "\n")
    
    if (thres_new < 1e-15) {
      break  # Exit loop if threshold becomes 0
      cat("no threshold can make SG_inv p.d.", "\n")
      return(NULL)  # Return NULL if the threshold becomes too small
    }
    
    cov_mat_thres <- spress_cov(cov_mat, threshold = thres_new)
    Tst_sym_pd(cov_mat_thres) # cat new test result
    
    thres <- thres_new # this iteration new become next iter thres to further multiply 0.1
  }
  
  
    return(list(SIGMA = as.matrix(SIGMA), 
           #SIGMA_inv = as.matrix(forceSymmetric(SIGMA_inv))
           SIGMA_inv = as.matrix(cov_mat_thres))
    )
}




#=======
# Test
#=======
#sg_inv_ini <- sg_inv * (abs(sg_inv) > 1e-3)
#SG_SG_inv_thres <- Thres_tune_cov(thres_ini = thres_ini, 
#                                  cov_mat_thres = sg_inv_ini, 
#                                  cov_mat = sg_inv)
# ini thres: 0.001 
# new thres: 1e-04 
#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"




#length(which(SG_SG_inv_thres$SIGMA_inv == 0)) # [1] 822





