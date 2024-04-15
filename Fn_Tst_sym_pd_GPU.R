#=============
# 15 Apr. 2024
#=============

# Aim:
  # GPUmatrix version of Fn_Tst_sym_pd

Tst_sym_pd_gpu <- function(gpu_mat){
  
  mat_cpu <- as.array(gpu_mat, device = "cpu")
  
  if (isSymmetric(mat_cpu)) {
    print("sym: Yes")
    
    eig_val <- eigen(gpu_mat, symmetric = T, only.values = T)$val
    eig_val_re <- Re(eig_val) # gpu matrix
    
    eig_val_cpu <- as.array(eig_val_re, device = "cpu")
    
    
    if (all(eig_val_cpu > 0)){
      print("p.d.: Yes")
    } else {print ("p.d.: No")}
  
  }else{
    print("sym: No")
  }
}





