#=============
# 15 Apr. 2024
#=============

# Aim:
  # Tensor version of the function Tst_sym_pd 

# Method:
  # 1. matrix computation such as sym verification, eigen value
    # will be done on "cuda"
  # 2. use if else condition then print out result will be
    # move back to "cpu" 


Tst_sym_pd_Tsr <- function(Tsr_mat) {
  
  Tsr_dim <- dim(Tsr_mat)
  
  if (length(Tsr_dim) != 2){
    stop("Input tensor must be a 2D matrix.")
  }
  
  Tsr_mat_t <- torch_t(Tsr_mat)
  
  eq_mat <- torch_eq(Tsr_mat, Tsr_mat_t)
  eq_prod <- torch_prod(eq_mat)
  
  eq_prod_cpu <- eq_prod$cpu() # move to cpu for as.array for if else
  
  if (as.array(eq_prod_cpu) == 1){
    print("sym: Yes")
    
    Tsr_mat_cpu <- Tsr_mat$cpu()
    
    mat_cpu <- as.matrix(as.array(Tsr_mat_cpu))
    
    gpu_mat <- as.gpu.matrix(mat_cpu, device = "cuda")
    
    eig_val <- eigen(gpu_mat, symmetric = T, only.values = T)$val
    
    eig_val_re <- Re(eig_val) # gpu matrix
    
    eig_val_cpu <- as.array(eig_val_re, device = "cpu")
    
    
    if (all(eig_val_cpu > 0)){
      print("p.d.: Yes")
    } else {print ("p.d.: No")}
    
  } else {
    print("sym: No")
  }
 
}




