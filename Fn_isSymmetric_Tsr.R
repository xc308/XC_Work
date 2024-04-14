#=============
# 14 Apr. 2024
#=============
# Aim:
  # a function that can check a 2D tensor matrix's symmetricity directly
    # on GPU device. 


isSymmetric_Tsr <- function(Tsr_mat) {
  
  Tsr_dim <- dim(Tsr_mat) # a vector
  
  if (length(Tsr_dim) != 2){
    stop("Input tensor must be a 2D matrix.")
  }
  
  Tsr_mat_t <- torch_t(Tsr_mat)
  
  eq_mat <- torch_eq(Tsr_mat, Tsr_mat_t)
  eq_prod <- torch_prod(eq_mat)
  
  
  if (as.array(eq_prod) == 1){
    print("sym: Yes")
  } else {
    print("sym: NO")
  }
  
 #torch_where(torch_eq(eq_prod, 1), print("sym: Yes"), print("sym: NO"))
 #torch_where(eq_prod > 0, print("sym: Yes"), print("sym: NO"))
}





