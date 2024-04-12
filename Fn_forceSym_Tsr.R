#==============
# 12 Apr. 2024
#==============

# Aim:
  # A function that will force a 2D tensor matrix to be symmetric

# Arg:
  # a 2D tensor matrix on GPU

# Return:
  # a tensor forced symmetric


forceSym_Tsr <- function(Tensor_mat) {
  
  sym_tsr_mat <- (torch_t(Tensor_mat) + Tensor_mat) / 2
  
  return(sym_tsr_mat)
}







