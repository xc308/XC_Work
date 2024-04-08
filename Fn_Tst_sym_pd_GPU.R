#============
# 8 Apr. 2024
#============
# Aim:
  # This function modifies "Fn_Tst_sym_pd.R" to be able to Test GPU matrix obj


# Arg:
  # tst.M.GPU


Tst_sym_pd_gpu <- function(tst.M.GPU){
  
  if (isSymmetric(as.matrix(tst.M.GPU)) == TRUE){
    print("Symmetric: Yes")
    
    eig_gpuM <- eigen(tst.M.GPU)$val # eigen value, real & img
    
    if (all(as.array(eig_gpuM@gm$real > 0))){
      print("p.d.: Yes")
    } else {print("p.d.: No")}
    
  }else{print("Symmetric: No")}
}


