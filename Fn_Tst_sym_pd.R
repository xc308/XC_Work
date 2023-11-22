#=============
# 22 Nov. 2023
#=============

# Function:
  # This is a function to test the symmetry and p.d. of a matrix
  # If a matrix is symmetry, 
    # keep testing if all eigen values are positive; 

  # If a matrix is not symmetry, 
    # it must not be p.d.


Tst_sym_pd <- function(tst.M){
  
  if (isSymmetric(tst.M) == TRUE){
    print("Symmetric: Yes")
    
    if (all(eigen(tst.M, symmetric = T, only.values = T)$val > 0)){
      print("p.d.: Yes")
    } else {print("p.d.: No")}
    
  }else{print("Symmetric: No")}
}


#======
# Test
#======

#Tst_sym_pd(SIGMA_inv_5_pert)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

#Tst_sym_pd(SG_SG_inv_5_TST2$SIGMA_inv)
# [1] "Symmetric: Yes"
# [1] "p.d.: No"

