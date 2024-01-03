#=============
# 22 Nov. 2023
#=============
# Function:
  # This is a function creates Identity matrix using sparse matrix
  # Arg:
      # size: the dimension of the Identity matrix
      # value: the value on the main diagonal

I_sparse <- function(size, value) {
  
  require("Matrix")
  sparseMatrix(i = 1:size, j = 1:size, x = value)
  
}

