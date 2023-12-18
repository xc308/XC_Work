#=============
# 18 Dec.2023
#=============

# Function:
  # Aim: turn a non-neg matrix into row-stochastic
        # by dividing each element of the matrix with
        # its row sum

  # Arg:
        # A_mat


row_stoc <- function(A_mat) {
  scaling_factor <- 1 / rowSums(A_mat)
  return(scaling_factor * A)
}


#=====#
# Test
#=====#

#A = matrix(c(1, 2, 3, 4), nrow = 2)
#row_stoc(A)
#eigen(row_stoc(A), only.values = T)
# $values
# [1]  1.00000000 -0.08333333

