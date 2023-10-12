#==============
# 11 Oct. 2023
#==============

# Aim:
  # learn to create band sparse matrix in Matrix package
  # for the construction of bandsparse of SIGMA
  # ref:Matrix.pdf
  # 


bandSparse(n, m = n, k, diagonals, symmetric = FALSE,
    repr = "C", giveCsparse = (repr == "C"))

# Arguments:
  # n, m : the matrix dimension (n;m) = (nrow; ncol).
  # k: integers specifying the diagonals that are not set to zero. 
      # These are interpreted relative to the main diagonal, 
      # which is k=0. 
      # Positive and negative values of k indicate diagonals 
      # above and below the main diagonal, respectively.


  # diagonals: optional list of sub-/super- diagonals; 
      # if missing, the result will be a pattern
      # matrix, i.e., inheriting from class nMatrix

  # symmetric: logical; if true the result will be symmetric (inheriting from class symmetricMatrix)
      # and only the upper or lower triangle must be specified (via k and diagonals).

  # repr: character string, one of "C", "T", or "R", 
      # specifying the sparse representation
      # to be used for the result, i.e., one from the super classes CsparseMatrix, TsparseMatrix,
      # or RsparseMatrix.


# Values:
  # a sparse matrix (of class CsparseMatrix) of dimension n by m 
  # with diagonal “bands” as specified


diags <- list(1:30, 10*(1:20), 100*(1:20))
s1 <- bandSparse(13, k = -c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
s1

s2 <- bandSparse(13, k = c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
s2

# 13 x 13 sparse Matrix of class "dsCMatrix"

#[1,]   1  10 100   .   .   .  10   .   .    .    .    .    .
#[2,]  10   2  20 200   .   .   .  20   .    .    .    .    .
#[3,] 100  20   3  30 300   .   .   .  30    .    .    .    .
#[4,]   . 200  30   4  40 400   .   .   .   40    .    .    .
#[5,]   .   . 300  40   5  50 500   .   .    .   50    .    .
#[6,]   .   .   . 400  50   6  60 600   .    .    .   60    .
#[7,]  10   .   .   . 500  60   7  70 700    .    .    .   70
#[8,]   .  20   .   .   . 600  70   8  80  800    .    .    .
#[9,]   .   .  30   .   .   . 700  80   9   90  900    .    .
#[10,]   .   .   .  40   .   .   . 800  90   10  100 1000    .
#[11,]   .   .   .   .  50   .   .   . 900  100   11  110 1100
#[12,]   .   .   .   .   .  60   .   .   . 1000  110   12  120
#[13,]   .   .   .   .   .   .  70   .   .    . 1100  120   13

stopifnot(identical(s1, t(s2)), is(s1,"dsCMatrix"))


m5 <- Matrix(rnorm(25), ncol = 5)
band(m5, -1, 2) # general band

m4 <- matrix(rnorm(25), nrow = 5)
band(m4, -2, 2)

# 5 x 5 Matrix of class "dgeMatrix"
#[,1]       [,2]       [,3]         [,4]       [,5]
#[1,] -0.4538689 -0.2731564  0.7484807  0.000000000  0.0000000
#[2,] -0.4464357 -1.6621134 -0.6667239  0.233845293  0.0000000
#[3,] -0.7723462 -1.3448900 -0.2402593  0.005502958 -1.8770703
#[4,]  0.0000000  1.0750961 -0.1702701 -0.242091354 -0.1405182
#[5,]  0.0000000  0.0000000 -1.6317667  1.242963777 -0.9708745


#--------------------
# Try on SIGMA7_Chain
#--------------------
SIGMA_Chain_7_copy <- SIGMA_Chain_7

ncol(SIGMA_Chain_7_copy) # [1] 140

SIGMA_Chain_7_copy_Band <- band(SIGMA_Chain_7_copy, 
                                -0.6*ncol(SIGMA_Chain_7_copy), 
                                0.6*ncol(SIGMA_Chain_7_copy))

eigen(SIGMA_Chain_7_copy_Band, only.values = T)
all(eigen(SIGMA_Chain_7_copy_Band, only.values = T)$val > 0)
# [1] FALSE

Test_sym_pd(SIGMA_Chain_7_copy_Band)
# [1] "Symmetric: Yes"
[1] "p.d.: No"

plt_Sig(as.matrix(SIGMA_Chain_7_copy_Band), p = 7)



SIGMA_Chain_7_copy_Band2 <- band(SIGMA_Chain_7_copy, 
                                -0.9*ncol(SIGMA_Chain_7_copy), 
                                0.9*ncol(SIGMA_Chain_7_copy))

eigen(SIGMA_Chain_7_copy_Band2, only.values = T)
all(eigen(SIGMA_Chain_7_copy_Band2, only.values = T)$val > 0)


Test_sym_pd(SIGMA_Chain_7_copy_Band2)



#---------------
# Try band on SIGMA_Chain_7_large
#---------------

SIGMA_Chain_7_large_CP <- SIGMA_Chain_7_large
Test_sym_pd(SIGMA_Chain_7_large_CP)

SIGMA_Chain_7_large_CP_Bd <- band(SIGMA_Chain_7_large_CP, -0.6*ncol(SIGMA_Chain_7_large_CP), 
     0.6*nrow(SIGMA_Chain_7_large_CP))

Test_sym_pd(SIGMA_Chain_7_large_CP_Bd)
# [1] "Symmetric: Yes"
#[1] "p.d.: Yes"


for(bdw in seq(0.1, 0.9, by = 0.1)) {
  SIGMA_Chain_7_large_CP_Bd <- band(SIGMA_Chain_7_large_CP, -bdw*ncol(SIGMA_Chain_7_large_CP), 
                                    bdw*nrow(SIGMA_Chain_7_large_CP))
  
  print(bdw)
  Test_sym_pd(SIGMA_Chain_7_large_CP_Bd)
  
}

# [1] 0.1
# [1] "Symmetric: Yes"
# [1] "p.d.: No"
# [1] 0.2
# [1] "Symmetric: Yes"
# [1] "p.d.: No"
# [1] 0.3
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.4
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.5
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.6
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.7
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.8
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"
# [1] 0.9
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


plt_Sig(SIGMA_Chain_7_large, p = 7)

diag(SIGMA_Chain_7_large)
diag(SIGMA_Chain_7)




#----------------------------------------
# Compare the speed of cholesky inverting
#----------------------------------------

Fn_Chol_Inv <- function(SIGMA) {
  
  chol2inv(chol(SIGMA))
}

library(microbenchmark)

microbenchmark(
  
  Fn_Chol_Inv(SIGMA_Chain_7_copy_Band2),
  Fn_Chol_Inv(SIGMA_Chain_7_copy)
)


# Unit: microseconds
#                              expr     min       lq     mean   median      uq
#Fn_Chol_Inv(SIGMA_Chain_7_copy_Band2) 454.334 471.2300 486.2277 480.4585 500.459
#      Fn_Chol_Inv(SIGMA_Chain_7_copy) 694.751 726.5005 763.5182 754.5425 788.605
#     max neval
#544.959   100
#980.375   100

# 

486.2277 / 763.5182 # [1] 0.6368253

#-------------
# Conclusions
#-------------
# with bandwidth 0.9, i.e., only the last 10% of off-diag elements
  # away from the main diagonal elements are set to exact zeroes, 

# the quantile for computational time for cholesky inversion upon 
# 100 runs of evaluations is 

# Unit: microseconds
#                              expr     min       lq     mean   median      uq
#Fn_Chol_Inv(SIGMA_Chain_7_copy_Band2) 454.334 471.2300 486.2277 480.4585 500.459
#      Fn_Chol_Inv(SIGMA_Chain_7_copy) 694.751 726.5005 763.5182 754.5425 788.605

#     max neval
#544.959   100
#980.375   100





