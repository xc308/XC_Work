##******************************************#
## Test matrix symetric and Postive definite
##******************************************#

Test_sym_pd <- function(test.mat) {
  if (isSymmetric(test.mat) == TRUE) {
    print("Symmetric: Yes")
  } else print("Symmetric: No")
  
  if (all(eigen(test.mat)$values > 0)) {
    print("p.d.: Yes")
  } else print("p.d.: No")
}






