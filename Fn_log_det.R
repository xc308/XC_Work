#=============
# 19 Jan.2024
#=============

# Aim:
  # Function to calculate log det (SIGMA) part in log likelihood
  # Arg: Chol_U a cholesky factor of SIGMA

log_det <- function(Chol_U) {
  Uii <- diag(Chol_U)
  return(2 * sum(log(Uii)))
}






