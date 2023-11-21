#=============
# 21 Nov. 2023
#=============
# Funtion:
    # Aim: for plotting the image of a constructed SIGMA or SIGMA_inv
    # Args: 
        # Sigma: a matrix form of SIGMA or SIGMA_inv
        # p: the number of variate fields

plt_Sig <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = paste("p = ", p))
}