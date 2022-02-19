## =========================================  
## Wishart and Inverse Wishart distribution 
## =========================================  

# ref: https://github.com/johnnyzhz/wishartprior/blob/main/R/wishart.R

stat.iWishart <- function(sigma, m){
  p <- nrow(sigma)
  m.iwishart <- sigma/(m-p-1)
  diag.sigma <- diag(sigma)
  v.iwishart <- (m-p+1)*sigma^2 + (m-p-1)*diag.sigma%*%t(diag.sigma)
  v.iwishart <- v.iwishart/((m-p)*(m-p-1)^2*(m-p-3))
  list(mean=m.iwishart, variance=v.iwishart)
}


#----------------#
# stats.iWishart
#----------------#

# this function calculates the mean and variance of inverse.Wishart dist
## sigma: the scale matrix
## df: degrees of freedom

stats.iWishart <- function(sigma, df) {
  p <- nrow(sigma)
  m.iWishart <- sigma / (df - p - 1)
  
  diag.sigma <- diag(sigma)
  v.iWishart <- (df - p + 1) * sigma^2 + (df - p - 1) * diag.sigma %*% t(diag.sigma)
  v.iWishart <- v.iWishart / ((df - p) * (df - p - 1)^2 * (df - p - 3))
  
  list(Mean = m.iWishart, Variance = v.iWishart)
}


# try
## A = identity 2 by 2
A = diag(1, 2)
stats.iWishart(sigma = A, df = 1e6)



#----------------#
# stats.Wishart
#----------------#

# this function calculates the mean and var of Wishart distribution
# sigma : scale matrix
# df : degrees of freedom
stats.Wishart <- function(sigma, df) {
  m.Wishart <- df * sigma
  
  diag.sigma <- diag(sigma)
  var.Wishart <- df * (sigma^2 + diag.sigma %*% t(diag.sigma))
  
  list(Mean = m.Wishart, Variance = var.Wishart)
}

# try
stats.Wishart(sigma = A, df = 1e6)



#--------------------#
# posterior iWishart
#--------------------#
# given the MV data with var-cov matrix Sigma
# and a iWishart prior for this Sigma, i.e. iWishart(V0, m0)
# the posterior [Sigam | D] ~ iWishart(nS + V0, n + m0) 
# n: sample data size
# S: sample data var-covariance matrix
# m0: df for prior iWishart distribution
# V0: var-cov for prior iWishart distribution

posterior.iWishart <-function(n, S, m0, V0) {
  m1 <- n + m0
  V1 <- n * S + V0
  
  stats.iWishart(V1, m1)
}



#---------------------#
# Posterior of Wishart
#---------------------#

# when MVN data use precision matrix P in its density,
# the prior set for this precision matrix P is Wishart (U0.inv, w0)
# then the posterior is still Wishart 
# [P | D] ~ Wishart (U1.inv, w1)


posterior.Wishart <- function(n, S, U0, w0) {
  w1 <- n + w0
  U1.inv <- solve(n * S + solve(U0))
  
  stats.Wishart(U1.inv, w1)
}




