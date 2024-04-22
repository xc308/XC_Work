#=============
# 22 Apr. 2024
#=============

# Aim:
  # try optimParallel() 

# Ref:
  # optimParallel 


#install.packages("optimParallel")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library("optimParallel")


x <- rnorm(n = 1e7, mean = 5, sd = 2)
negll <- function(par, x) -sum(dnorm(x = x, mean = par[1], sd = par[2], log = TRUE))

cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)

o2 <- optimParallel(par = c(1, 1), 
                    fn = negll, 
                    x = x, 
                    lower = c(-Inf, 0.0001))

o2

