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


x <- rnorm(n = 500, mean = 5, sd = 2)
negll <- function(par, x) -sum(dnorm(x = x, mean = par[1], sd = par[2], log = TRUE))

# 2 parameters, p = 2
# evaluation 2p + 1 = 5
# so want 5 tasks in parallel

cl <- makeCluster(5); setDefaultCluster(cl = cl)

o2 <- optimParallel(par = c(1, 1), 
                    fn = negll, 
                    x = x, 
                    lower = c(-Inf, 0.0001))

o2



