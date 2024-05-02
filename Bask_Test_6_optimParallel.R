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

source("Fn_I_sparse.R")
I_sparse <- get("I_sparse")



cl <- makeCluster(2, outfile = "Bask_Test6_oP.log")
setDefaultCluster(cl = cl)

clusterExport(cl, "I_sparse")


x <- rnorm(n = 500, mean = 5, sd = 2)
clusterExport(cl, "x")


negll <- function(par, x){
  scale_2 <- I_sparse(size = 2, value = 1)
  -scale_2 * sum(dnorm(x = x, mean = par[1], sd = par[2], log = TRUE))
} 



# 2 parameters, p = 2
# evaluation 2p + 1 = 5
# so want 5 tasks in parallel

o2 <- optimParallel(par = c(1, 1), 
                    fn = negll, 
                    x = x, 
                    method = "L-BFGS-B",
                    lower = c(-Inf, 0.0001),
                    control = list(trace = 0, 
                                   maxit = 100,
                                   pgtol = 1e-4))

stopCluster(cl)
o2



