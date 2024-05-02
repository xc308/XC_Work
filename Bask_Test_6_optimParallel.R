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


cl <- makeCluster(2, outfile = "Bask_Test6_oP.log")
setDefaultCluster(cl = cl)

# Specify the library path on each worker
clusterEvalQ(cl, {
  .libPaths("/bask/projects/v/vjgo8416-xchen")
})

# Load the packages on each worker
clusterEvalQ(cl, {
  #library(torch)
  #library(GPUmatrix)
  library(optimParallel)
})


# Load the core R package (no need to specify the path)
clusterEvalQ(cl, {
  library(Matrix)
})

library(Matrix)
source("Fn_I_sparse.R")
I_sparse <- get("I_sparse")
clusterExport(cl, "I_sparse")


x <- rnorm(n = 500, mean = 5, sd = 2)
clusterExport(cl, "x")


negll <- function(par, x){
  scale_2 <- as.numeric(I_sparse(size = 1, value = 2))
  -scale_2 * sum(dnorm(x = x, mean = par[1], sd = par[2], log = TRUE))
} 



# 2 parameters, p = 2
# evaluation 2p + 1 = 5
# so want 5 tasks in parallel

#o2 <- optimParallel(par = c(1, 1), 
#                    fn = negll, 
#                    x = x, 
#                    method = "L-BFGS-B",
#                    lower = c(-Inf, 0.0001),
#                    control = list(trace = 0, 
#                                   maxit = 100,
#                                   pgtol = 1e-4))



# Perform parallel optimization
o2 <- tryCatch({
  optimParallel(par = c(1, 1), 
                fn = negll, 
                x = x, 
                method = "L-BFGS-B",
                lower = c(-Inf, 0.0001),
                control = list(trace = 0, 
                               maxit = 100,
                               pgtol = 1e-4))
}, error = function(e) {
  # Log any errors that occur during optimization
  cat("Error during optimization:", e$message, "\n", file = "Bask_Test6_oP.log", append = TRUE)
  NULL
})



stopCluster(cl)
o2

# $par
#[1] 5.099278 2.125822

#$value
#[1] 2173.097

#$counts
#function gradient 
#16       16 

#$convergence
#[1] 0

#$message
#[1] "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"

##----------
# conclusion
##----------

# when function contains other function, 
# still can converge

# although log file still catch the same error messages
# starting worker pid=3226805 on localhost:11490 at 21:53:54.340
#starting worker pid=3226804 on localhost:11490 at 21:53:54.340
#Error in unserialize(node$con) : error reading from connection
#Calls: <Anonymous> ... doTryCatch -> recvData -> recvData.SOCKnode -> unserialize
#Execution halted

# but these logged error does NOT affect the optimization
# and its convergence. 


