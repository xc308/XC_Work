#=============
# 23 Apr. 2024
#=============

# Aim:
  # First test whether enviornment for parallle makeCluster and others is there. 



library(parallel)    
cores <- detectCores()
cr <- cores - 1
cl <- makeCluster(cr)

setDefaultCluster(cl = cl) 

parSapply(cl, 1:5, sqrt)
stopCluster(cl)