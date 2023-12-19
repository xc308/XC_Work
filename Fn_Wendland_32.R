#============
# 3 Dec. 2023
#============

# Function:
  # Wendland function 
  
# Arg:
  # r: displacement of pairs of locations


# Parameters:
  # mannually set:
    # 1. k: order determine the smoothness (continuously differentiable)
          # the higher the smoother
          # k = 3/2 to match matern nu = 3/2

    # 2. R: compact support radius, beyond which the function
          # value is exact 0. 


  # leave for inference:
    # A: amplitude
    # dlt: translation horizontally


WendLd_32 <- function(r, R, dlt, A){
  
  ifelse(abs(r - dlt) <= R, A*(1 - abs(r - dlt)/R)^4 * (1 + 4 * abs(r - dlt)/R), 0)
  
}


#======
# Plot
#======

#r <- seq(-1, 1, by = 0.01)
#R <- 0.5

#A = 1
#dlt = 0.4

#par(mfrow = c(1,1))
#Wdld_val <- WendLd_32(r = r, R = R, dlt = 0.4, A = 1)
#plot(r, Wdld_val, type = "l",
     #xlab = "r", ylab = "Wend vals", 
     #main = "Wendland Function (k = 3/2, dlt = 0.4)")

