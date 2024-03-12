#==========================
# 3 Dec. 2023; 12 Mar. 2024
#===========================



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


#===============
# 2D Wendland_32
#===============

WendLd32_2D <- function(shft_dst_mat, A){
  
  R <- round(quantile(shft_dst_mat)[3], 1) # find the 50% quantile of shft_dst
  ifelse(shft_dst_mat <= R, A*(1 - shft_dst_mat/R)^4 * (1 + 4 * shft_dst_mat/R), 0)
  
}


#WL_2D_dlt_0204 <- WendLd32_2D(shft_dst_mat = Shft_dst_mat_0204, A = 1)
#quantile(WL_2D_dlt_0204)

#    0%          25%          50%          75%         100% 
#0.000000e+00 0.000000e+00 1.838579e-05 1.913940e-01 8.219221e-01





#===============
# 1D Wendland_32
#===============

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

