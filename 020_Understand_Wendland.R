#========================
# 25, 30 Nov. 3 Dec. 2023
#========================

# Aim:
  # Understand Wendland function

# Arg:
  # r: displacement
  # R: support radius beyond which the Wendland value is exact 0


#--------
# setting
#--------

image.path <- "./Results/"



#--------------
# Basic version
#--------------

WendLd_32 <- function(r, R) {
  
  ifelse( abs(r) < R, (1 - abs(r)/R)^4 * (1 + 4*abs(r)/R), 0)
  
}


r <- seq(-1, 1, by = 0.01)
str(r)
R <- 0.5

WendL_vals <- WendLd_32(r = r, R = R)
str(WendL_vals) # num [1:21]



plot(r, WendL_vals, type = "l",
     xlab = "r", ylab = "Wendland",
     main = "Wendland Function (k = 3/2)")


#--------------------
# with skew parameter
#--------------------

WendLd_32_dlt <- function(r, R, dlt){
  
  ifelse(abs(r - dlt) < R, (1 - abs(r - dlt)/R)^4 * (1 + 4 * abs(r - dlt)/R), 0)
  
}


#------
# plot
#------

r <- seq(-1, 1, by = 0.01)
R <- 0.5

jpeg(paste0(image.path, "Wendland.jpeg"), 
     width = 8, height = 7, units = "in", res = 300)
par(mfrow = c(1, 2))

WendL_vals <- WendLd_32(r = r, R = R)
plot(r, WendL_vals, type = "l",
     xlab = "r", ylab = "Wendland vals",
     main = "Wendland Function (k = 3/2, dlt = 0)")

WendL_vals_dlt <- WendLd_32_dlt(r = r, R = R, dlt = 0.4)
plot(r, WendL_vals_dlt, type = "l", 
     xlab = "r", ylab = "Wendland vals",
     main = "Wendland (k = 3/2, dlt = 0.4)")

dev.off()


