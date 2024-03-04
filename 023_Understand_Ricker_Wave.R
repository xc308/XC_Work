#=============
# 5 Dec. 2023
#=============

# Aim:
  # Understand Ricker Wavelet


# Define the Ricker wavelet function
ricker_wavelet <- function(t, f0) {
  term1 <- (1 - 2 * pi^2 * f0^2 * t^2)
  term2 <- exp(-pi^2 * f0^2 * t^2)
  result <- term1 * term2
  return(result)
}

# Time values
t_values <- seq(-0.2, 0.2, length.out = 1000)  # Adjust time range and resolution as needed
t_values <- seq(-1, 1, by = 0.01)
t_values <- seq(-1, 1, by = 0.1)
t_values <- seq(-10, 10, by = 0.1)
t_values <- seq(-100, 100, by = 1)
t_values <- seq(-100, 100, by = 10)
t_values <- seq(-10000, 10000, by = 80)

# Peak frequency
f0 <- 100
f0 <- 50
f0 <- 25  # Adjust the peak frequency as desired
f0 <- 5
f0 <- 10
f0 <- 1
f0 <- 0.5
f0 <- 0.2
f0 <- 0.1
f0 <- 0.05
f0 <- 0.01
f0 <- 0.001
f0 <- 0.0001


# Calculate the Ricker wavelet values
wavelet_values <- ricker_wavelet(t_values, f0)

# Plot the Ricker wavelet
plot(t_values, wavelet_values, type = "l", col = "blue",
     xlab = "Time", ylab = "Amplitude", 
     main = paste("Ricker Wavelet with f0 =", f0))



#-----------------
# Try combinations
#-----------------
t_values <- seq(-1, 1, by = 0.01)
f0 <- 1


t_values <- seq(-100, 100, by = 1)
f0 <- 0.01


t_values <- seq(-1, 1, by = 0.1)
f0 <- 1


#---------------
# With translate
#---------------
ricker_wavelet_dlt <- function(t, f0, A, dlt) {
  term1 <- (1 - 2 * pi^2 * f0^2 * (abs(t - dlt))^2)
  term2 <- exp(-pi^2 * f0^2 * (abs(t - dlt))^2)
  result <- A * term1 * term2
  return(result)
}


t_values <- seq(-1, 1, by = 0.1)
t_values <- seq(-1, 1, by = 0.01)
f0 <- 1

dlt <- 0.5
dlt <- 0
R_w_dlt <- ricker_wavelet_dlt(t = t_values, f0 = f0,  A = 1, dlt = dlt)
R_w_dlt <- ricker_wavelet_dlt(t = t_values, f0 = f0,  A = 1, dlt = dlt)


jpeg(paste0(image.path, "Ricker_wave_dlt05.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
plot(t_values, R_w_dlt, type = "l", col = "blue",
     xlab = "t", ylab = "Ricker Wave Value", 
     main = paste("Ricker Wavelet", "\n", "f0 =", f0,
                  ", dlt =", dlt))
dev.off()


jpeg(paste0(image.path, "Ricker_wave_dlt0.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
plot(t_values, R_w_dlt, type = "l", col = "blue",
     xlab = "t", ylab = "Ricker Wave Value", 
     main = paste("Ricker Wavelet", "\n", "f0 =", f0,
                  ", dlt =", dlt))
dev.off()


#-------------
# side by side
#-------------

jpeg(paste0(image.path, "Ricker_wave_SBS.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
par(mfrow = c(1, 2))

dlt <- 0
R_w_dlt <- ricker_wavelet_dlt(t = t_values, f0 = f0,  A = 1, dlt = dlt)
plot(t_values, R_w_dlt, type = "l", 
     xlab = "t", ylab = "Ricker Wave Value", 
     main = paste("Ricker Wavelet", "\n", "f0 =", f0,
                  ", dlt =", dlt))


dlt <- 0.5
R_w_dlt <- ricker_wavelet_dlt(t = t_values, f0 = f0,  A = 1, dlt = dlt)
plot(t_values, R_w_dlt, type = "l", 
     xlab = "t", ylab = "Ricker Wave Value", 
     main = paste("Ricker Wavelet", "\n", "f0 =", f0,
                  ", dlt =", dlt))
dev.off()




#--------------------
# decay f0 not square
#--------------------
ricker_wavelet_dlt_2 <- function(t, f0, A, dlt) {
  term1 <- (1 - 2 * pi^2 * f0 * (abs(t - dlt))^2)
  term2 <- exp(-pi^2 * f0 * (abs(t - dlt))^2)
  result <- A * term1 * term2
  return(result)
}

t_values <- seq(-1, 1, by = 0.01)
t_values <- seq(-100, 100, by = 1)
f0 <- 0.01

R_w_value <- ricker_wavelet_dlt_2(t = t_values, f0 = f0, A = 1, dlt = 0)

plot(t_values, R_w_value, type = "l")


#----------
# Discover
#----------
# if turn f0^2 to f0, then when grid becomes large
  # and the decay f0 needs to lower to echo the low
  # oscillation, f0^2 will speed up the decrease of oscialtion
  # the curve will be broader 







