#==============
# 23 Nov. 2023
#==============

# Aim:
  # to visualize wave_5 function in 1D and 3D to better understand its property


#---------
# Settings
#---------
library(plot3D)
getwd()
img_path <- "./Results/"


#--------
# Domain
#--------
ds <- 0.01
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

# for 3D requires displacement
H <- t(outer(s, s, "-"))


#------
# Range check
#------
source("Fn_Waves.R")

DELTA <- c(-0.7, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 1)
R <- matrix(NA, ncol = 3)
for (delta in DELTA){
  w5 <- wave_v5(h = H, delta = delta, A = 1)
  
  r <- range(w5)
  dr <- t(c(delta, r))
  
  R <- rbind(R, dr)
}
R
colnames(R) <- c("delta", "min", "max")
R[-1, ]
#     delta        min max
# [1,]  -0.7 -0.4693878   1
#[2,]  -0.3 -1.0000000   1
#[3,]  -0.1 -1.0000000   1
#[4,]   0.1 -1.0000000   1
#[5,]   0.3 -1.0000000   1
#[6,]   0.5 -1.0000000   1
#[7,]   0.7 -0.4693878   1
#[8,]   1.0 -0.6200000   1


#-----
# 3D
#----
## positive dlt
jpeg(paste0(img_path, "wave_v5_3D_pos_dlt.jpeg"),
    width = 8, height = 7, units = "in", res = 300)

par(mfrow = c(3, 2), mar = c(1.5, 1, 3.5, 2.5))
DELTA <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
for (delta in DELTA) {
  w5 <- wave_v5(h = H, delta = delta, A = 1)
  w5_mat <- matrix(w5, nrow = nrow(H))
  
  main_title <- as.expression(paste("Tri_Wave (A = 1) \n delta = ", delta))
  persp3D(s, s, w5_mat, theta = 135, 
          main = main_title)
}

dev.off()


## Neg dlt
jpeg(paste0(img_path, "wave_v5_3D_neg_dlt.jpeg"),
     width = 8, height = 7, units = "in", res = 300)

par(mfrow = c(3, 2), mar = c(1.5, 1, 3.5, 2.5))
DELTA <- -1 * c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
for (delta in DELTA) {
  w5 <- wave_v5(h = H, delta = delta, A = 1)
  w5_mat <- matrix(w5, nrow = nrow(H))
  
  main_title <- as.expression(paste("Tri_Wave (A = 1) \n delta = ", delta))
  persp3D(s, s, w5_mat, theta = 135, 
          main = main_title)
}
dev.off()


#---------
# 1D plots
#---------
length(H) # 400
length(s) # 20
plt_w56 <- function(dlt){
  
  # plot wave_v5 at each grid of s
  w_5 <- wave_v5(h = s, delta = dlt, A = 1)
  #w_6 <- wave_v6(h = s, delta = dlt, A = 1)
  plot(s, w_5, type = "l", 
       main = bquote(atop("Tri_wave_v5 (A = 1)", 
                          "delta = "~ .(dlt))))
  
}

jpeg(paste0(img_path, "Wave_v5_1D_pos_dlt.jpeg"),
     width = 8, height = 7, units = "in", res = 300)

jpeg(paste0(img_path, "Wave_v6_1D_pos_dlt.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(2.5, 2.5, 3.5, 2.5))
DELTA <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)

for (dlt in DELTA){
  plt_w56(dlt)
} 
dev.off()



## neg-delta
jpeg(paste0(img_path, "Wave_v5_1D_neg_dlt.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
par(mfrow = c(3, 2), mar = c(2.5, 2.5, 3.5, 2.5))
DELTA <- -1*c(0.1, 0.3, 0.5, 0.7, 0.9, 1)

for (dlt in DELTA){
  plt_w56(dlt)
} 
dev.off()


#=======================
# Settings 2: large grid
#=======================

ds <- 10
s <- seq(-100 + ds/2, 100 - ds/2, by = ds)

jpeg(paste0(img_path, "Wave_v5_1D_pos_dlt_big_grid.jpeg"),
     width = 8, height = 7, units = "in", res = 300)

DELTA <- c(10, 30, 50, 70, 90, 100)
par(mfrow = c(3, 2), mar = c(2.5, 2.5, 3.5, 2.5))

for (dlt in DELTA){
  plt_w5(dlt)
} 
dev.off()



#------
# large grid but very small dlt, what happens
#------
ds <- 10
s <- seq(-100 + ds/2, 100 - ds/2, by = ds)

#jpeg(paste0(img_path, "Wave_v5_1D_pos_dlt_big_grid.jpeg"),
     #width = 8, height = 7, units = "in", res = 300)

DELTA <- c(0.1, 0.5, 1, 3, 5)
par(mfrow = c(3, 2), mar = c(2.5, 2.5, 3.5, 2.5))

for (dlt in DELTA){
  plt_w5(dlt)
} 



wave_v5(h = s, delta = 0.1, A = 1)


#==============
# side by side
#==============
install.packages("gridExtra")
library(gridExtra)

install.packages("jpeg")
library(jpeg) # for readJPEG

install.packages("grid")
library(grid) # for rasterGrob


getwd()
image.path
# [1] "./Results/"


# Read the JPEG files
image1 <- readJPEG(paste0(image.path, "wave_v5_3D_pos_dlt.jpeg"))
image2 <- readJPEG(paste0(image.path, "wave_v5_3D_neg_dlt.jpeg"))


# Create plots from the images
plot1 <- rasterGrob(image1, interpolate = TRUE)
plot2 <- rasterGrob(image2, interpolate = TRUE)


# Arrange the plots side by side
SBS_plots <- grid.arrange(plot1, plot2, ncol = 2)


# Save the arranged plots as JPEG
jpeg(paste0(image.path, "wave_v5_SBS.jpeg"),
     width = 8, height = 7, units = "in", res = 300)
grid.draw(SBS_plots)
dev.off()
