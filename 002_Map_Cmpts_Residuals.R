#==========================#
# Map the Cmpts' Residuals
#==========================#

# if residuals exhibit spatial patterns, 
# then justify further spatial modelling



install.packages("sp")
install.packages("raster")

install.packages("lattice", .Library)
install.packages("latticeExtra")
install.packages("rasterVis") # ref:https://gis.stackexchange.com/questions/342488/error-in-implementing-levelplot-in-r-for-a-categorical-raster-with-rat

install.packages("RColorBrewer")


library(raster)
library(sp)

library(lattice)
library(latticeExtra)
library(rasterVis)

library(RColorBrewer)


#---------------------#
# Rasterize data frame
#---------------------#

# ref:https://stackoverflow.com/questions/19627344/how-to-create-a-raster-from-a-data-frame-in-r

df_Cmpts_Res_wide_log_16 <- df_Cmpts_Res_wide_log %>%
  filter(Year == 2016)

str(df_Cmpts_Res_wide_log_16)  
# 'data.frame':	27384 obs. of  9 variables


Rst_Cmpt_Res_log_16 <- list()
df_Cmpt_Res_log_16 <- list()
NM <- paste(c("BC", "DU", "OM", "SS", "SU", "PM25"), "_Residuals_log", sep = "")
for (i in seq_along(NM)) {
  df_Cmpt_Res_log_16[[i]] <- df_Cmpts_Res_wide_log_16 %>%
    select(Lon, Lat, NM[i])
  
  Rst_Cmpt_Res_log_16[[i]] <- rasterFromXYZ(df_Cmpt_Res_log_16[[i]])
  
}


Bk_Res_log_16 <- brick(Rst_Cmpt_Res_log_16[[1]], Rst_Cmpt_Res_log_16[[2]],
      Rst_Cmpt_Res_log_16[[3]], Rst_Cmpt_Res_log_16[[5]],
      Rst_Cmpt_Res_log_16[[4]], Rst_Cmpt_Res_log_16[[6]])


#class      : RasterBrick 
#dimensions : 186, 480, 89280, 6  (nrow, ncol, ncell, nlayers)
#resolution : 0.75, 0.75  (x, y)
#extent     : -179.625, 180.375, -55.875, 83.625  (xmin, xmax, ymin, ymax)
#crs        : NA 
#source     : memory
#names      : BC_Residuals, DU_Residuals, OM_Residuals, SS_Residuals, SU_Residuals, PM25_Residuals 
#min values :        -3.510039,        -5.210456,        -2.779548,        -3.233773,        -2.921836,          -2.457869 
#max values :         4.026947,         7.231854,         3.832994,         3.050802,         4.886131,           3.515366 

names(Bk_Res_log_16) <- c("BC_Residuals", "DU_Residuals",
                          "OM_Residuals", "SU_Residuals",
                          "SS_Residuals", "PM25_Residuals")


Bk_Res_log_16


#--------#
# Ticks
#--------#

tick_at_log <- c(-5, -3, -2, 0, 1, 3, 4, 7)
exp(tick_at_log)
# [1] 6.737947e-03 4.978707e-02
# [3] 1.353353e-01 1.000000e+00
# [5] 2.718282e+00 2.008554e+01
# [7] 5.459815e+01 1.096633e+03

labels_expback_total <- c(0.005, 0.05, 0.1, 1, 3, 5, 50, 1000)


#-------------#
# chose color
#-------------#

cols <- rev(brewer.pal(11, name = "RdBu"))
brewer.div <- colorRampPalette(cols, interpolate = "spline")


#------#
# plot
#------#

plt <- levelplot(Bk_Res_log_16, cuts = 499,
                 col.regions = brewer.div(500),
                 layout = c(2, 3),
                 colorkey = list(labels = list(labels = labels_expback_total),
                                 width = 0.7))

plt + latticeExtra::layer(sp.polygons(WHO_map, col = "black", lwd = 0.05))











#----------#
# Quantile
#----------#

quantile(df_Cmpts_Res_wide_log$BC_Residuals_log)
#           0%          25%          50%          75%         100% 
# -3.641150795 -0.656487485 -0.006686259  0.599218713  5.976576612 

quantile(df_Cmpts_Res_wide_log$DU_Residuals_log)
#         0%        25%        50%        75%       100% 
# -5.6810601 -0.8726389 -0.0434291  0.7836720  7.3503436 

quantile(df_Cmpts_Res_wide_log$OM_Residuals_log)
#          0%         25%         50%         75%        100% 
# -3.06485639 -0.58965103  0.01771825  0.55726757  4.53982242 

quantile(df_Cmpts_Res_wide_log$SS_Residuals_log)
#          0%         25%         50%         75%        100% 
# -3.36549056 -0.61581700 -0.03256249  0.61464714  3.25795270 

quantile(df_Cmpts_Res_wide_log$SU_Residuals_log)
#          0%         25%         50%         75%        100% 
# -3.16071072 -0.53046114 -0.05759097  0.51741907  4.90779859 

quantile(df_Cmpts_Res_wide_log$PM25_Residuals_log)
#          0%         25%         50%         75%        100% 
#      -2.4578694 -0.4279401 -0.0403321  0.4013001  4.3661661 

# -5, -0.5, -0.05, 0.4, 0.5, 0.5, 5, 7






