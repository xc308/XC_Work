#==========================#
# Map the Cmpts' Residuals
#==========================#

# if residuals exhibit spatial patterns, 
# then justify further spatial modelling


install.packages(c("latticeExtra", "lattice","rasterVis"))
install.packages("RColorBrewer")
install.packages("raster")


library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)
library(raster)



#---------------------#
# Rasterize data frame
#---------------------#

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
      Rst_Cmpt_Res_log_16[[3]], Rst_Cmpt_Res_log_16[[4]],
      Rst_Cmpt_Res_log_16[[5]], Rst_Cmpt_Res_log_16[[6]])






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


#--------#
# Ticks
#--------#

tick_at_log <- c(-5, -0.5, -0.05, 0.4, 0.5, 0.6, 5, 7)
exp(tick_at_log)
# [1] 6.737947e-03 6.065307e-01
# [3] 9.512294e-01 1.491825e+00
# [5] 1.648721e+00 1.822119e+00
# [7] 1.484132e+02 1.096633e+03

labels_expback_total <- c(0.005, 0.5, 1.0, 1.5, 2, 10, 100, 1000)


#-------#
# color
#-------#

cols_ryb <- rev(brewer.pal(11, name = "RdYlBu"))
brewer.div.ryb <- colorRampPalette(cols_ryb, interpolate = "spline")



plt <- levelplot(B_log, cuts = 499,
                 col.regions = brewer.div(500),
                 layout = c(2, 3),
                 colorkey = list(labels = list(labels = labels_expback_total),
                                 width = 0.7))

plt + latticeExtra::layer(sp.polygons(WHO_map, col = "black", lwd = 0.05))






