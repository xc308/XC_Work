#**********************#
# Conditional Modelling
#**********************#

head(df_Cmpts_Res_wide_log)
str(df_Cmpts_Res_wide_log)
# 'data.frame':	136920 obs. of  9 variables:
  #$ Lon               : num  -42.8 -42 -39.8 -38.2 -37.5 ...
  #$ Lat               : num  83.2 83.2 83.2 83.2 83.2 ...
  #$ Year              : num  2016 2016 2016 2016 2016 ...
  #$ BC_Residuals_log  : num  -1.37 -1.51 -1.54 -1.54 -1.51 ...
  #$ DU_Residuals_log  : num  0.27611 0.11173 0.16029 -0.00793 -0.02871 ...
  #$ OM_Residuals_log  : num  -1.31 -1.45 -1.49 -1.5 -1.48 ...
  #$ SS_Residuals_log  : num  -1.08 -1.2 -1.21 -1.39 -1.46 ...
  #$ SU_Residuals_log  : num  -0.672 -0.795 -0.806 -0.824 -0.804 ...
  #$ PM25_Residuals_log: num  -0.77 -0.903 -0.938 -1.003 -1.009 ...


#=========#
# Packages
#=========#

## core operation
library(dplyr)
library(tidyr)
library(Matrix)

## plot
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

## sp
library(sp)
library(maptools)
library(mapproj)
library(RandomFields)

## computing
library(foreach)
library(doParallel)

# triangualtion and tesselation
library(INLA)
install.packages("deldir")
library(deldir) 















