#**********************#
# Conditional Modelling
#**********************#

## contiue from 001
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



#=======#
# Data
#=======#

## Work on 2016 Cmpts Res log data

head(df_Cmpts_Res_wide_log)

df_Cmpt_Res_log_16 <- df_Cmpts_Res_wide_log %>%
  filter(Year == 2016)

str(df_Cmpt_Res_log_16)

##'data.frame':	27384 obs. of  9 variables:
    #$ Lon               : num  -42.8 -42 -39.8 -38.2 -37.5 ...
    #$ Lat               : num  83.2 83.2 83.2 83.2 83.2 ...
    #$ Year              : num  2016 2016 2016 2016 2016 ...
    #$ BC_Residuals_log  : num  -1.37 -1.51 -1.54 -1.54 -1.51 ...
    #$ DU_Residuals_log  : num  0.27611 0.11173 0.16029 -0.00793 -0.02871 ...
    #$ OM_Residuals_log  : num  -1.31 -1.45 -1.49 -1.5 -1.48 ...
    #$ SS_Residuals_log  : num  -1.08 -1.2 -1.21 -1.39 -1.46 ...
    #$ SU_Residuals_log  : num  -0.672 -0.795 -0.806 -0.824 -0.804 ...
    #$ PM25_Residuals_log: num  -0.77 -0.903 -0.938 -1.003 -1.009 ...

## the unit for log_scale residuals are all mug/cm^3



#================#
# Form Z function
#================#

## This function concatenate two variables of a df into one vector
  # and is able to create a reverse version 

source("F003_01_Form_Z.R")

Z_BC_OM <- Form_Z(df = df_Cmpt_Res_log_16, col_num1 = 4, col_num2 = 6, model_num = 1)

str(Z_BC_OM)  # num [1:54768, 1]



#=========#
# Obs size
#=========#

m1 <- m2 <- nrow(df_Cmpt_Res_log_16)
m  <-  m1 + m2  # total obs



#=====================#
# Process Discretization
#=====================#

## Each Cmpt has its own true process
## Follow the Finite Element Method's idea, 
  ## we dicretize each of the process into a collection of simple shape elements (or subdomains), altogether termed as mesh
  ## and the process over each subdomain can itself be approximated by a 2D planar surface
  
## e.g. using INLA to create a mesh of triangles with 3 nodes each
  ## we then use deldir() fn to tesselate base on each of tri node, and 
  ## by tri and v tess propoerty, the tri nodes are exactly the centriods of v tess areas
  ## which is easier to identify the simple element using just 1 point each and
  ## from the output summary, the area of each tesslation is ready-made. 


#-----------------#
# Mesh Construction
#-----------------#

# base on the grid dots of df, create a meash that connects the dots onto a larger net
  # so as to enable the integration later

str(df_Cmpt_Res_log_16[c("Lon", "Lat")])

mesh <- inla.mesh.2d(loc = df_Cmpt_Res_log_16[c("Lon", "Lat")],
             cutoff = 0,
             max.edge = 0.75,
             offset = 4)

str(mesh)


#-------------------#
# mesh locations
#-------------------#

str(mesh$loc)  # num [1:218377, 1:3]
mesh_locs <- mesh$loc[, 1:2]

str(mesh_locs) # num [1:218377, 1:2]


#---------------------------#
# Distance btw pairs of locs
#---------------------------#

## The grid dots (Lon, Lat) we created are based on degree so earth ref system
## need to transform it as a Cartesian system such that the process can be approximated on a 2D planar 

## use Gneiting(2010) -- greate-circle distance
  # the unit is "miles"
  # The function RFearth2cartesian returns a matrix in one-to-one correspondence with coord assuming that the earth is an ellipsoid.
  # The function RFearth2dist calculates distances, cf. dist, assuming that the earth is an ellipsoid.


D <- as.matrix(RFearth2dist(coord = as.matrix(mesh_locs))) # dist in Cartisian sys

D_vec <- as.double(c(D))

rm(Rst_Cmpt_Res_log_16)
rm(WHO_map)
rm(WHO_regions)
rm(fitted_BC_Mean_log)
rm(fit_BC_log)
rm(fitted_BC_Residual_log)
rm(lm_fit_Cmpts)
rm(dfr)
rm(Bk_Res_log_16)
rm(df_all)
rm(df_all_log)
rm(df_all_long_log)
rm(df_Cmpts)
rm(df_Cmpts_Mn_Res_long_log)
rm(df_Cmpts_Mn_Res_wide_log)
rm(form_Compts)
rm(plt)
rm(Z_BC_OM)
rm(df_Cmpts_Res_wide_log)
rm(labels_expback_total)

