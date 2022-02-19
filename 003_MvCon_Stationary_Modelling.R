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
rm(mes)

#-------------------#
# mesh locations
#-------------------#

str(mesh$loc)  # num [1:218377, 1:3]
mesh_locs <- mesh$loc[, 1:2]

str(mesh_locs) # num [1:218377, 1:2]  # ellipsoid

save(mesh_locs, file = "/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/mesh_locs.rda")


#---------------------------#
# Distance btw pairs of locs
#---------------------------#

## The grid dots (Lon, Lat) we created are based on degree so earth ref system
## need to transform it as a Cartesian system such that the process can be approximated on a 2D planar 

## use Gneiting(2010) -- greate-circle distance
  # the unit is "miles"
  # The function RFearth2cartesian returns a matrix in one-to-one correspondence with coord assuming that the earth is an ellipsoid.
  # The function RFearth2dist calculates distances, cf. dist, assuming that the earth is an ellipsoid.


## Practical problem:
  ## 
      # The matrix size it can create is limited by .Machine$integer.max [1] 2147483647
      sqrt(2147483647) 
      # so the dimension of the mesh_loc matrix is capped at 46340 * 46340

 
      
.Machine$integer.max      
           


mesh_locs_cart <- as.matrix(RFearth2cartesian(coord = as.matrix(mesh_locs)))
str(mesh_locs_cart)  # num [1:218377, 1:3]

M_loc_cart <- mesh_locs_cart[, 1:2]   # cartisian
str(M_loc_cart) # num [1:218377, 1:2] 

dist(M_loc_cart)  # fail




#### History ####


## Try using sparse matrix

# https://stackoverflow.com/questions/5171593/r-memory-management-cannot-allocate-vector-of-size-n-mb

str(mesh_locs) # num [1:218377, 1:2]
mesh_locs_spars <- mesh_locs %>% as("dgeMatrix")
str(mesh_locs_spars)

try <- as.matrix(mesh_locs_spars )
rm(try)
rm(mesh_locs_spars)
## Not useful



## try bigmemory package (Warning: R cannot save the file)
install.packages("bigmemory")
library(bigmemory)
#detach("package:bigmemory")

bg_mesh_locs <- as.big.matrix(mesh_locs)
str(bg_mesh_locs)
object.size(bg_mesh_locs)


## try snow package
install.packages("snow")
library(snow)


cl <- makeCluster(1)
D <- parRapply(cl, x = as.matrix(mesh_locs), fun = RFearth2dist)
str(D)



## To calculate the Euclidean distance for large matrix
install.packages("distances")
library(distances)


D <- distances(M_loc_cart)
str(D)  # 'distances' num [1:2, 1:218377]

# turn distances obj into dist obj (fail)
D_dist <- distance_matrix(D)



## Try to calculate the distance by myself
n1 <- n2 <- nrow(M_loc_cart)

h <- matrix(0, n1 * n2, 2)   
for (i in 1:n2) {
  h[((i - 1) * n1) : (i * n1), ] <- t(t(mesh_locs_cart) - mesh_locs_cart[i, ])
}

# Error in matrix(0, n1 * n2, 2) : invalid 'nrow' value (too large or NA)
#In addition: Warning message:
#  In n1 * n2 : NAs produced by integer overflow
