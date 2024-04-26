rm(list = ls())
getwd()

#===================================#
# Code Chunk 0: Preliminaries (2016)
#===================================#

load_path <- "~/OneDrive - University of Exeter/XC_PhD/Data/Processed/"

#load_path <- "~/Downloads/XC_PhD/Data/Processed/"

load_files <- c("bc/", "du/", "om/", "ss/", "su/", "pm25/")

load_names <- c("bc_rean_2016", "du_rean_2016", "om_rean_2016",
                "ss_rean_2016", "su_rean_2016", "pm25_rean_2016")  

load_geo_names <- c("shapefiles", "WHO_lines", "WHO_regions") 

full_path <- paste0(load_path, load_files, load_names, ".RData", sep = "")
full_path[1]

full_path_geo <- paste0(load_path, load_geo_names, ".RData", sep = "")


for (i in seq_along(full_path)) {
  load(full_path[i])
}


install.packages("sp")
yes
library(sp)

for (i in seq_along(full_path_geo)) {
  load(full_path_geo[i])
}

getwd()
wk_path <- "~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/"

#wk_path <- "~/Downloads/XC_PhD/Data/Processed/XC_WORK/"


setwd(wk_path)
setwd(paste0(wk_path, "Graph_preliminary", sep = ""))

#load("Data_for_chunk(0-6).RData")
#load("Data_for_all_map_chunks.RData")


#=================================#
# Code Chunk 1: change data units 
#=================================#

install.packages("stringr")
library(stringr)


comp_names <- lapply(load_names, function(x) str_to_upper(x))

comp_data <- list(bc_rean_2016, du_rean_2016, om_rean_2016,
                  ss_rean_2016, su_rean_2016, pm25_rean_2016)

for (i in seq_along(comp_data)) {
  assign(comp_names[[i]], lapply(comp_data, function(x) x * 1e9)[[i]]) 
}


assign(comp_names[[1]], lapply(comp_data, function(x) x * 1e9)[[1]]) 
comp_names[[1]]



#==============================================#
# Code Chunk 2: Range of Each Scaled Component 
#==============================================#

data_rescaled <- list(BC_REAN_2016, DU_REAN_2016, OM_REAN_2016,
                      SS_REAN_2016, SU_REAN_2016, PM25_REAN_2016)

data_range <- lapply(data_rescaled, range)

range_name <- paste0(comp_names, "_range", sep = "")


for (i in seq_along(range_name)) {
  assign(range_name[i], data_range[[i]])
  #print(c(range_name[i], data_range[[i]]))
}


BC_REAN_2016_range # [1] 9.273717e-04 2.098888e+01
DU_REAN_2016_range # [1] 7.934834e-12 7.955341e+02
OM_REAN_2016_range # [1]   0.04656196 216.79821210
SS_REAN_2016_range # [1]  0.002575607 22.977516924
SU_REAN_2016_range # [1]  0.005663878 85.693526563
PM25_REAN_2016_range # [1]   0.08783205 797.70192111

Range <- data.frame(BC_REAN_2016_range, DU_REAN_2016_range,
                    OM_REAN_2016_range, SS_REAN_2016_range,
                    SU_REAN_2016_range, PM25_REAN_2016_range)

head(Range)



#=================================#
# Code Chunk 3: Rasterize Data
#=================================#
install.packages("raster")
yes
library(raster)

# resolution 
lon_re <- seq(from = 0, to = 359.25, by = 0.75) # 480
lat_re <- seq(from = 90, to = -90, by = -0.75) # 241

# think  CAMS does 0->360,  map this to -180->180
lon_re[lon_re > 180] <- lon_re[lon_re > 180] - 360

head(lon_re) # [1] 0.00 0.75 1.50 2.25 3.00 3.75
head(order(lon_re)) # [1] 242 243 244 245 246 247


# raster framework
#R_PM25 <- R_SU <- R_SS <- R_OM <- R_DU <- R_BC <- raster(ncol = length(lon_re), nrow = length(lat_re))

# GS try
#PM25_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(pm25_rean_2016[order(lon_re), ]))
#BC_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(bc_rean_2016[order(lon_re), ]))
#DUST_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(du_rean_2016[order(lon_re), ]))
#OM_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(om_rean_2016[order(lon_re), ]))
#SS_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(ss_rean_2016[order(lon_re), ]))
#SU_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(su_rean_2016[order(lon_re), ]))


R_PM25 <- R_SU <- R_SS <- R_OM <- R_DU <- R_BC <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75)



#------------
# Try new res
#------------

#R_BC2 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 7.5)
#R_BC2
# class      : RasterLayer 
#dimensions : 24, 48, 1152  (nrow, ncol, ncell)
#resolution : 7.5, 7.5  (x, y)
#extent     : -179.625, 180.375, -89.625, 90.375  (xmin, xmax, ymin, ymax)
#crs        : NA 



# GS
#PM25_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(pm25_rean_2016[order(lon_re), ]))
#BC_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(bc_rean_2016[order(lon_re), ]))
#DUST_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(du_rean_2016[order(lon_re), ]))
#OM_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(om_rean_2016[order(lon_re), ]))
#SS_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(ss_rean_2016[order(lon_re), ]))
#SU_2016 <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75, vals = t(su_rean_2016[order(lon_re), ]))



R_BC
#class      : RasterLayer 
#dimensions : 241, 480, 115680  (nrow, ncol, ncell)
#resolution : 0.75, 0.75  (x, y)
#extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
#crs        : NA

dim(BC_REAN_2016) # [1] 480 241
dim(t(BC_REAN_2016)) # [1] 241 480


# want to assign values to each raster layer
#values(R_BC) <- t(BC_REAN_2016)

#Data_trans <- list(t(BC_REAN_2016), t(DU_REAN_2016), t(OM_REAN_2016),
#  t(SS_REAN_2016), t(SU_REAN_2016), t(PM25_REAN_2016))

Data_trans <- list(t(BC_REAN_2016[order(lon_re), ]), t(DU_REAN_2016[order(lon_re), ]), t(OM_REAN_2016[order(lon_re), ]),
                   t(SS_REAN_2016[order(lon_re), ]), t(SU_REAN_2016[order(lon_re), ]), t(PM25_REAN_2016[order(lon_re), ]))


dim(Data_trans[[1]]) #  [1] 241 480

R <- c(R_BC, R_DU, R_OM, R_SS, R_SU, R_PM25) # list
for (i in seq_along(R)) {
  values(R[[i]]) <- Data_trans[[i]]
}

names(R) <- c("R_BC", "R_DU", "R_OM", "R_SS", "R_SU", "R_PM25")
R



#=================================#
# Code Chunk 4: Trim Raster Values
#=================================#
WHO_map # SpatialPolygonsDataFrame

R_mask <- lapply(R, function(x) mask(x, WHO_map))

hasValues(R_mask[[1]]) # [1] TRUE
inMemory(R_mask[[1]]) # [1] TRUE

R_mask[[1]]
# class      : RasterLayer 
# dimensions : 241, 480, 115680  (nrow, ncol, ncell)
# resolution : 0.75, 0.75  (x, y)
# extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : memory
# names      : layer 
# values     : 0.003380918, 20.98888  (min, max)

#~~~~~~~~OLD lon_res~~~~~~~~~~~~~~~~
# class      : RasterLayer 
#dimensions : 241, 480, 115680  (nrow, ncol, ncell)
#resolution : 0.75, 0.746888  (x, y)
#extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#source     : memory
#names      : layer 
#values     : 0.003862587, 16.7525  (min, max)
#~~~~~~~~~~~~~~~


#=================================#
# Code Chunk 5: Name Raster Layers
#=================================#
#Layer_bc <- R_mask[[1]]
#sqrt(Layer_bc)

Layer_names  <- paste("Layer_", c("bc", "du", "om", "ss", "su", "pm25"), sep = "")

layer_list <- vector("list", length = length(Layer_names))
for (i in seq_along(Layer_names)) {
  layer_list[[i]] <- assign(Layer_names[i], R_mask[[i]])
}

names(layer_list) <- Layer_names



#============================================#
# Code Chunk 6: Sqrt, cubrt, log Layer Values
#============================================#

#-----#
#sqrt
#-----#

Layer_names_sqrt <- paste("Layer_", c("bc", "du", "om", "ss", "su", "pm25"), 
                          "_sqrt", sep = "")

layer_sqrt_list <- vector("list", length = length(Layer_names_sqrt))
for (i in seq_along(Layer_names_sqrt)) {
  layer_sqrt_list[[i]] <- assign(Layer_names_sqrt[i], lapply(layer_list, sqrt)[[i]])
}

names(layer_sqrt_list) <- Layer_names_sqrt


#----------#
# cubic rt
#----------#

Layer_names_cbrt <- paste("Layer_", c("bc", "du", "om", "ss", "su", "pm25"), 
                          "_cbrt", sep = "")

layer_cbrt_list <- vector("list", length = length(Layer_names_cbrt))
for (i in seq_along(Layer_names_cbrt)) {
  layer_cbrt_list[[i]] <- assign(Layer_names_cbrt[i], 
                                 lapply(layer_list, function(x) x ^ (1/3))[[i]])
}


names(layer_cbrt_list) <- Layer_names_cbrt

#27 ^ (1/3)


#------#
# log
#------#

Layer_names_log <- paste("Layer_", c("bc", "du", "om", "ss", "su", "pm25"), 
                         "_log", sep = "")

layer_log_list <- vector("list", length = length(Layer_names_log))
for (i in seq_along(Layer_names_log)) {
  layer_log_list[[i]] <- assign(Layer_names_log[i], lapply(layer_list, log)[[i]])
}

names(layer_log_list) <- Layer_names_log



#log(795.5341)
#exp(6.679014)



#============================================#
# Code Chunk 7-1: Brick Sqrt Raster Layers List
#============================================#

#-----------------#
# Total sqrt Brick
#-----------------#

B_sqrt <- brick(layer_sqrt_list[[1]], layer_sqrt_list[[2]],
                layer_sqrt_list[[3]], layer_sqrt_list[[4]],
                layer_sqrt_list[[5]], layer_sqrt_list[[6]])


names(B_sqrt) <- Layer_names_sqrt
B_sqrt

#class      : RasterBrick 
#dimensions : 241, 480, 115680, 6  (nrow, ncol, ncell, nlayers)
#resolution : 0.75, 0.746888  (x, y)
#extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#source     : memory
#names      : Layer_bc_sqrt, Layer_du_sqrt, Layer_om_sqrt, Layer_ss_sqrt, Layer_su_sqrt, Layer_pm25_sqrt 
#min values :    0.06214971,    0.01386099,    0.32680529,    0.12974081,    0.10282345,      0.41162294 
#max values :      4.092981,     28.205213,     14.724069,      4.509582,      9.257080,       28.243617 


#----------------#
# 1, 4 sqrt Brick
#----------------#
B_sqrt_14 <- brick(layer_sqrt_list[[1]], layer_sqrt_list[[4]])
names(B_sqrt_14) <- c("BC", "SS")

cellStats(B_sqrt_14, 'range')
#              BC        SS
# [1,] 0.06214971 0.1297408
# [2,] 4.09298126 4.5095817


#-----------------#
# 5, 3 sqrt Brick
#------------------#

B_sqrt_53 <- brick(layer_sqrt_list[[5]], layer_sqrt_list[[3]])
names(B_sqrt_53) <- c("SU", "OM") 

cellStats(B_sqrt_53, 'range')
#             SU         OM
# [1,] 0.1028234  0.3268053
# [2,] 9.2570798 14.7240691


#----------------#
# 2, 6 sqrt Brick
#----------------#

B_sqrt_26 <- brick(layer_sqrt_list[[2]], layer_sqrt_list[[6]])
names(B_sqrt_26) <- c("DU", "PM25")

cellStats(B_sqrt_26, 'range')

#               DU       PM25
# [1,]  0.01386099  0.4116229
# [2,] 28.20521322 28.2436174




#=============================================#
# Code Chunk 7-2: Brick Cbrt Raster Layers List
#=============================================#

B_cbrt <- brick(layer_cbrt_list[[1]], layer_cbrt_list[[2]],
                layer_cbrt_list[[3]], layer_cbrt_list[[4]],
                layer_cbrt_list[[5]], layer_cbrt_list[[6]])


names(B_cbrt) <- Layer_names_cbrt


#----------------#
# 1, 4 cbrt Brick
#----------------#
B_cbrt_14 <- brick(layer_cbrt_list[[1]], layer_cbrt_list[[4]])
names(B_cbrt_14) <- c("BC", "SS")

cellStats(B_cbrt_14, 'range')

#       BC        SS
#[1,] 0.1569011 0.2562818
#[2,] 2.5587420 2.7295486


#-----------------#
# 5, 3 cbrt Brick
#------------------#

B_cbrt_53 <- brick(layer_cbrt_list[[5]], layer_cbrt_list[[3]])
names(B_cbrt_53) <- c("SU", "OM") 

cellStats(B_cbrt_53, 'range')

#            SU        OM
# [1,] 0.2194799 0.4744525
# [2,] 4.4087554 6.0073818


#----------------#
# 2, 6 cbrt Brick
#----------------#

B_cbrt_26 <- brick(layer_cbrt_list[[2]], layer_cbrt_list[[6]])
names(B_cbrt_26) <- c("DU", "PM25")

cellStats(B_cbrt_26, 'range')

#             DU      PM25
#[1,] 0.05770269 0.5533497
#[2,] 9.26587118 9.2742802



#=============================================#
# Code Chunk 7-3: Brick log Raster Layers List
#=============================================#

#B_bc_log <- brick(layer_log_list[[1]])
#cellStats(B_bc_log, 'range')
#           [,1]
# [1,] -5.556418
# [2,]  2.818547

#-------------#
# B log total
#-------------#

B_log <- brick(layer_log_list[[1]], layer_log_list[[2]],
               layer_log_list[[3]], layer_log_list[[4]],
               layer_log_list[[5]], layer_log_list[[6]])

names(B_log) <- c("BC", "DU", "OM", "SS", "SU", "PM25")


#-----------#
# B log 14
#----------#

B_log_14 <- brick(layer_log_list[[1]], layer_log_list[[4]])

names(B_log_14) <- c("BC", "SS")

cellStats(B_log_14, 'range')
#             BC        SS
# [1,] -5.556418 -4.084433
# [2,]  2.818547  3.012409


#----------#
# B log 53
#----------#

B_log_53 <- brick(layer_log_list[[5]], layer_log_list[[3]])

names(B_log_53) <- c("SU", "OM") 

cellStats(B_log_53, 'range')
#             SU        OM
# [1,] -4.549484 -2.236781
# [2,]  4.450777  5.378967


#----------#
# B log 26
#----------#

B_log_26 <- brick(layer_log_list[[2]], layer_log_list[[6]])

names(B_log_26) <- c("DU", "PM25")

cellStats(B_log_26, 'range')

#            DU      PM25
# [1,] -8.557354 -1.775295
# [2,]  6.679014  6.681735



#==================================#
# Code Chunk 8: Lattice-Style Plot
#==================================#

#install.packages(c("latticeExtra", "lattice","rasterVis"))
#install.packages("RColorBrewer")

library(lattice)
library(latticeExtra)
library(rasterVis)
library(RColorBrewer)

getwd()


#-------------#
# chose color
#-------------#

cols <- rev(brewer.pal(11, name = "RdBu"))
brewer.div <- colorRampPalette(cols, interpolate = "spline")

cols_ryb <- rev(brewer.pal(11, name = "RdYlBu"))
brewer.div.ryb <- colorRampPalette(cols_ryb, interpolate = "spline")

col_yor <- rev(brewer.pal(9, "YlOrRd"))
brewer.div.yor <- colorRampPalette(col_yor, interpolate = "spline")



#========================================#
# Code Chunk 8-1: Plot Block on log scale
#========================================#

#---------------------#
# plot log block total
#---------------------#
pretty(c(-4, 4))

tick_at_log <- c(-2.5, -1.5, -1, 0, 1, 2, 3, 4, 5, 6.5)
exp(tick_at_log)
# [1]   0.0820850   0.2231302   0.3678794
# [4]   1.0000000   2.7182818   7.3890561
# [7]  20.0855369  54.5981500   148.413  665.1416330

tick_at_log <- c(-4, -2, 0, 2, 3, 4, 5, 6.5)
exp(tick_at_log)
# [1]   0.01831564   0.13533528   1.00000000
# [4]   7.38905610  20.08553692  54.59815003
# [7] 148.41315910 665.14163304

#labels_expback_total <- c(0.01, 0.14, 1, 7, 20, 55, 148, 665)
labels_expback_total <- c(0.01, 1, 2, 7, 20, 55, 200, 700)

plt <- levelplot(B_log, cuts = 499,
                 col.regions = brewer.div(500),
                 layout = c(2, 3),
                 colorkey = list(labels = list(labels = labels_expback_total,
                                               cex = 1),
                                 width = 0.7))

plt + latticeExtra::layer(sp.polygons(WHO_map, col = "black", lwd = 0.05))



plt_ryb <- levelplot(B_log, cuts = 499,
                     col.regions = brewer.div.ryb(500),
                     layout = c(2, 3), 
                     colorkey = list(labels = list(labels = labels_expback_total),
                                     width = 0.7))

plt + latticeExtra::layer(sp.polygons(WHO_map, col = "black", lwd = 0.05))
plt_ryb + latticeExtra::layer(sp.polygons(WHO_map, col = "black", lwd = 0.05))


#--------------#
# plot B log 14
#--------------#

cellStats(B_log_14, 'range')
par(mfrow = c(1, 1))
hist(log(BC_REAN_2016))
hist(log(SS_REAN_2016))


targ_tick <- c(0.005, 0.01, 0.05, 0.15, 0.35, 0.5, 1, 2.5, 7.5, 20.5)
log(targ_tick)
log(20.5)

tick_at_log_14 <- c(-5.2983174, -4.6051702, -2.9957323, 
                    -1.8971200, -1.0498221, -0.6931472,
                    0.0000000, 0.9162907,
                    2.0149030, 3.0204249)
exp(tick_at_log_14)
labels_expback_14 <- c(0.005, 0.01, 0.05, 0.15, 0.35, 0.5,
                       1, 3, 7, 20)

levelplot(B_log_14, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), 
          colorkey = list(labels = list(labels = labels_expback_14),
                          width = 0.7,
                          space = "bottom"))



#--------------#
# plot B log 53
#--------------#

cellStats(B_log_53, 'range')
log(85.7)  #[1] 4.450853
log(216.8) #[1] 5.378975
exp(2.5) # [1] 12.18249

hist(log(SU_REAN_2016))
hist(log(OM_REAN_2016))


tick_at_log_53 <-  c(0, 2.5, 4, 4.5, 5, 5.4)
exp(tick_at_log_53)


levelplot(B_log_53, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), 
          colorkey = list(labels = list(labels = round(exp(tick_at_log_53)))))


#--------------#
# plot B log 26
#--------------#

cellStats(B_log_26, 'range')
hist(log(DU_REAN_2016))
hist(log(PM25_REAN_2016))

pretty(c(-8.5, 6))

tick_at_log_26 <- c(-5, -4, -2, 0, 3, 4, 5, 6.5)


#tick_at_log_26 <- c(0.0, 1.5, 2.5, 3.0, 4.0, 4.5, 5.5, 6.5)
exp(tick_at_log_26)

exp(-5) # [1] 0.006737947
exp(1)
exp(1.5) # [1] 4.481689
exp(2) # [1] 7.389056
exp(2.5) # [1] 12.18249
exp(4) # [1] 54.59815
exp(3) # [1] 20.08554

labels_expback <- c(0.01, 0.02, 0.14, 1, 20, 55, 148, 665)


levelplot(B_log_26, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), 
          colorkey = list(labels = list(labels = labels_expback)))


#------------------------#
# try PM25 independently
#------------------------#
seq(0, 7, by = 0.5)
tick_at_log_6 <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.5)
exp(tick_at_log_6)

R_pm_log <- layer_log_list[[6]]
names(R_pm_log) <- "PM25"

levelplot(R_pm_log, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), margin = F,
          colorkey = list(labels = list(labels = round(exp(tick_at_log_6)))))


cellFromXY(R_pm_log, c(30, 0)) # 57881 afric
exp(extract(R_pm_log, 57881: 57888))

hist(log(PM25_REAN_2016))



#========================================#
# Code Chunk 8-2: Plot Block on cub scale
#========================================#


#---------------------#
# Plot cub total block
#---------------------#
# not distinct enough
tick_at_cub <- pretty(c(0, 9), 10)


levelplot(B_cbrt, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 3), 
          colorkey = list(labels = list(labels = tick_at_cub ^ 3)))



#-------------------#
# plot block cbrt 14
#------------------#

tick_at_14_cb <- pretty(c(0, 3), 6)

levelplot(B_cbrt_14, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), 
          colorkey = list(labels = list(labels = tick_at_14_cb ^ 3)))


#-------------------#
# plot block cbrt53
#-------------------#

tick_at_53_cb <- pretty(c(0, 6), 5)

levelplot(B_cbrt_53, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1),
          colorkey = list(labels = list(labels = tick_at_53_cb ^ 3)))




#------------------#
# plot Block cbrt26
#------------------#

tick_at_26_cb <- pretty(c(0, 10), 5)

levelplot(B_cbrt_26, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1),
          colorkey = list(labels = list(labels = tick_at_26_cb ^ 3)))



#========================================#
# Code Chunk 8-3: Plot Block on sqrt scale
#========================================#

#-------------------#
# plot block sqrt 14
#------------------#
tick_at_14 <- pretty(c(0, 5), 5)

levelplot(B_sqrt_14, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1), 
          colorkey = list(labels = list(labels = tick_at_14 ^ 2)))



#-------------------#
# plot block sqrt53
#-------------------#

tick_at_53 <- pretty(c(0, 15), 10)

levelplot(B_sqrt_53, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1),
          colorkey = list(labels = list(labels = tick_at_53 ^ 2)))



#------------------#
# plot Block sqrt26
#------------------#

tick_at_26 <- pretty(c(0, 30), 10)

levelplot(B_sqrt_26, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 1),
          colorkey = list(labels = list(labels = tick_at_26 ^ 2)))



#-----------------------#
# Plot total brick sqrt
#----------------------#
# (not distinct enough)
tick_at <- pretty(c(0, 30))

levelplot(B_sqrt, cuts = 499, 
          col.regions = brewer.div.grn(500),
          layout = c(2, 3),
          colorkey = list(labels = list(labels = tick_at ^ 2)))



#---------------------#
# Plot 1345 sqrt block
#---------------------#
# still not distince enough 
tick_at_1345 <- pretty(c(0, 15), 10)

levelplot(B_sqrt_1345, cuts = 499,
          col.regions = brewer.div(500),
          layout = c(2, 2),
          colorkey = list(labels = list(labels = tick_at_1345 ^ 2)))


levelplot(B_sqrt_1345, cuts = 499,
          col.regions = brewer.div.yor(500),
          layout = c(2, 2),
          colorkey = list(labels = list(labels = tick_at_1345 ^ 2)))




