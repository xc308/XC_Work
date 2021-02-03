#========================================#
# Code Chunk 0: Preliminaries (2012-2016)
#========================================#

load_path <- "~/OneDrive - University of Exeter/XC_PhD/Data/Processed/"

#load_path <- "~/Downloads/XC_PhD/Data/Processed/"

load_files <- c("bc/", "du/", "om/", "ss/", "su/", "pm25/")

load_names <- c("bc_rean_", "du_rean_", "om_rean_",
                "ss_rean_", "su_rean_", "pm25_rean_")  

load_years <- c("2016", "2015", "2014", "2013", "2012")

full_path <- vector("list", length = length(load_years))
for (i in seq_along(load_years)) {
  full_path[[i]] <- paste0(load_path, load_files, load_names, load_years[i], ".RData", sep = "")
}

full_path[[1]]

full_path <- unlist(full_path) # vector 1:30


#full_path <- paste0(load_path, load_files, load_names, ".RData", sep = "")
#full_path[1]


load_geo_names <- c("shapefiles", "WHO_lines", "WHO_regions") 
full_path_geo <- paste0(load_path, load_geo_names, ".RData", sep = "")


for (i in seq_along(full_path)) {
  load(full_path[i])
}


#install.packages("sp")
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

install.packages("parallel")
library(parallel)


load_names_yr <- vector("list", length = length(load_years))
for (i in seq_along(load_years)) {
  load_names_yr[[i]] <- paste(load_names, load_years[i], sep = "")
}

#load_names_yr[[1]]

comp_names <- lapply(load_names_yr, function(x) str_to_upper(x))
comp_names
s <- comp_names[[1]] # a chr vector 1:6

comp_data <- list(list(bc_rean_2016, du_rean_2016, om_rean_2016,
                  ss_rean_2016, su_rean_2016, pm25_rean_2016),
                  list(bc_rean_2015, du_rean_2015, om_rean_2015,
                  ss_rean_2015, su_rean_2015, pm25_rean_2015),
                  list(bc_rean_2014, du_rean_2014, om_rean_2014,
                  ss_rean_2014, su_rean_2014, pm25_rean_2014),
                  list(bc_rean_2013, du_rean_2013, om_rean_2013,
                  ss_rean_2013, su_rean_2013, pm25_rean_2013),
                  list(bc_rean_2012, du_rean_2012, om_rean_2012,
                  ss_rean_2012, su_rean_2012, pm25_rean_2012))
# a list of 5 lists of 6 matrices

comp_data[[1]] 
length(comp_data)  # 5
length(comp_data[[1]]) # 6


for (i in seq_along(1:5)) {
  for(j in seq_along(1:6))
    assign(comp_names[[i]][j], mclapply(comp_data[[i]], function(x) x * 1e9, mc.cores = 3)[[j]]) 
}


#assign(comp_names[[1]], lapply(C, function(x) x * 1e9)[[1]]) 
#comp_names[[1]]



#=================================#
# Code Chunk 2: Rasterize Data
#=================================#
#install.packages("raster")
library(raster)

#----------#
# resolution 
#----------#
lon_re <- seq(from = 0, to = 359.25, by = 0.75) # 480
lat_re <- seq(from = 90, to = -90, by = -0.75) # 241

# think  CAMS does 0->360,  map this to -180->180
lon_re[lon_re > 180] <- lon_re[lon_re > 180] - 360

head(lon_re) # [1] 0.00 0.75 1.50 2.25 3.00 3.75
head(order(lon_re)) # [1] 242 243 244 245 246 247


#---------------------#
# Define Raster Names
#---------------------#

rev(R_names)

R_names <- c("R_PM25", "R_SU", "R_SS", "R_OM", "R_DU", "R_BC")
R_yrs <- c("_2016", "_2015", "_2014", "_2013", "_2012")

R_all <- vector("list", length = length(R_yrs))
for (i in seq_along(R_yrs)) {
  R_all[[i]] <- paste(rev(R_names), R_yrs[i], sep = "")
}

R_all[[1]]
R_all


#------------#
# Rasterize 
#------------#
#a <- assign(R_all[[1]][1], raster(xmx = 180.375, xmn = -179.625, 
                       ymx = 90.375, ymn = -90.375, res = 0.75))

#a

for (i in seq_along(1:5)) {
  for (j in seq_along(1:6)) {
    assign(R_all[[i]][j], raster(xmx = 180.375, xmn = -179.625, 
                            ymx = 90.375, ymn = -90.375, res = 0.75))
  }
}


R_BC_2012
# class      : RasterLayer 
# dimensions : 241, 480, 115680  (nrow, ncol, ncell)
# resolution : 0.75, 0.75  (x, y)
# extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
# crs        : NA 



#----------------------------#
# Single year rasterize (ref)
#----------------------------#

#R_PM25 <- R_SU <- R_SS <- R_OM <- R_DU <- R_BC <- raster(xmx = 180.375, xmn = -179.625, ymx = 90.375, ymn = -90.375, res = 0.75)

R_BC
#class      : RasterLayer 
#dimensions : 241, 480, 115680  (nrow, ncol, ncell)
#resolution : 0.75, 0.75  (x, y)
#extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
#crs        : NA

dim(BC_REAN_2016) # [1] 480 241
dim(t(BC_REAN_2016)) # [1] 241 480



#======================================#
# Code Chunk 3: Assign Values to raster
#======================================#

# want to assign values to each raster layer
#values(R_BC) <- t(BC_REAN_2016)

#Data_trans <- list(t(BC_REAN_2016), t(DU_REAN_2016), t(OM_REAN_2016),
#  t(SS_REAN_2016), t(SU_REAN_2016), t(PM25_REAN_2016))

Data_all_trans <- list(list(t(BC_REAN_2016[order(lon_re), ]), t(DU_REAN_2016[order(lon_re), ]), 
                        t(OM_REAN_2016[order(lon_re), ]), t(SS_REAN_2016[order(lon_re), ]), 
                        t(SU_REAN_2016[order(lon_re), ]), t(PM25_REAN_2016[order(lon_re), ])),
                   list(t(BC_REAN_2015[order(lon_re), ]), t(DU_REAN_2015[order(lon_re), ]), 
                        t(OM_REAN_2015[order(lon_re), ]), t(SS_REAN_2015[order(lon_re), ]), 
                        t(SU_REAN_2015[order(lon_re), ]), t(PM25_REAN_2015[order(lon_re), ])),
                   list(t(BC_REAN_2014[order(lon_re), ]), t(DU_REAN_2014[order(lon_re), ]), 
                        t(OM_REAN_2014[order(lon_re), ]), t(SS_REAN_2014[order(lon_re), ]), 
                        t(SU_REAN_2014[order(lon_re), ]), t(PM25_REAN_2014[order(lon_re), ])),
                   list(t(BC_REAN_2013[order(lon_re), ]), t(DU_REAN_2013[order(lon_re), ]), 
                        t(OM_REAN_2013[order(lon_re), ]), t(SS_REAN_2013[order(lon_re), ]), 
                        t(SU_REAN_2013[order(lon_re), ]), t(PM25_REAN_2013[order(lon_re), ])),
                   list(t(BC_REAN_2012[order(lon_re), ]), t(DU_REAN_2012[order(lon_re), ]), 
                        t(OM_REAN_2012[order(lon_re), ]), t(SS_REAN_2012[order(lon_re), ]), 
                        t(SU_REAN_2012[order(lon_re), ]), t(PM25_REAN_2012[order(lon_re), ])))


length(Data_all_trans) # 5
length(Data_all_trans[[1]]) #  6
dim(Data_all_trans[[1]][[1]]) # [1] 241 480


R_comp_yr <- list(list(R_BC_2016, R_DU_2016, R_OM_2016, R_SS_2016, R_SU_2016, R_PM25_2016),
                  list(R_BC_2015, R_DU_2015, R_OM_2015, R_SS_2015, R_SU_2015, R_PM25_2015),
                  list(R_BC_2014, R_DU_2014, R_OM_2014, R_SS_2014, R_SU_2014, R_PM25_2014),
                  list(R_BC_2013, R_DU_2013, R_OM_2013, R_SS_2013, R_SU_2013, R_PM25_2013),
                  list(R_BC_2012, R_DU_2012, R_OM_2012, R_SS_2012, R_SU_2012, R_PM25_2012))

for (i in 1:5) { # year
  for (j in 1:6) { # each comp with year i
    values(R_comp_yr[[i]][[j]]) <- Data_all_trans[[i]][[j]]
  }
}


R_comp_yr[[1]][[1]]  # 2016 BC R
# class      : RasterLayer 
# dimensions : 241, 480, 115680  (nrow, ncol, ncell)
# resolution : 0.75, 0.75  (x, y)
# extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : memory
# names      : layer 
# values     : 0.0009273717, 20.98888  (min, max)


#------------------#
# Name each raster
#------------------#

for (i in 1:5) {
  for (j in 1:6) {
    names(R_comp_yr[[i]][[j]]) <- R_all[[i]][[j]]
  }
}

R_comp_yr[[1]] # all components in 2016


#------------------#
# single year (ref)
#------------------#

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

R_mask_all <- vector("list", length = 5) # 5 yrs
for (i in 1:5) {
  R_mask_all[[i]] <- mclapply(R_comp_yr[[i]], function(x) mask(x, WHO_map), mc.cores = 5)
}

R_mask_all[[1]][[1]]
# class      : RasterLayer 
# dimensions : 241, 480, 115680  (nrow, ncol, ncell)
# resolution : 0.75, 0.75  (x, y)
# extent     : -179.625, 180.375, -90.375, 90.375  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : memory
# names      : R_BC_2016 
# values     : 0.003380918, 20.98888  (min, max)


R_mask_all[[1]]



#=================================#
# Code Chunk 5: Stack Raster Layers
#=================================#
# by year

stack_i <- function(i = i) {
  stack(R_mask_all[[i]][[1]], R_mask_all[[i]][[2]],
        R_mask_all[[i]][[3]], R_mask_all[[i]][[4]],
        R_mask_all[[i]][[5]], R_mask_all[[i]][[6]])
}

R_mask_all[[1]]
stack_i(1)


stk_lst <- list(stk_2016 = stack_i(1), stk_2015 = stack_i(2), 
     stk_2014 = stack_i(3), stk_2013 = stack_i(4),
     stk_2012 = stack_i(5))
stk_lst

