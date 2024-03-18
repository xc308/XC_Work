#=============
# 18 Mar. 2024
#=============

# Aim:  
  # Prepare data for 2D TST12 and parallelization


# Data:
  # ref: Sp_T_Exp_detrend_global_log.R
  # ref: df_Cmpts_Res_long_log.RData
  
  # df_Cmpt_Res_log_16.rda


#=============
# Pre-settings
#=============
library(dplyr)

head(df_Cmpt_Res_log_16)
# Lon Lat Year BC_Residuals_log DU_Residuals_log OM_Residuals_log
  # SS_Residuals_log SU_Residuals_log PM25_Residuals_log


Coods_Res_log_16 <- cbind(df_Cmpt_Res_log_16$Lon, df_Cmpt_Res_log_16$Lat)

head(Coods_Res_log_16)

range(Coods_Res_log_16[, 2]) # Lat
# [1] -55.50  83.25

range(Coods_Res_log_16[, 1]) # Lon
# [1] -179.25  180.00



## Coods arrange rule: 
  # given a latitude, the longitude follows decreasing order
  # the latitude is also decreasing order

## so can divide into 4 groups (for parallel) by Lat
  # also can divide into 4 groups by Lon

  # 1st need to get cut and get the corresponding group index


#==================================
# cut the data into 4 groups by Lat
#==================================
# save to a new copy
df_Res_log_16 <- df_Cmpt_Res_log_16


# Define the latitude range and step
start_lat <- 83.25
end_lat <-  -55.50
step <- 0.75


# The number of division
num_div <- 4 

# division size
div_size = (end_lat - start_lat) / num_div
# [1] -34.6875

breaks = seq(start_lat, end_lat, length.out = num_div + 1)
length(breaks) # [1] 5
str(breaks)
# num [1:5] 83.2 48.6 13.9 -20.8 -55.5

diff(breaks)
# [1] -34.6875 -34.6875 -34.6875 -34.6775



df_Res_log_16 <- df_Res_log_16 %>%
  mutate(Lat_div_indx = cut(Lat, breaks = breaks, labels = F))


head(df_Res_log_16)
tail(df_Res_log_16)

# assigne Lat indx 1 to the last row that has the only indx catogray NA 
df_Res_log_16$Lat_div_indx[nrow(df_Res_log_16)] <- 1

df_div_Lat_lst <- split(df_Res_log_16, df_Res_log_16$Lat_div_indx)
str(df_div_Lat_lst)


#===================================
# cut the data into 4 groups by Lon
#===================================

diff(df_Res_log_16$Lon) # most are 0.75, but exists skip, ocean btw land

start_lon <- -179.25
end_lon <- 180

# The number of division
num_div <- 4 

# division size
div_size_lon = (end_lon - start_lon) / num_div
# [1] 89.8125

breaks_lon = seq(start_lon, end_lon, length.out = num_div + 1)
length(breaks_lon) # [1] 5
breaks_lon
# [1] -179.2500  -89.4375    0.3750   90.1875  180.0000


diff(breaks_lon)
# [1] 89.8125 89.8125 89.8125 89.8125


df_Res_log_16 <- df_Res_log_16 %>%
  mutate(Lon_div_indx = cut(Lon, breaks = breaks_lon, labels = F))


#--------------------------------
# check NA for the boundary cut
#--------------------------------

which(is.na(df_Res_log_16$Lon_div_indx)) 
# [1] 1772 2782 3114 3474

df_Res_log_16[1770:1779,]
df_Res_log_16[2780:2784, ]
df_Res_log_16[3112:3116, ]
df_Res_log_16[3472:3475,]
## all corresponds to Lon: -179.25
# so shall categorize to Lon_div_indx = 1

df_Res_log_16$Lon_div_indx[1772] <- 1
df_Res_log_16$Lon_div_indx[2782] <- 1
df_Res_log_16$Lon_div_indx[3114] <- 1
df_Res_log_16$Lon_div_indx[3474] <- 1

nrow(df_Res_log_16) # [1] 27384

df_div_Lon_lst <- split(df_Res_log_16, df_Res_log_16$Lon_div_indx)
str(df_div_Lon_lst)

# List of 4
# $ 1:'data.frame':	3793 obs. of  11 variables:
# $ 2:'data.frame':	6581 obs. of  11 variables:
# $ 3:'data.frame':	10155 obs. of  11 variables:
# $ 4:'data.frame':	6855 obs. of  11 variables:



#-----
# sort
#-----

# Within each Lon list, sort the data frame 
  # Lon from small to large (ascending), -180 ~ -90, -90~0, 0~90, 90~180
  # while within each Lon, Lat is sorted from large to small (desc) # 83.25 ~ 50, 50~0,0~-50


df_div_Lon_lst_sorted <- list()
for (i in 1:length(df_div_Lon_lst)) {
  
  df_div_Lon_lst_sorted[[i]] <- arrange(df_div_Lon_lst[[i]], Lon, desc(Lat))
}

str(df_div_Lon_lst_sorted)
# List of 4
#$ :'data.frame':	3793 obs. of  11 variables:
# $ :'data.frame':	3793 obs. of  11 variables:
#..$ Lon               : num [1:3793] -179 -179 -179 -179 -178 ...
#..$ Lat               : num [1:3793] 71.2 68.2 67.5 66.8 71.2 ...
#..$ Lon_div_indx      : int [1:3793] 1 1 1 1 1 1 1 1 1 1 ...

#$ :'data.frame':	6581 obs. of  11 variables:
#  ..$ Lon               : num [1:6581] -89.2 -89.2 -89.2 -89.2 -89.2 ...
#  ..$ Lat               : num [1:6581] 81.8 80.2 79.5 78.8 76.5 ...
#  ..$ Lon_div_indx      : int [1:6581] 2 2 2 2 2 2 2 2 2 2 ...

# $ :'data.frame':	10155 obs. of  11 variables:
# ..$ Lon               : num [1:10155] 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 ...
# ..$ Lat               : num [1:10155] 52.5 51.8 51 49.5 48.8 ...
# ..$ Lon_div_indx      : int [1:10155] 3 3 3 3 3 3 3 3 3 3 ...

# $ :'data.frame':	6855 obs. of  11 variables:
# ..$ Lon               : num [1:6855] 90.8 90.8 90.8 90.8 90.8 ...
# ..$ Lat               : num [1:6855] 75 74.2 73.5 72.8 72 ...
# ..$ Lon_div_indx      : int [1:6855] 4 4 4 4 4 4 4 4 4 4 ...


df_Res_log_16_sorted <- do.call(rbind, df_div_Lon_lst_sorted)
str(df_Res_log_16_sorted)
# 'data.frame':	27384 obs. of  11 variables:
#$ Lon               : num  -178 -178 -178 -178 -178 ...
#$ Lat               : num  71.2 68.2 67.5 66.8 66 ...
#$ Year              : num  2016 2016 2016 2016 2016 ...
#$ BC_Residuals_log  : num  -0.835 -0.304 -0.139 -0.236 -0.527 ...
#$ DU_Residuals_log  : num  0.338 0.854 1.063 0.892 0.43 ...
#$ OM_Residuals_log  : num  -0.6494 -0.1025 0.0297 -0.1527 -0.4518 ...
#$ SS_Residuals_log  : num  1.28 1.01 1.07 1.19 1.43 ...
#$ SU_Residuals_log  : num  -0.3935 -0.1554 0.0187 -0.0218 -0.2455 ...
#$ PM25_Residuals_log: num  0.442 0.414 0.499 0.425 0.374 ...
#$ Lat_div_indx      : num  4 4 4 4 4 4 4 4 4 4 ...
#$ Lon_div_indx      : int  1 1 1 1 1 1 1 1 1 1 ...



#=======================================================
# Re-organize the BC,DU,etc.order to align with the DAG
#=======================================================

# 1: DU; 2: SU; 3: BC; 4: OM; 5:SS


df_Res_log_16_sorted_DAG_V1 <- data.frame(Lon = df_Res_log_16_sorted$Lon,
           Lat = df_Res_log_16_sorted$Lat, 
           Year = df_Res_log_16_sorted$Year, 
           Z1 = df_Res_log_16_sorted$DU_Residuals_log, 
           Z2 = df_Res_log_16_sorted$SU_Residuals_log, 
           Z3 = df_Res_log_16_sorted$BC_Residuals_log, 
           Z4 = df_Res_log_16_sorted$OM_Residuals_log, 
           Z5 = df_Res_log_16_sorted$SS_Residuals_log, 
           PM25_res = df_Res_log_16_sorted$PM25_Residuals_log, 
           Lat_div_indx = df_Res_log_16_sorted$Lat_div_indx,
           Lon_div_indx = df_Res_log_16_sorted$Lon_div_indx)




#======
# Save
#======

saveRDS(df_Res_log_16_sorted, file = "df_Res_log_16_sorted.rds")
#try <- readRDS("df_Res_log_16_sorted.rds")

#head(try)

saveRDS(df_Res_log_16_sorted_DAG_V1, file = "df_Res_log_16_sorted_DAG_V1.rds")





