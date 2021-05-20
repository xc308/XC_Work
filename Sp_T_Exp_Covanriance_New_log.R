#**********************#
# Sp_T_Exp_Covariance 2 (log scale)
#**********************#

str(df_Cmpts_Res_long_log)
# Res get from stepwise selection ~ Lon + Lat + Lat^2

# 'data.frame':	821520 obs. of  7 variables:
#$ Lon     : num  -42.8 -42 -39.8 -38.2 -37.5 ...
#$ Lat     : num  83.2 83.2 83.2 83.2 83.2 ...
#$ Year    : num  2016 2016 2016 2016 2016 ...
#$ ID      : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Cmpts   : chr  "BC" "BC" "BC" "BC" ...
#$ Values  : num  0.0079 0.00689 0.00667 0.0067 0.00693 ...
#$ Residual: num  0.307 0.305 0.301 0.299 0.298 ...



#===========================#
# Get the spatial locations
#===========================#
install.packages("dplyr")
install.packages("tidyr")
library(dplyr)
library(tidyr)


sp_loc_df <- filter(df_Cmpts_Res_long_log, Cmpts == "BC") %>%
  filter(Year == 2016) %>%
  select(Lon, Lat) %>%
  arrange(Lon, Lat)


head(sp_loc_df) #
str(sp_loc_df)  # 'data.frame':	27384 obs. of  2 variables:
tail(sp_loc_df)



#=======================#
# Put df into space-wide
#=======================#

X_log <- list()
df_Cmpts <- list()
for(i in unique(df_Cmpts_Res_long_log$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long_log, Cmpts == i)
  
  X_log[[i]] <- select(df_Cmpts[[i]], Lon, Lat, Year, Residual) %>%
    spread(key = Year, value = Residual) %>%
    select(-Lon, -Lat) %>% 
    t()
  
}

str(X_log)



#=======================#
# Univariate correlation
#=======================#

# Uni-var
# Lag0 only

Lag0_corr_BC_log <- cor(X_log[["BC"]])
Lag0_corr_DU_log <- cor(X_log[["DU"]])
Lag0_corr_OM_log <- cor(X_log[["OM"]])
Lag0_corr_SU_log <- cor(X_log[["SU"]])
Lag0_corr_SS_log <- cor(X_log[["SS"]])
Lag0_corr_PM25_log <- cor(X_log[["PM25"]])

save(Lag0_corr_BC, Lag0_corr_DU, Lag0_corr_OM,
     Lag0_corr_SU, Lag0_corr_SS, Lag0_corr_PM25, 
     file = "Lag0_corr_uni.RData")



#==============================#
# cross-covariance/correlation
#==============================#

# cross-correlation btw pairs of Cmpts
# lag-0 only

str(X_log[["BC"]])

## BC VS DU
Lag0_corr_BC_DU_log <- cor(X_log[["BC"]], X_log[["DU"]])

## DU VS SU
Lag0_corr_DU_SU_log <- cor(X_log[["DU"]], X_log[["SU"]])

## SU VS PM25
Lag0_corr_SU_PM25_log <- cor(X_log[["SU"]], X_log[["PM25"]])

## PM VS OM
Lag0_corr_PM25_OM_log <- cor(X_log[["PM25"]], X_log[["OM"]])

## OM VS SS
Lag0_corr_OM_SS_log <- cor(X_log[["OM"]], X_log[["SS"]])


Lag0_cor <- c(Lag0_corr_BC_DU, Lag0_corr_DU_SU, Lag0_corr_SU_PM25, Lag0_corr_PM25_OM, Lag0_corr_OM_SS)
for(i in unique(Lag0_cor)) {
  save(i, file = paste0(i, ".RData"))
}




#=========================#
# Plot of cross/-covariance
#=========================#

#======#
# Plot
#======#

install.packages("fields")
install.packages("sp")
library(fields)
library(sp)


#----------------------------#
# Divid the Lon into 4 strips
#----------------------------#

sp_loc_df$n <- 1:nrow(sp_loc_df) # index each loaction
lim_lon <- range(sp_loc_df$Lon)  # -179.25  180.00
lon_strips <- seq(lim_lon[1], lim_lon[2], length.out = 5)

sp_loc_df$lon_strip <- cut(sp_loc_df$Lon, lon_strips, 
                           labels = FALSE, include.lowest = TRUE)

head(sp_loc_df)
tail(sp_loc_df)
nrow(sp_loc_df) # 27384



#---------------#
# Plot function
#---------------#

EMP_cor_lat_plt <- function(C, sp_loc_df) {
  require(fields)
  
  for (i in seq_along(unique(sp_loc_df$lon_strip))) {
    sp_strip <-  filter(sp_loc_df,  lon_strip  == i) %>% 
      arrange (Lat)
    
    idx <- sp_strip$n
    jitter <- seq(0, 1e-4, length = length(idx))
    
    image.plot(x = sp_strip$Lat + jitter, 
               y = sp_strip$Lat + jitter,
               z = C[idx, idx],
               zlim = c(-1, 1),
               #zlim = c(-0.1, 0.1),
               xlab = "Latitude", ylab = "Latitude",
               col = tim.colors(100), 
               #col = brewer_AC.div(500),
               cex = 200)
  }
}



#========================#
# Plot of uni-correlation
#========================#

png("Lag0_corr_BC-%02d.png")
png("Lag0_corr_BC-04.png")
EMP_cor_lat_plt(Lag0_corr_BC, sp_loc_df)

png("Lag0_corr_DU-%02d.png")
png("Lag0_corr_DU-04.png")
EMP_cor_lat_plt(Lag0_corr_DU, sp_loc_df)

png("Lag0_corr_OM-%02d.png")
png("Lag0_corr_OM-04.png")
EMP_cor_lat_plt(Lag0_corr_OM, sp_loc_df)

png("Lag0_corr_SU-%02d.png")
EMP_cor_lat_plt(Lag0_corr_SU, sp_loc_df)

png("Lag0_corr_SS-%02d.png")
EMP_cor_lat_plt(Lag0_corr_SS, sp_loc_df)

png("Lag0_corr_PM25-%02d.png")
EMP_cor_lat_plt(Lag0_corr_PM25, sp_loc_df)



#==========================#
# Plot of cross-correlation
#==========================#

#jpeg("Lag0_corr_BC_DU-%02d.jpeg")
png("Lag0_corr_BC_DU-%02d.png")
EMP_cor_lat_plt(Lag0_corr_BC_DU, sp_loc_df)
dev.off()

png("Lag0_corr_DU_SU-%02d.png")
EMP_cor_lat_plt(Lag0_corr_DU_SU, sp_loc_df)
dev.off()

png("Lag0_corr_SU_PM25-%02d.png")
EMP_cor_lat_plt(Lag0_corr_SU_PM25, sp_loc_df)
dev.off()

png("Lag0_corr_PM25_OM-%02d.png")
EMP_cor_lat_plt(Lag0_corr_PM25_OM, sp_loc_df)
dev.off()

png("Lag0_corr_OM_SS-%02d.png")
EMP_cor_lat_plt(Lag0_corr_OM_SS, sp_loc_df)
dev.off()



