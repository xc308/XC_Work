#load("/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/df_Cmpts_Res_long_log.RData")


#load("/Users/gavin/Desktop/XC_WORK/XC_WORK/df_Cmpts_Res_long_log.RData")

load("df_Cmpts_Res_long_log.RData")

df_Cmpts_Res_long <- df_Cmpts_Res_long_log

rm(df_Cmpts_Res_long_log)

#===========================#
# Get the spatial locations
#===========================#
#install.packages("dplyr")
#install.packages("tidyr")
library(dplyr)
library(tidyr)


sp_loc_df <- filter(df_Cmpts_Res_long, Cmpts == "BC") %>%
  filter(Year == 2016) %>%
  select(Lon, Lat) %>%
  arrange(Lon, Lat)


head(sp_loc_df) #
str(sp_loc_df)  # 'data.frame':	27384 obs. of  2 variables:
tail(sp_loc_df)


#=======================#
# Put df into space-wide
#=======================#

X <- list()
df_Cmpts <- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long, Cmpts == i)
  
  X[[i]] <- select(df_Cmpts[[i]], Lon, Lat, Year, Residual) %>%
    spread(key = Year, value = Residual) %>%
    select(-Lon, -Lat) %>%
    t()
  
}

str(X)

save(X, file = "5yr_Residual_matrix_all_Cmpts.rda")


#install.packages("fields")
#install.packages("sp")
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

EMP_cor_lat_plt2 <- function(C, sp_loc_df, strip_num = NA) {
  
  #for (i in seq_along(unique(sp_loc_df$lon_strip))) {
    sp_strip <-  filter(sp_loc_df,  lon_strip  == strip_num) %>% 
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
               cex = 200, 
               cex.lab = 2,
               cex.axis = 2)
  #}
}
