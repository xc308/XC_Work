#**********************#
# Sp_T_Exp_Covariance 2
#**********************#


str(df_Cmpts_Res_long)

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

save(X, file = "5yr_Residual_matrix_all_Cmpts")



#==============================#
# cross-covariance/correlation
#==============================#

# cross-correlation btw pairs of Cmpts
# lag-0 only

str(X[["BC"]])

## BC VS DU
Lag0_corr_BC_DU <- cor(X[["BC"]], X[["DU"]])

## DU VS SU
Lag0_corr_DU_SU <- cor(X[["DU"]], X[["SU"]])

## SU VS PM25
Lag0_corr_SU_PM25 <- cor(X[["SU"]], X[["PM25"]])

## PM VS OM
Lag0_corr_PM25_OM <- cor(X[["PM25"]], X[["OM"]])

## OM VS SS
Lag0_corr_OM_SS <- cor(X[["OM"]], X[["SS"]])

save(Lag0_corr_BC_DU, Lag0_corr_DU_SU, Lag0_corr_SU_PM25, 
     Lag0_corr_PM25_OM, Lag0_corr_OM_SS, file = "Lag0_cross_corr.RData")



