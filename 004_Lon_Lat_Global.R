#************************************
# create the lon lat for all location
#************************************
# ref: Data_process_for_df_all.R
# ref: XC_Map_GS.R

library(tidyr)
library(dplyr)


Lon <- seq(-179.625, 180.37, by = 0.75)
Lat <- seq(-90.375, 90.375, by = 0.75)

coords_glob <- expand.grid(Lon = Lon, Lat = Lat)
head(coords_glob)

coords_glob %>%
  filter(Lat == -54.75)



Cmpts <- c("BC", "DU", "OM", "SS", "SU")
CC <- expand.grid(Cmpts, Cmpts)





