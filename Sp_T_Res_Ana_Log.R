#******************************#
# Residual Analysis (log scale)
#******************************#

str(df_Cmpts_Res_long_log)
# 'data.frame':	821520 obs. of  7 variables:
#$ Lon     : num  -42.8 -42 -39.8 -38.2 -37.5 ...
#$ Lat     : num  83.2 83.2 83.2 83.2 83.2 ...
#$ Year    : num  2016 2016 2016 2016 2016 ...
#$ ID      : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Cmpts   : chr  "BC" "BC" "BC" "BC" ...
#$ Values  : num  -4.84 -4.98 -5.01 -5.01 -4.97 ...
#$ Residual: num  0.307 0.305 0.301 0.299 0.298 ...

install.packages("fields") # for dist()
library(fields)

install.packages("ape") # for moran's I
library(ape)

install.packages("lmtest") # for DW test
library(lmtest)

library(dplyr)
library(tidyr)
library(ggplot2)
library(STRbook)


#========================================#
# Temporal dependence: Durbin-Watson Test
#========================================#

# Test residual's temp depd at every location for 6 componets

install.packages('tidyr')
library(tidyr)

install.packages("purrr")
library(purrr)

install.packages("broom")
library(broom)


#-----------------------------------------#
# Test on the residuals after global detrend
#-------------------------------------------#

dwtest_one_station <- function(data) {
  dwtest(Residual ~ 1, data = data)
}


DW_Cmpts <- list()
nested_df_Cmpts <- list()
df_Cmpts<- list()
for(i in unique(df_Cmpts_Res_long_log$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long_log, Cmpts == i)
  
  nested_df_Cmpts[[i]] <- group_by(df_Cmpts[[i]], Lon, Lat) %>% nest()
  
  DW_Cmpts[[i]] <- nested_df_Cmpts[[i]] %>% 
    mutate(dwtest = purrr::map(data, dwtest_one_station)) %>%
    mutate(test_df = purrr::map(dwtest, tidy)) %>%
    unnest(test_df)
  
}



#~~~~~~~~~~#
# P-values
#~~~~~~~~~~#

DW_results_Cmpts <- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  DW_results_Cmpts[[i]] <- mean(DW_Cmpts[[i]]$p.value < (0.05 / nrow(DW_Cmpts[[i]]))) * 100
}


DW_results_Cmpts_log_2 <- list()
for(i in unique(df_Cmpts_Res_long_log$Cmpts)) {
  DW_results_Cmpts_log_2[[i]] <- mean(DW_Cmpts[[i]]$p.value < 0.05 ) * 100
}

# $BC
#[1] 7.5409

#$DU
#[1] 5.466696

#$OM
#[1] 6.668127

#$SS
#[1] 5.306018

#$SU
#[1] 8.442886

#$PM25
#[1] 6.56953



#----------------------------#
# Test on the original Values 
#----------------------------#

# test on values after removing a constant temporal trend
dwtest_one_station_2 <- function(data) {
  dwtest(Values ~ 1, data = data)
}


DW_Cmpts_log_2 <- list()
nested_df_Cmpts <- list()
df_Cmpts<- list()
for(i in unique(df_Cmpts_Res_long_log$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long_log, Cmpts == i)
  
  nested_df_Cmpts[[i]] <- group_by(df_Cmpts[[i]], Lon, Lat) %>% nest()
  
  DW_Cmpts_log_2[[i]] <- nested_df_Cmpts[[i]] %>% 
    mutate(dwtest2 = purrr::map(data, dwtest_one_station_2)) %>%
    mutate(test_df2 = purrr::map(dwtest2, tidy)) %>%
    unnest(test_df2)
  
}


DW_results_Cmpts_log_Val <- list() # test on original values (log)
for(i in unique(df_Cmpts_Res_long_log$Cmpts)) {
  DW_results_Cmpts_log_Val[[i]] <- mean(DW_Cmpts_log_2[[i]]$p.value < 0.05 ) * 100
}



# $BC
#[1] 7.97546

#$DU
#[1] 5.514169

#$OM
#[1] 6.836109

#$SS
#[1] 5.229331

#$SU
#[1] 8.578002

#$PM25
#[1] 6.595092





