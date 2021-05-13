#*****************#
# Residual Analysis
#*****************#

str(df_Cmpts_Res_long)

# 'data.frame':	821520 obs. of  7 variables:
#$ Lon     : num  -42.8 -42 -39.8 -38.2 -37.5 ...
#$ Lat     : num  83.2 83.2 83.2 83.2 83.2 ...
#$ Year    : num  2016 2016 2016 2016 2016 ...
#$ ID      : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Cmpts   : chr  "BC" "BC" "BC" "BC" ...
#$ Values  : num  0.0079 0.00689 0.00667 0.0067 0.00693 ...
#$ Residual: num  0.307 0.305 0.301 0.299 0.298 ...

install.packages("fields") # for dist()
library(fields)

install.packages("ape") # for moran's I
library(ape)

install.packages("lmtest") # for DW test
library(lmtest)



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



dwtest_one_station <- function(data) {
  dwtest(Residual ~ 1, data = data)
}


DW_Cmpts <- list()
nested_df_Cmpts <- list()
df_Cmpts<- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long, Cmpts == i)
  
  nested_df_Cmpts[[i]] <- group_by(df_Cmpts[[i]], Lon, Lat) %>% nest()
  
  DW_Cmpts[[i]] <- nested_df_Cmpts[[i]] %>% 
    mutate(dwtest = purrr::map(data, dwtest_one_station)) %>%
    mutate(test_df = purrr::map(dwtest, tidy)) %>%
    unnest(test_df)
  
}

DW_Cmpts[["BC"]]
# # A tibble: 27,384 x 8
# Groups:   Lon, Lat [27,384]
# Lon   Lat data        dwtest statistic p.value method     alternative         
#<dbl> <dbl> <list>      <list>     <dbl>   <dbl> <chr>      <chr>               
#  1 -42.8  83.2 <tibble[,5… <htes…      1.95   0.474 Durbin-Wa… true autocorrelatio…
#  2 -42    83.2 <tibble[,5… <htes…      2.07   0.533 Durbin-Wa… true autocorrelatio…

DW_results_Cmpts <- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  DW_results_Cmpts[[i]] <- mean(DW_Cmpts[[i]]$p.value < (0.05 / nrow(DW_Cmpts[[i]]))) * 100
}


save(DW_results_Cmpts, file = "DW_test_results_Cmpts")



head(DW_Cmpts[["BC"]]$p.value < 0.05 / nrow(DW_Cmpts[["BC"]]))
# [1] FALSE FALSE FALSE FALSE FALSE FALSE

mean(DW_Cmpts[["BC"]]$p.value < 0.05 / nrow(DW_Cmpts[["BC"]])) * 100
