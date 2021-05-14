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

DW_results_Cmpts_2 <- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  DW_results_Cmpts_2[[i]] <- mean(DW_Cmpts[[i]]$p.value < 0.05 ) * 100
}


save(DW_results_Cmpts, file = "DW_test_results_Cmpts.RData")


load(DW_test_results_Cmpts)

head(DW_Cmpts[["BC"]]$p.value < 0.05 / nrow(DW_Cmpts[["BC"]]))
# [1] FALSE FALSE FALSE FALSE FALSE FALSE

mean(DW_Cmpts[["BC"]]$p.value < 0.05 / nrow(DW_Cmpts[["BC"]])) * 100

mean(DW_Cmpts[["BC"]]$p.value < 0.05) * 100 
mean(DW_Cmpts[["OM"]]$p.value < 0.05) * 100 

head(DW_Cmpts[["BC"]]$p.value)
hist(DW_Cmpts[["BC"]]$p.value)

#=====================#
# SP-T semivariogram
#=====================#

install.packages("sp")
library(sp)

install.packages("spacetime")
library(spacetime)

install.packages("gstat")
library(gstat)

library(dplyr)

#-------#
# STFDF
#-------#

head(df_Cmpts_Res_long)

STFDF_Cmpts <- list()
df_Cmpts <- list()
for(i in unique(df_Cmpts_Res_long$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_Cmpts_Res_long, Cmpts == i)
  
  STFDF_Cmpts[[i]] <- STFDF(sp = SpatialPoints(coords = unique(data.frame(df_Cmpts[[i]]$Lon, df_Cmpts[[i]]$Lat))),
        time = xts(unique(df_Cmpts[[i]]$Year), order.by = as.Date(as.character(unique(df_Cmpts[[i]]$Year)), format = "%Y")),
        data = select(df_Cmpts[[i]], -Lon, -Lat, -Year))
  
}

rownames(df_Cmpts_Res_long) <- NULL

head(df_Cmpts_Res_long)


sp_part <- SpatialPoints(coords = unique(data.frame(df_Cmpts_bc$Lon, df_Cmpts_bc$Lat)))
str(sp_part) # @ coords     : num [1:27384, 1:2]


# ref: https://stackoverflow.com/questions/23224142/converting-data-frame-to-xts-order-by-requires-an-appropriate-time-based-object
install.packages("xts")
library(xts)

t_part <- xts(unique(df_Cmpts_bc$Year), order.by = as.Date(as.character(unique(df_Cmpts_bc$Year)), format = "%Y"))
str(t_part)  # Data: num [1:5, 1] 2012 2013 2014 2015 2016


str(select(df_Cmpts_bc, -Lon, -Lat, -Year)) # data.frame':	136920 obs. of  4 variables

df_Cmpts_bc <- filter(df_Cmpts_Res_long, Cmpts == "BC")
STFDF_bc <- STFDF(sp = sp_part, 
                  time = t_part, 
                  data = select(df_Cmpts_bc, Residual))

class(STFDF_bc)

vv_bc <- variogram(object = Residual ~ 1,
          data = STFDF_bc,
          width = 80,
          cutoff = 1000,
          tlags = 0:4)

plot(vv_bc)

ts <- xts(df_Cmpts_bc$Year, order.by = as.Date(as.character(df_Cmpts_bc$Year), format = "%Y"))
str(ts) 
# An ‘xts’ object on 2012-05-13/2016-05-13 containing:
# Data: num [1:136920, 1] 2012 2012 2012 2012 2012 ...
# Indexed by objects of class: [Date] TZ: UTC
unique(ts) # 2012 2013 2014 2015 2016


var_one_yr <- function(data) {
  variogram(Residual ~ 1, 
            data = data)
}

str(df_bc_yr_nest)

df_bc_yr_nest <- group_by(df_Cmpts_bc, Year) %>% nest()
head(df_bc_yr_nest)

df_bc_yr_nest %>%
  map(var_one_yr)


