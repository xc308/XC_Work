#***********************#
# Grid-wise Tmp Detrend
#***********************#

library(dplyr)
library(tidyr)

library(broom)
library(purrr)


head(df_all)
str(df_all) # 136920 = 27384 * 5

mod <- list()
Cmpts <- c("BC", "DU", "OM", "SS", "SU", "PM25")
for(Cmpt in Cmpts) {
  mod[[Cmpt]] <- as.formula(paste0(Cmpt, " ~ ", "1 + Year"))
  
}


res_one_grid <- function(model) {
  residuals(model)
}


Tmpmod_one_grid_lst <- list()
for (Cmpt in Cmpts) {
  Tmpmod_one_grid_lst[[Cmpt]] <- function(data, form) {
    lm(form, data = data)
  }
}


grd_Tmp_fit <- list()
Cmpts <- c("BC", "DU", "OM", "SS", "SU", "PM25")
for(Cmpt in Cmpts){
  grd_Tmp_fit[[Cmpt]] <- dplyr::select(df_all, Lon, Lat, Cmpt, Year) %>% group_by(Lon, Lat) %>% nest() %>%
    mutate(Model = map(data, Tmpmod_one_grid_lst[[Cmpt]], form = mod[[Cmpt]])) %>%
    mutate(Res = map(Model, res_one_grid)) %>%
    mutate(Mean_Res = map(Res, mean)) 
    
}
## Note: map(data, Tmpmod_one_grid_lst[[Cmpt]], form = mod[[Cmpt]])
  # supply the arg when call the function, as looping index won't loop over function body

save(grd_Tmp_fit, file = "Grid-wise_Tmp_detrend_fit_res.RData")


head(grd_Tmp_fit[["BC"]])
# # A tibble: 27,384 x 5
# Groups:   Lon, Lat [27,384]
#Lon   Lat data                 Model  Res      
#<dbl> <dbl> <list>               <list> <list>   
#  1 -42.8  83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]>
#  2 -42    83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]>
#  3 -39.8  83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]>  



## Unlist Mean_Res
for(Cmpt in Cmpts) {
  
  grd_Tmp_fit[[Cmpt]]$Mean_Res_unlst <- unlist(grd_Tmp_fit[[Cmpt]]$Mean_Res)
              
}

grd_Tmp_fit[['BC']]
# # A tibble: 27,384 x 7
# Groups:   Lon, Lat [27,384]
#     Lon   Lat data                 Model  Res       Mean_Res  Mean_Res_unlst
#<dbl> <dbl> <list>               <list> <list>    <list>             <dbl>
#  1 -42.8  83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]> <dbl [1]>       0       
#  2 -42    83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]> <dbl [1]>       2.60e-19
#  3 -39.8  83.2 <tibble[,2] [5 × 2]> <lm>   <dbl [5]> <dbl [1]>      -8.67e-20


## Ref
grd_Tmp_fit[["BC"]]$Mean_Res_unlst <- unlist(grd_Tmp_fit[["BC"]]$Mean_Res)
head(grd_Tmp_fit[["BC"]])
          


#=======================#
# Gridwise Tmp Residuals 
#=======================#


# a function to get Lon_Lat_Res
Lon_lat_res <- function(Cmpt) {
  grd_Tmp_fit[[Cmpt]] %>% select(Lon, Lat, Res) %>% unnest(Res)
}

res_df_lst <- list()
for(Cmpt in Cmpts) {
  res_df_lst[[Cmpt]] <- Lon_lat_res(Cmpt)
}

length(res_df_lst[[1]]$Lon) # 136920 = 27384 * 5
length(grd_Tmp_fit[["BC"]]$Lon) # 27384



Year <- rep(seq(2016, 2012), length(grd_Tmp_fit[["BC"]]$Lon))
length(Year) # 136920

for(Cmpt in Cmpts){
  res_df_lst[[Cmpt]]$Year <- Year
}


res_df_lst[["BC"]]
# # A tibble: 136,920 x 4
# Groups:   Lon, Lat [27,384]
#     Lon   Lat      Res  Year
#<dbl> <dbl>    <dbl> <int>
#1 -42.8  83.2 -0.00626  2016
#2 -42.8  83.2  0.00260  2015
#3 -42.8  83.2  0.0106   2014
#4 -42.8  83.2 -0.00400  2013
#5 -42.8  83.2 -0.00296  2012


Lon_Lat_res_yr_lst <-list()
for(Cmpt in Cmpts){
  Lon_Lat_res_yr_lst[[Cmpt]] <- res_df_lst[[Cmpt]] %>% ungroup()
}

head(Lon_Lat_res_yr_lst[["BC"]])
## A tibble: 6 x 4
#    Lon   Lat      Res  Year
#   <dbl> <dbl>    <dbl> <int>
#1 -42.8  83.2 -0.00626  2016
#2 -42.8  83.2  0.00260  2015
#3 -42.8  83.2  0.0106   2014

#head(grd_resd %>% ungroup())
# ref:https://stackoverflow.com/questions/38511743/adding-missing-grouping-variables-message-in-dplyr-in-r


save(Lon_Lat_res_yr_lst, file = "gd_ws_tmp_fit_res_yr.RData")



#=================#
# Space-wide resd
#=================#

Lon_Lat_res_yr_lst[["BC"]]

X_lst <- list()
for(Cmpt in Cmpts) {
  X_lst[[Cmpt]] <- Lon_Lat_res_yr_lst[[Cmpt]] %>% 
    spread(key = Year, value = Res) %>%
    select(-Lon, -Lat) %>%
    t()
}

dim(X_lst[["BC"]]) #  5 * 27384

save(X_lst, file = "Sp_wd_res_lst_grd_wise.RData")

## X_lst:
  # is grid-wise tmp detrended res in space-wide format 5 * 27384
  # can be used to calculate cov in file "Sp_T_Exp_Covariance_New"



##### can be done further in ACF, lag-tau cov######
           ## but just tried on BC #####
#============================================#
# ACF of grid-wise temporal detrended residual 
#============================================#

# group by year
# averaging over space grids 

TS_bc <- rowMeans(X_bc)
str(TS_bc)
acf(TS_bc)

# no auto-correlation left in the residuals of grid-wise temporally detrend



#============#
# Lag0_corr_bc
#============#

Lag0_cor_tp_detrend_bc <- cor(X_bc)
save(Lag0_cor_tp_detrend_bc, file = "Lag0_cor_tp_detrend_bc.RData")






#=======#
# ref
#=======#

a <- dplyr::select(df_all, Lon, Lat, "BC", Year) %>% group_by(Lon, Lat) %>% nest() %>%
  mutate(Model = map(data, Tmpmod_one_grid_lst[["BC"]], form = mod[["BC"]])) 
  

b <- dplyr::select(df_all, Lon, Lat, BC, Year) %>% group_by(Lon, Lat)
head(b, 2)  


mod <- list()
Cmpts <- c("BC", "DU", "OM", "SS", "SU", "PM25")
for(Cmpt in Cmpts) {
  mod[[Cmpt]] <- as.formula(paste0(Cmpt,  " ~ ", "1 + Year"))
  
}

mod[["BC"]] 

Tmpmod_one_grid_lst <- list()
#mod_lm_lst <- list()
for (Cmpt in Cmpts) {
  Tmpmod_one_grid_lst[[Cmpt]] <- function(data, form) {
      lm(form, data = data)
  }
}

Tmpmod_one_grid_lst[["BC"]]



#-----------#
# Suggestions
#-----------#


##
function(Cmpt, ){
  
  #Tmpmod_one_grid_lst
  function(data, form) {
    lm(form, data = data)
  }
  
  
  #mod
  
  
  dplyr::select(df_all, Lon, Lat, Cmpt, Year) %>% group_by(Lon, Lat) %>% nest() %>%
    mutate(Model = map(data, Tmpmod_one_grid_lst[[Cmpt]], form = mod[[Cmpt]])) 
}
