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
          


#------#
# ref
#------#

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
