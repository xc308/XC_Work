#********#
# Recap 
#********#

## Continue working from 24 May Sp_T_Exp_detrend_global_log.R

getwd()
# /Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/XC_Work
save(df_all, file = "/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/df_all.rda")
load(file = "/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/df_all.rda")


#==========#
# Packages
#==========#

install.packages("dplyr")
install.packages("tidyr")


library(dplyr)
library(tidyr)



#=====================#
# df_all & df_all_log
#=====================#

head(df_all)
# at each give location, the sum of 5 Cmpts is PM25
# 0.007900851 + 0.003339927 + 0.11029157 + 0.05850084+ 0.01452385
# = 0.194557


log(5) + log(3) + log(7)  # 4.65396
log(15)  # 2.70805




df_all_lg <- df_all
head(df_all_lg)


df_all_log <- df_all_lg %>% 
  mutate(across(.col = c(BC:PM25), .fns = log))

rm(df_all_lg)
head(df_all_log)
str(df_all_log)
# 'data.frame':	136920 obs. of  10 variables:
  #$ Lon : num  -42.8 -42 -39.8 -38.2 -37.5 ...
  #$ Lat : num  83.2 83.2 83.2 83.2 83.2 ...
  #$ BC  : num  -4.84 -4.98 -5.01 -5.01 -4.97 ...
  #$ DU  : num  -5.7 -5.86 -5.8 -5.95 -5.97 ...
  #$ OM  : num  -2.2 -2.34 -2.38 -2.39 -2.37 ...
  #$ SS  : num  -2.84 -2.96 -2.97 -3.15 -3.23 ...
  #$ SU  : num  -4.23 -4.35 -4.36 -4.37 -4.35 ...
  #$ PM25: num  -1.64 -1.77 -1.8 -1.86 -1.86 ...
  #$ Year: num  2016 2016 2016 2016 2016 ...
  #$ ID  : int  1 2 3 4 5 6 7 8 9 10 ...

# 5 yrs data 27384 * 5 = 136920


#------------#
# fit BC lm
#------------#
fit_BC_log <- lm(BC ~ Lon + Lat + I(Lat^2), data = df_all_log)

fitted_BC_Mean_log <- fit_BC_log$fitted.values
fitted_BC_Residual_log <- fit_BC_log$residuals

str(fitted_BC_Mean_log)
# Named num [1:136920] -3.47 -3.47 -3.47 -3.46 
str(fitted_BC_Residual_log)
# Named num [1:136920] -1.37 -1.51 -1.54 -1.54 -1.51 

#rm(fitted_BC_Mean)
#rm(fitted_BC_Residual)



#======================#
# df_Cmpts_Res_long_log (spot a mistake)
#======================#

load("/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/df_Cmpts_Res_long_log.RData")

str(df_Cmpts_Res_long_log)
# 'data.frame':	821520 obs. of  7 variables:
  #$ Lon     : num  -42.8 -42 -39.8 -38.2 -37.5 ...
  #$ Lat     : num  83.2 83.2 83.2 83.2 83.2 ...
  #$ Year    : num  2016 2016 2016 2016 2016 ...
  #$ ID      : int  1 2 3 4 5 6 7 8 9 10 ...
  #$ Cmpts   : chr  "BC" "BC" "BC" "BC" ...
  #$ Values  : num  -4.84 -4.98 -5.01 -5.01 -4.97 ...
  #$ Residual: num  0.307 0.305 0.301 0.299 0.298 ...

# 821520 / 6 = 136920 = 27384 * 5yr

head(df_Cmpts_Res_long_log)
# Residual: after global detrending i, "~", "Lon + Lat + I(Lat ^ 2)"

tail(df_Cmpts_Res_long_log)



#--------------------#
# Detrending using lm
#--------------------#



form_Compts <- list()
for (i in unique(df_all_long_log$Cmpts)) {
  form_Compts <- filter(df_Cmpts_Res_long_log, Cmpts == i) %>%
    as.formula(paste0(i, "~", "Lon + Lat + I(Lat^2)")) %>%
  lm_fit_cmpts <- lm(formula = form_Compts[[i]], data = df_all_log)
  
}






df_Cmpts <-list()
for(i in unique(df_all_long_log$Cmpts)) {
  
  df_Cmpts[[i]] <- filter(df_all_long_log, Cmpts == i)
  form_cmpts[[i]] <- as.formula(paste0(i, "~", "Lon + Lat + I(Lat ^ 2)"))
  lm_fit_cmpts[[i]] <- lm(form_cmpts[[i]], data = df_all_log)
  df_Cmpts[[i]]$Residual <- lm_fit_cmpts[[i]]$residuals
  
}

df_Cmpts_Res_long_log <- do.call("rbind", df_Cmpts) # 821520 * 7 = 136920 * 6 components * 7 var




#-------
# Wide
#------

df_Cmpts_Res_wide_log <- df_Cmpts_Res_long_log %>% 
  select(Lon, Lat, Year, Cmpts, Residual) %>%
  spread(key = Cmpts, value = Residual)

head(df_Cmpts_Res_wide_log)

save(df_Cmpts_Res_wide_log, file = "/Users/xchen/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/RData/df_Cmpts_Res_wide_log.RData")



range(df_Cmpts_Res_wide_log$BC)  # [1] -1.050601 53.809133
quantile(df_Cmpts_Res_wide_log$BC)
#         0%         25%         50%         75%        100% 
# -1.05060097 -0.44997465 -0.19906984  0.09433469 53.80913259 

quantile(df_Cmpts_Res_wide_log$DU)
#        0%        25%        50%        75%       100% 
# -21.240140 -10.537969  -5.738078   2.681979 953.412098 

quantile(df_Cmpts_Res_wide_log$)



