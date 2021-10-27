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



#======================#
# df_Cmpts_Res_long_log
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



#-------
# Wide
#------

df_Cmpts_Res_wide_log <- df_Cmpts_Res_long_log %>% 
  select(Lon, Lat, Year, Cmpts, Residual) %>%
  spread(key = Cmpts, value = Residual)

head(df_Cmpts_Res_wide_log)

