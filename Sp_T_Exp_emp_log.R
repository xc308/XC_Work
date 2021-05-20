#********************************************#
# Space-time Exploratory Analysis (on log scale)
#********************************************#

#===============#
# Empirical ACF
#===============#

str(df_Cmpts_Res_long_log)
# 'data.frame':	821520 obs. of  7 variables:
#$ Lon     : num  -42.8 -42 -39.8 -38.2 -37.5 ...
#$ Lat     : num  83.2 83.2 83.2 83.2 83.2 ...
#$ Year    : num  2016 2016 2016 2016 2016 ...
#$ ID      : int  1 2 3 4 5 6 7 8 9 10 ...
#$ Cmpts   : chr  "BC" "BC" "BC" "BC" ...
#$ Values  : num  -4.84 -4.98 -5.01 -5.01 -4.97 ...
#$ Residual: num  0.307 0.305 0.301 0.299 0.298 ...

df_all_long_log <- df_Cmpts_Res_long_log %>% select(-Residual)
str(df_all_long_log )
# 'data.frame':	821520 obs. of  6 variables:



#---------#
# 6 Cmpts
#---------#

df_Cmpts <- list()
year_mean_Cmpt_log <- list()
for(i in unique(df_all_long_log$Cmpts)) {
  df_Cmpts[[i]] <- filter(df_all_long_log, Cmpts == i)
  
  year_mean_Cmpt_log[[i]] <- df_Cmpts[[i]] %>% 
    select(Year, Values) %>%
    group_by(Year) %>%
    nest() %>%
    mutate(YrMean = map(data, colMeans)) %>%
    ungroup()
  
}


YrMean_Cmpts_log <- list()
for(i in seq_along(year_mean_Cmpt_log)) {
  YrMean_Cmpts_log[[i]] <- do.call("rbind", year_mean_Cmpt_log[[i]]$YrMean)
}


save(YrMean_Cmpts_log, file = "YrMean_avg_all_grids_Cmpts_log.RData")


## add labels
Labs <-  c("BC", "DU", "OM", "SS", "SU", "PM25")

Labels <- rep(Labs, each = 5)
str(Labels)

str(YrMean_Cmpts_log)
YrMean_Cmpts_long_log <- do.call("rbind", YrMean_Cmpts_log)
str(YrMean_Cmpts_long_log) # num [1:30, 1] -1.39 -1.32 -1.28 -1.3 -1.16 ...
head(YrMean_Cmpts_long_log)


TS_Cmpts_log <- data.frame(YrMean_Cmpts_long_log, Labels, stringsAsFactors = F)
str(TS_Cmpts_log)



#-----------------#
# Plot ACF 6 Cmpts
#-----------------#

install.packages("ggfortify")
library(ggfortify)
library(ggplot2)

install.packages("gridExtra")
library(gridExtra)

library(dplyr)


Plt_ts_log <- list()
TS_log <- list()
for(i in unique(TS_Cmpts_log[,2])) {
  TS_log[[i]] <- filter(TS_Cmpts_log, Labels == i)
  Plt_ts_log[[i]] <- autoplot(acf(TS_log[[i]]$Values, plot = F), main = i)
  
}

do.call(grid.arrange, Plt_ts_log)




