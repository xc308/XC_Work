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
head(TS_Cmpts_log)


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



#==================#
# Empirical Tmp Mean
#==================#

head(TS_Cmpts_log)
nrow(TS_Cmpts_log)

Year <- rep((2016:2012), 6)
length(Year)

Tp_Mean_df <- data.frame(TS_Cmpts_log, Year, stringsAsFactors = F)
head(Tp_Mean_df)


P_tp <- list()
plot_tp_df <- list()
for (i in unique(Tp_Mean_df$Label)) {
  plot_tp_df[[i]] <- filter(Tp_Mean_df, Labels == i)
  
  P_tp[[i]] <- ggplot(plot_tp_df[[i]]) + 
    geom_line(aes(x = Year, y = as.numeric(Values)),
              colour = "LightSkyBlue") +
    xlab("Year") + 
    ylab("Temporal Mean \n (log)") + 
    facet_wrap(~Labels) +
    theme_light()
} 

do.call(grid.arrange, P_tp)




#=======================#
# Empirical Spatial Mean
#=======================#

Sp_mean_Cmpt_log <- list()
for(i in unique(df_all_long_log$Cmpts)) {
  Sp_mean_Cmpt_log[[i]] <- filter(df_all_long_log, Cmpts == i) %>%
    select(Lon, Lat, Values) %>%
    group_by(Lon, Lat) %>%
    nest() %>%
    mutate(SpMean = map(data, colMeans)) %>%
    ungroup()
}

str(Sp_mean_Cmpt_log)
head(Sp_mean_Cmpt_log[[1]])
nrow(Sp_mean_Cmpt_log[[1]]) # 27384

Sp_mean_Cmpt_log_long <- do.call("rbind", Sp_mean_Cmpt_log)
str(Sp_mean_Cmpt_log_long)
# $ SpMean:List of 164304


## add labels
Labs <-  c("BC", "DU", "OM", "SS", "SU", "PM25")

Labels <- rep(Labs, each = 27384)
str(Labels) # chr [1:164304] "BC" "BC" "BC" "B


SpMean_Cmpts_long_df <- data.frame(Sp_mean_Cmpt_log_long, Labels, stringsAsFactors = F)
head(SpMean_Cmpts_long_df, 10)


SpMean_Cmpts_long_df2 <- SpMean_Cmpts_long_df %>% select(-data)
head(SpMean_Cmpts_long_df2)
tail(SpMean_Cmpts_long_df2)

rownames(SpMean_Cmpts_long_df2) <- NULL



#======#
# Plot 
#======#

# to see how empirical mean varies according to the lon, lat
# Ref: https://stackoverflow.com/questions/54373807/multipanel-plot-with-ggplot2

install.packages("grid")
install.packages("gridExtra")
install.packages("lattice")
install.packages("latticeExtra")
install.packages("cowplot")

library(grid)
library(gridExtra)
library(lattice)
library(latticeExtra)
library(cowplot)


#--------------#
# Latitude wise
#--------------#

P_lat <- list()
plot_df <- list()
for (i in unique(SpMean_Cmpts_long_df2$Labels)) {
  plot_df[[i]] <- filter(SpMean_Cmpts_long_df2, Labels == i)
  
  P_lat[[i]] <- ggplot(plot_df[[i]]) + 
    geom_point(aes(x = Lat, y = as.numeric(SpMean)),
               colour = "LightSkyBlue", alpha = 0.15) +
    xlab("Lattitude (deg)") + 
    ylab("Empirical Mean") + 
    facet_wrap(~ Labels) +
    theme_light()
} 

do.call(grid.arrange, P_lat) # plot


#--------------#
# Longitude wise
#--------------#

P_lat <- list()
plot_df <- list()
for (i in unique(SpMean_Cmpts_long_df2$Labels)) {
  plot_df[[i]] <- filter(SpMean_Cmpts_long_df2, Labels == i)
  
  P_lat[[i]] <- ggplot(plot_df[[i]]) + 
    geom_point(aes(x = Lon, y = as.numeric(SpMean)),
               colour = "LightSkyBlue", alpha = 0.15) +
  #geom_point(aes(x = Lon, y = as.numeric(SpMean)),
   #          colour = "#b3daff", alpha = 0.15) +
    xlab("Longitude (deg)") + 
    ylab("Empirical Mean") + 
    facet_wrap(~ Labels) +
    theme_light()
} 

do.call(grid.arrange, P_lat) # plot



