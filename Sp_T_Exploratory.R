#********************************#
# Space-time Exploratory Analysis
#********************************#

getwd()
source("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_process_for_explore.RData")


install.packages("fields")
install.packages("CCA")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("gstat")
install.packages("sp")
install.packages("spacetime")

library(fields)
library(CCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gstat)  
library(sp)
library(spacetime) 


install.packages("devtools")
install.packages("usethis")
library(usethis)
library(devtools)

install_github("andrewzm/STRbook", force = TRUE)
library(STRbook)



#=======================#
# Empirical Spatial Mean
#=======================#

# Group by (lon, lat), averaging over years

#----------------------#
# Single Component Ref
#----------------------#

sp_av_BC <- df_all %>% select(Lon, Lat, BC, Year, ID) %>% group_by(Lat, Lon) %>% 
  summarise(mu_emp_BC = mean(BC))

count(sp_av_BC) # each Lat corresponds to how many different Lon, 186 different lat 
str(count(df_all, Lat, Lon)) # at each location, 5 datapoints averaging over
head(count(df_all, Lat, Lon))


#s <- df_all[c("Lon", "Lat", "BC", "Year", "ID")]


#-----------------------------#
# Select and group on the fly
#-----------------------------#

comps <-  c("BC", "DU", "OM", "SS", "SU", "PM25")

df_comps_lst <- vector("list", length = length(comps))
df_gp_lst <- vector("list", length = length(comps))
for (i in seq_along(comps)) {
  df_comps_lst[[i]] <-  df_all[c("Lon", "Lat", comps[i], "Year", "ID")]
  
  df_gp_lst[[i]] <- df_comps_lst[[i]] %>% group_by(Lat, Lon) 
    
}

df_gp_lst[[6]]


#--------------------------------#
# Calcuate spatial empirical mean
#--------------------------------#

emp_mu_BC <- df_gp_lst[[1]] %>% summarise(emp_mu_BC = mean(BC))
head(emp_mu_BC) # 27384 *3

emp_mu_DU <- df_gp_lst[[2]] %>% summarise(emp_mu_DU = mean(DU))
head(emp_mu_DU)

emp_mu_OM <- df_gp_lst[[3]] %>% summarise(emp_mu_OM = mean(OM))
head(emp_mu_OM)

emp_mu_SS <- df_gp_lst[[4]] %>% summarise(emp_mu_SS = mean(SS))
head(emp_mu_SS)

emp_mu_SU <- df_gp_lst[[5]] %>% summarise(emp_mu_SU = mean(SU))
head(emp_mu_SU)

emp_mu_PM25 <- df_gp_lst[[6]] %>% summarise(emp_mu_PM25 = mean(PM25))
head(emp_mu_PM25)


#----------#
# Question
#----------#
sp_av_comp <- function(df = df_all, comp = comp) {
  select(df, Lon, Lat, all_of(comp), Year, ID) %>% 
    group_by(Lat, Lon) %>% 
    summarise(mu_emp = mean(comp))
}

head(df_all)
str(df_all)

sp_av_BC <- sp_av_comp(df = df_all, comp = BC)

sp_av_lst <- for (i in seq_along(Comps)) {
  sp_av_comp(comp = Comps[[i]])
}

#-------------------#
# Wrong number why?
#-------------------#
Comps <- list(df_all$BC, df_all$DU, df_all$OM, df_all$SS,
              df_all$SU, df_all$PM25)

emp_mu_lst <- vector("list", length = length(Comps))
for (i in seq_along(Comps)) {
  emp_mu_lst[[i]] <- df_gp_lst[[i]] %>% summarise(mu_emp = mean(Comps[[i]]))
}

head(emp_mu_lst[[1]])



#==================#
# Add Label columns 
#==================#
# prepare for ggplot facet_wrap

nrow(emp_mu_BC) # [1] 27384

Lab_bc <- rep("BC", nrow(emp_mu_BC))
emp_mu_BC$Label <- rep("BC", nrow(emp_mu_BC))


Labs <- list("BC", "DU", "OM", "SS", "SU", "PM25")
emp_mu_lst <- list(emp_mu_BC, emp_mu_DU, emp_mu_OM,
                   emp_mu_SS, emp_mu_SU, emp_mu_PM25)

for (i in seq_along(Labs)) {
  emp_mu_lst[[i]]$Label <- rep(Labs[[i]], nrow(emp_mu_lst[[i]]))
}

head(emp_mu_lst[[1]])
head(emp_mu_lst[[2]])



#------------#
# Rename df
#-----------#

names(emp_mu_lst[[1]]) <- c("Lat", "Lon", "emp_mu", "Label") 

name_lst <- list(c("Lat", "Lon", "emp_mu", "Label"))

for(i in seq(1:6)) {
  names(emp_mu_lst[[i]]) <- name_lst[[1]]
}

head(emp_mu_lst[[1]])
#     Lat   Lon  emp_mu Label
#    <dbl> <dbl>   <dbl> <chr>
#  1 -55.5 -68.2 0.00321 BC   


#-----------------------#
# Row bind different df 
#-----------------------#

# to make long format

emp_mu_df_lng <- bind_rows(emp_mu_lst[[1]], emp_mu_lst[[2]],
          emp_mu_lst[[3]], emp_mu_lst[[4]],
          emp_mu_lst[[5]], emp_mu_lst[[6]])

# 164304 * 4
head(emp_mu_df_lng, 2)
tail(emp_mu_df_lng)



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
for (i in unique(emp_mu_df_lng$Label)) {
  plot_df[[i]] <- filter(emp_mu_df_lng, Label == i)
  
  P_lat[[i]] <- ggplot(plot_df[[i]]) + 
    geom_point(aes(x = Lat, y = emp_mu),
               colour = "LightSkyBlue", alpha = 0.15) +
    xlab("Lattitude (deg)") + 
    ylab("Empirical Mean") + 
    facet_wrap(~Label) +
    theme_light()
} 

do.call(grid.arrange, P_lat) # plot


#--------------#
# Longitude wise
#--------------#

P_lon <- list()
plot_df <- list()
for (i in unique(emp_mu_df_lng$Label)) {
  plot_df[[i]] <- filter(emp_mu_df_lng, Label == i)
  
  P_lon[[i]] <- ggplot(plot_df[[i]]) + 
    geom_point(aes(x = Lon, y = emp_mu),
               colour = "LightSkyBlue", alpha = 0.15) +
    xlab("Longitude (deg)") + 
    ylab("Empirical Mean") + 
    facet_wrap(~Label) +
    theme_light()
} 

do.call(grid.arrange, P_lon)


#-------------#
# Alternative
#-------------#

# NOT GOOD since BC, SS, SU's flunctuation is not distince due to the same y axes
P_lat_2 <- ggplot(emp_mu_df_lng) +
  geom_point(aes(x = Lat, y = emp_mu),
             colour = "SkyBlue", alpha = 0.04) + 
  facet_wrap(~Label) +
  xlab("Latitude (deg)") + 
  ylab("Empirical Mean") +
  theme_classic()
  
print(P_lat_2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~#



#========================#
# Empirical Temporal Mean 
# (Not very informative)
#========================#

head(df_all)

emp_tp_mu_BC <- df_all %>% select(Lon, Lat, BC, Year, ID) %>%
  group_by(Year) %>% summarise(emp_temp_mean = mean(BC))

emp_tp_mu_BC


#-------------------#
# Select on the fly
#-------------------#

comps <-  c("BC", "DU", "OM", "SS", "SU", "PM25")

df_comps_lst <- vector("list", length = length(comps))
df_yr_lst <- vector("list", length = length(comps))
for (i in seq_along(comps)) {
  df_comps_lst[[i]] <-  df_all[c("Lon", "Lat", comps[i], "Year", "ID")]
  
  df_yr_lst[[i]] <- df_comps_lst[[i]] %>% group_by(Year) 
  
}

df_yr_lst[[6]]


#------------------------#
# Calculate Temporal Mean
#------------------------#

emp_tp_mu_BC <- df_yr_lst[[1]] %>% summarise(emp_tp_mu_BC = mean(BC))
head(emp_tp_mu_BC) # 27384 *3

emp_tp_mu_DU <- df_yr_lst[[2]] %>% summarise(emp_tp_mu_DU = mean(DU))
head(emp_tp_mu_DU)

emp_tp_mu_OM <- df_yr_lst[[3]] %>% summarise(emp_tp_mu_OM = mean(OM))
head(emp_tp_mu_OM)

emp_tp_mu_SS <- df_yr_lst[[4]] %>% summarise(emp_tp_mu_SS = mean(SS))
head(emp_tp_mu_SS)

emp_tp_mu_SU <- df_yr_lst[[5]] %>% summarise(emp_tp_mu_SU = mean(SU))
head(emp_tp_mu_SU)

emp_tp_mu_PM25 <- df_yr_lst[[6]] %>% summarise(emp_tp_mu_PM25 = mean(PM25))
head(emp_tp_mu_PM25)


#------------#
# Add Labels
#------------#

Lab_bc <- rep("BC", nrow(emp_mu_BC))
emp_mu_BC$Label <- rep("BC", nrow(emp_mu_BC))


Labs <- list("BC", "DU", "OM", "SS", "SU", "PM25")
emp_tp_mu_lst <- list(emp_tp_mu_BC, emp_tp_mu_DU, emp_tp_mu_OM,
                   emp_tp_mu_SS, emp_tp_mu_SU, emp_tp_mu_PM25)

for (i in seq_along(Labs)) {
  emp_tp_mu_lst[[i]]$Label <- rep(Labs[[i]], nrow(emp_tp_mu_lst[[i]]))
}

head(emp_tp_mu_lst[[1]])



#==================#
# Add Label columns 
#==================#
# prepare for ggplot facet_wrap

nrow(emp_tp_mu_BC) # [1] 5

Lab_bc <- rep("BC", nrow(emp_tp_mu_BC))
emp_tp_mu_BC$Label <- rep("BC", nrow(emp_tp_mu_BC))


Labs <- list("BC", "DU", "OM", "SS", "SU", "PM25")
emp_tp_mu_lst <- list(emp_tp_mu_BC, emp_tp_mu_DU, emp_tp_mu_OM,
                   emp_tp_mu_SS, emp_tp_mu_SU, emp_tp_mu_PM25)

for (i in seq_along(Labs)) {
  emp_tp_mu_lst[[i]]$Label <- rep(Labs[[i]], nrow(emp_tp_mu_lst[[i]]))
}

head(emp_tp_mu_lst[[1]])
head(emp_tp_mu_lst[[2]])


#------------#
# Rename df
#-----------#

names(emp_tp_mu_lst[[1]]) <- c("Year", "emp_tp_mu", "Label") 

name_lst <- list(c("Year", "emp_tp_mu", "Label"))

for(i in seq(1:6)) {
  names(emp_tp_mu_lst[[i]]) <- name_lst[[1]]
}
head(emp_tp_mu_lst[[2]])


#-----------------------#
# Row bind different df 
#-----------------------#

# to make long format

emp_tp_mu_df_lng <- bind_rows(emp_tp_mu_lst[[1]], emp_tp_mu_lst[[2]],
                           emp_tp_mu_lst[[3]], emp_tp_mu_lst[[4]],
                           emp_tp_mu_lst[[5]], emp_tp_mu_lst[[6]])


head(emp_tp_mu_df_lng)
tail(emp_tp_mu_df_lng)



#======#
# Plot
#======#

head(df_all)

#---------------------------------------------#
# emp_tp_mu_BC together with all stations (ref)
#---------------------------------------------#

plt_emp_tp_mu_BC <- df_all %>% select(Lon, Lat, BC, Year, ID) %>%
  ggplot() + 
  geom_line(aes(x = Year, y = BC, group = ID),
            colour = "blue", alpha = 0.04) +
  geom_line(data = emp_tp_mu_lst[[1]],
            aes(x = Year, y = emp_tp_mu)) + 
  xlab("Year") + 
  ylab("Empirical Temporal Mean") + 
  theme_bw()

print(plt_emp_tp_mu_BC)


#-------------------#
# emp_tp_mu_6_comps
#-------------------#


P_tp <- list()
plot_tp_df <- list()
for (i in unique(emp_tp_mu_df_lng$Label)) {
  plot_tp_df[[i]] <- filter(emp_tp_mu_df_lng, Label == i)
  
  P_tp[[i]] <- ggplot(plot_tp_df[[i]]) + 
    geom_line(aes(x = Year, y = emp_tp_mu),
               colour = "LightSkyBlue") +
    xlab("Year") + 
    ylab("Temporal Mean") + 
    facet_wrap(~Label) +
    theme_light()
} 

do.call(grid.arrange, P_tp)






