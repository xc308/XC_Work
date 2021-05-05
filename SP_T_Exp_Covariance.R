#*****************************#
# Sp_T_Exploratory:
# Empirical Covariance Matrix
#*****************************#
rm(list = ls())

load("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_Sp_T_Exploratory.RData")



#=====================#
# Separate comoponents
#=====================#

head(df_all)
tail(df_all)

range(df_all[, "PM25"]) # [1]   0.0974795 972.9597204


#------------------------------------#
# Separate df_all into component wise
#------------------------------------#

#head(df_all[c("Lon", "Lat", "BC", "Year", "ID")])

comps <- c("BC", "DU", "OM", "SS", "SU", "PM25")
Comps_lst <- list()
for (i in seq_along(comps)) {
  Comps_lst[[i]] <- df_all[c("Lon", "Lat", comps[i], "Year", "ID")]
}

head(Comps_lst[[1]])
str(Comps_lst)



head(Comps_lst[[1]][3])
#             BC
#4503 0.007900851
#4504 0.006894099



#===========================#
# Get the spatial locations
#===========================#
install.packages("dplyr")
install.packages("tidyr")
library(dplyr)
library(tidyr)


str(Comps_lst[[1]])

sp_loca_df <- Comps_lst[[1]] %>% 
  filter(Year == 2016) %>%
  select(Lon, Lat) %>% 
  arrange(Lon, Lat)

head(sp_loca_df) #
str(sp_loca_df)  # 'data.frame':	27384 obs. of  2 variables:
tail(sp_loca_df)


#----------------#
# Recall emp_mu
#----------------#

nrow(emp_mu_bc) # 27384


#~~~~~~~~~~~#
# Question
#~~~~~~~~~~~#
## BC
lm <- lm(Comps_lst[[1]][3] ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
     Year + I(Year) ^ 2, data = Comps_lst[[1]])

str(Comps_lst)


lm_lst <- list()
for (i in seq_along(Comps_lst)) {
  lm_lst[[i]] <- lm(names(Comps_lst[[i]][3]) ~ 1, 
                    data = Comps_lst[[i]])
}
#~~~~~~~~~~~~~~~~~~


#------------------#
# Detrend the data
#------------------#

lm_BC <- lm(BC ~ Lon + Lat + I(Lat ^ 2) + 
              Year + I(Year ^ 2) , data = Comps_lst[[1]])

lm_DU <- lm(DU ~ Lon  + Lat + I(Lat ^ 2) + 
              Year + I(Year ^ 2), data = Comps_lst[[2]])

lm_OM <- lm(OM ~ Lon + Lat + I(Lat ^ 2) + 
              Year + I(Year ^ 2), data = Comps_lst[[3]])

lm_SS <- lm(SS ~ Lon + Lat + I(Lat ^ 2) + 
              Year + I(Year ^ 2), data = Comps_lst[[4]])

lm_SU <- lm(SU ~ Lon + Lat + I(Lat ^ 2) + 
              Year + I(Year ^ 2), data = Comps_lst[[5]])

lm_PM25 <- lm(PM25 ~ Lon + Lat + I(Lat ^ 2) + 
                Year + I(Year ^ 2), data = Comps_lst[[6]])


lm_lst <- list(lm_BC, lm_DU, lm_OM, lm_SS, lm_SU, lm_PM25)


#-------------------#
# Extract Residuals
#-------------------#

res_lst <- list()
for (i in seq_along(lm_lst)) {
  res_lst[[i]] <- lm_lst[[i]]$residuals
  Comps_lst[[i]]$Residuals <- res_lst[[i]]
}

head(Comps_lst[[1]])

range(Comps_lst[[1]]$Residuals)
hist(Comps_lst[[1]]$Residuals)



#============================#
# Empirical Covariance Matrix
#============================#

head(Comps_lst[[1]])

#----------------------#
# Put df into space-wide
#----------------------#

#~~~~~~~~~~~#
# single ref
#~~~~~~~~~~~#

X_bc <- select(Comps_lst[[1]], Lon, Lat, Year, Residuals) %>%
  spread(key = Year, value = Residuals) %>%
  select(-Lon, -Lat) %>%
  t()

str(X_bc) # num [1:5, 1:27384]: 



#------------------------------------------
X_BC_2 <- select(Comps_lst[[1]], Lon, Lat, Year, BC) %>%
  arrange(Lon, Lat) %>%
  spread(key = Year, value = BC) %>%
  select(-Lon, -Lat) %>%
  t()


str(X_BC_2)
all(X_BC_2[1, ] == X_bc[1, ]) # TRUE

head(X_BC_2[1, ])
head(X_bc[1, ])

#-------------------------------------------


#~~~~~~~~~~~~~~~#
#all components
#~~~~~~~~~~~~~~~#

head(Comps_lst[[1]])


X_lst <- list()
for (i in seq_along(Comps_lst)) {
  X_lst[[i]] <- select(Comps_lst[[i]], Lon, Lat, Year, names(Comps_lst[[i]][6])) %>%
    spread(key = Year, value = names(Comps_lst[[i]][6])) %>%
    select(-Lon, -Lat) %>%
    t()
}


str(X_lst)
# List of 6
# $ : num [1:5, 1:27384]

range(X_lst[[1]]) # [1]   0.0974795 972.9597204



#============================#
# Covariance (Component-wise)
#============================#

load("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_Sp_T_Exp_Covariance.RData")


#------------------#
# Lag-0 covariance
#------------------#

## To examine for component-wise non-stationarity

#Lag_0_cov_BC <- cov(X_lst[[1]], use = "complete.obs")
#Lag_0_cov_DU <- cov(X_lst[[2]], use = "complete.obs")
#Lag_0_cov_OM <- cov(X_lst[[3]], use = "complete.obs")
#Lag_0_cov_SS <- cov(X_lst[[4]], use = "complete.obs")
#Lag_0_cov_SU <- cov(X_lst[[5]], use = "complete.obs")
#Lag_0_cov_PM25 <- cov(X_lst[[6]], use = "complete.obs")
#str(Lag_0_cov_PM25)
# num [1:27384, 1:27384] 0.04884 0.0439 0.035 0.00913 0.0


#------------------#
# Lag-0 correlation
#------------------#

Lag_0_corr_BC <- cor(X_lst[[1]], use = "complete.obs")
hist(Lag_0_corr_BC)



#-----------------------------#
# Lag-1 covariance/correlation
#-----------------------------#

# To examine if space and time are separable

nrow(X_lst[[1]]) # 5

#Lag_1_cov_BC <- cov(X_lst[[1]][-1, ], X_lst[[1]][-5, ], use = "complete.obs")
Lag_1_corr_BC <- cor(X_lst[[1]][-1, ], X_lst[[1]][-5, ], use = "complete.obs")

hist(Lag_1_corr_BC)



#Lag_1_cov_DU <- cov(X_lst[[2]][-1, ], X_lst[[2]][-5, ], use = "complete.obs")
#Lag_1_cov_OM <- cov(X_lst[[3]][-1, ], X_lst[[3]][-5, ], use = "complete.obs")

Lag_1_cov_PM25 <- cov(X_lst[[6]][-1, ], X_lst[[6]][-5, ], use = "complete.obs")



#==============================#
# cross-covariance/correlation
#==============================#

# To examine the stationarity of cross-covariance btw components

str(X_lst)

#~~~~~~~~~~~#
# PM25 & DU
#~~~~~~~~~~~#

cross_cov_PM25_DU <- cov(X_lst[[6]], X_lst[[2]], use = "complete.obs")
#cross_corr_PM25_BC <- cov2cor(cross_cov_PM25_BC)
range(cross_cov_PM25_DU)

all(diag(cross_cov_PM25_DU) > 0) # [1] FALSE

range(diag(cross_cov_PM25_DU)) # [1]   -14.19178 15118.69358
hist(cross_cov_PM25_DU)
quantile(cross_cov_PM25_DU)


## on log scale
cross_cov_PM25_DU_log <- cov(log(X_lst[[6]]), log(X_lst[[2]]), use = "complete.obs")
range(cross_cov_PM25_DU_log)  # [1] -1.814203  1.782933
quantile(cross_cov_PM25_DU_log)
range(diag(cross_cov_PM25_DU_log))

hist(cross_cov_PM25_DU_log)


#~~~~~~~~~~#
# PM25 & OM
#~~~~~~~~~~#

cross_cov_PM25_OM_log <- cov(log(X_lst[[6]]), log(X_lst[[3]]), use = "complete.obs")
range(cross_cov_PM25_OM_log) # [1] -2.132976  3.475127
quantile(cross_cov_PM25_OM_log)
hist(cross_cov_PM25_OM_log)


#~~~~~~~~~~#
# PM25 & BC
#~~~~~~~~~~#

cross_cov_PM25_BC_log <- cov(log(X_lst[[6]]), log(X_lst[[1]]), use = "complete.obs")
range(cross_cov_PM25_BC_log) # [1] -2.363203  3.830183
quantile(cross_cov_PM25_BC_log)
hist(cross_cov_PM25_BC_log)


#~~~~~~~~~~#
# PM25 & SU
#~~~~~~~~~~#

cross_cov_PM25_SU_log <- cov(log(X_lst[[6]]), log(X_lst[[5]]), use = "complete.obs")
range(cross_cov_PM25_SU_log) # [1] -0.5525769  0.7456199
hist(cross_cov_PM25_SU_log)


#~~~~~~~~~~#
# PM25 & SS
#~~~~~~~~~~#

cross_cov_PM25_SS_log <- cov(log(X_lst[[6]]), log(X_lst[[4]]), use = "complete.obs")
range(cross_cov_PM25_SS_log) # [1] -0.6979235  1.3595946
quantile(cross_cov_PM25_SS_log)



#======#
# Plot
#======#

#----------------------------#
# Divid the Lon into 4 strips
#----------------------------#

sp_loca_df$n <- 1:nrow(sp_loca_df) # index each loaction
lim_lon <- range(sp_loca_df$Lon)  # -179.25  180.00
lon_strips <- seq(lim_lon[1], lim_lon[2], length.out = 5)

sp_loca_df$lon_strip <- cut(sp_loca_df$Lon, lon_strips, 
                            labels = FALSE, include.lowest = TRUE)

head(sp_loca_df)
tail(sp_loca_df)
nrow(sp_loca_df) # 27384


#-----------------------------#
# Plot Covariance/ Correlation
#-----------------------------#
install.packages("fields")
install.packages("sp")
library(fields)
library(sp)


#~~~~~~~~~~~~~#
# Choose color
#~~~~~~~~~~~~~#

install.packages("RColorBrewer")
library(RColorBrewer)

cols_Ac <- brewer.pal(8, name = "Accent")
brewer_AC.div <- colorRampPalette(cols_Ac, interpolate = "spline")


#~~~~~~~~~~~~~~#
# Plot function
#~~~~~~~~~~~~~~#

EMP_cor_lat_plt <- function(C, sp_loca_df) {
  require(fields)
  
  for (i in seq_along(unique(sp_loca_df$lon_strip))) {
    sp_strip <-  filter(sp_loca_df,  lon_strip  == i) %>% 
      arrange (Lat)
    
    idx <- sp_strip$n
    jitter <- seq(0, 1e-4, length = length(idx))
    
    image.plot(x = sp_strip$Lat + jitter, 
               y = sp_strip$Lat + jitter,
               z = C[idx, idx],
               zlim = c(-1, 1),
               #zlim = c(-0.1, 0.1),
               xlab = "Latitude", ylab = "Latitude",
               col = tim.colors(100), 
               #col = brewer_AC.div(500),
               cex = 200)
  }
}



#--------------------------#
# Univairate Lag_o_corr plots
#--------------------------#

## Lag_0_corr_BC
EMP_cor_lat_plt(Lag_0_corr_BC, sp_loca_df)




## Lag_0_corr_PM25 for single component stationarity
EMP_cor_lat_plt(Lag_0_corr_PM25, sp_loca_df)

## Lag_0_cov_PM25_log for single component stationarity
#EMP_cor_lat_plt(Lag_0_cov_PM25_log, sp_loca_df)





#----------------#
# Lag_1_corr_ uni
#----------------#

EMP_cor_lat_plt(Lag_1_corr_BC, sp_loca_df)



## cross-covariance of PM25 and DU log
EMP_cor_lat_plt(cross_cov_PM25_DU_log, sp_loca_df)

## cross-covariance of PM25 and OM log
EMP_cor_lat_plt(cross_cov_PM25_OM_log, sp_loca_df)


## cross-covariance of PM25 and BC log
EMP_cor_lat_plt(cross_cov_PM25_BC_log, sp_loca_df)


## cross-covariance of PM25 and SU log
EMP_cor_lat_plt(cross_cov_PM25_SU_log, sp_loca_df)


## cross-covariance of PM25 and SS log
EMP_cor_lat_plt(cross_cov_PM25_SS_log, sp_loca_df)



#-------------------#
# Examine one strip
#------------------#

sp_strip <- filter(sp_loca_df,  lon_strip  == 1) %>% 
      arrange (Lat) # 3793
    
    idx <- sp_strip$n # 42, 34, 3610, 3611, 3671,....
    jitter <- seq(0, 1e-4, length = length(idx))
    
    image.plot(x = sp_strip$Lat + jitter, # ensure 3793 different Lat locations
               y = sp_strip$Lat + jitter, 
               z = C[idx, idx], #3793*3793
               zlim = c(-1, 1),
               #zlim = c(-0.1, 0.1),
               xlab = "Latitude", ylab = "Latitude",
               #col = tim.colors(500), 
               col = brewer_AC.div(500),
               cex = 800)


head(sp_strip$Lat, 50)
Lag_0_cov_PM25[42, 42] #  0.1065634

filter(sp_strip, n == 42)
#          Lon    Lat  n lon_strip
# 66734 -169.5 -14.25 42     

idx <- c(42, 34, 3610)

Lag_0_cov_BC[idx, idx]

Lag_0_cov_BC[42, 42]
Lag_0_cov_BC[34, 34]
1



#-----------------------#
# Too large, reach limit
#-----------------------#

library(parallel)
Lag_0_cov_lst<- mclapply(X_lst, cov, use = "complete.obs", mc.cores = 6)

Lag_0_cov_lst <- list()
for(i in seq_along(X_lst)) {
  Lag_0_cov_lst[[i]] <- cov(X_lst[[i]], use = "complete.obs")
}

#-------#
# Error
#-------#
# Error: vector memory exhausted (limit reached?) 
# ref: https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos






