#*****************************#
# Sp_T_Exploratory:
# Empirical Covariance Matrix
#*****************************#
rm(list = ls())

load("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_Sp_T_Exploratory.RData")



#==================#
# Remove all trends
#==================#

head(df_all)
tail(df_all)


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

str(Comps_lst[[1]])
sp_loca_df <- Comps_lst[[1]] %>% 
  filter(Year == 2016) %>%
  select(Lon, Lat) %>% 
  arrange(Lon, Lat)

head(sp_loca_df) #
str(sp_loca_df)  # 'data.frame':	27384 obs. of  2 variables:



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

X_bc <- select(Comps_lst[[1]], Lon, Lat, Year, BC) %>%
  spread(key = Year, value = BC) %>%
  select(-Lon, -Lat) %>%
  t()

str(X_bc) # num [1:5, 1:27384]: 


#~~~~~~~~~~~~~~~#
#all components
#~~~~~~~~~~~~~~~#

head(Comps_lst[[i]])

X_lst <- list()
for (i in seq_along(Comps_lst)) {
  X_lst[[i]] <- dplyr::select(Comps_lst[[i]], Lon, Lat, Year, names(Comps_lst[[i]][3])) %>%
    spread(key = Year, value = names(Comps_lst[[i]][3])) %>%
    dplyr::select(-Lon, -Lat) %>%
    t()
}




str(X_lst)
# List of 6
# $ : num [1:5, 1:27384]



#============================#
# Covariance (Component-wise)
#============================#

#load("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_Sp_T_Exp_Covariance.RData")


#------------------#
# Lag-0 covariance
#------------------#

Lag_0_cov_BC <- cov(X_bc, use = "complete.obs")

Lag_0_cov_DU <- cov(X_lst[[2]], use = "complete.obs")
Lag_0_cov_OM <- cov(X_lst[[3]], use = "complete.obs")
Lag_0_cov_SS <- cov(X_lst[[4]], use = "complete.obs")
Lag_0_cov_SU <- cov(X_lst[[5]], use = "complete.obs")
Lag_0_cov_PM25 <- cov(X_lst[[6]], use = "complete.obs")

Math.cbrt(quantile(Lag_0_cov_BC))
quantile(Lag_0_cov_PM25)


cor_PM25 <- cov2cor(Lag_0_cov_PM25)
quantile(cor_PM25)
hist(cor_PM25)


#------------------#
# Lag-1 covariance
#------------------#

nrow(X_lst[[1]]) # 5

#Lag_1_cov_BC <- cov(X_lst[[1]][-1, ], X_lst[[1]][-5, ], use = "complete.obs")
#Lag_1_cov_DU <- cov(X_lst[[2]][-1, ], X_lst[[2]][-5, ], use = "complete.obs")
#Lag_1_cov_OM <- cov(X_lst[[3]][-1, ], X_lst[[3]][-5, ], use = "complete.obs")

Lag_1_cov_PM25 <- cov(X_lst[[6]][-1, ], X_lst[[6]][-5, ], use = "complete.obs")



#======#
# Plot
#======#

#----------------------------#
# Divid the Lon into 4 strips
#----------------------------#

sp_loca_df$n <- 1:nrow(sp_loca_df) # index each loaction
lim_lon <- range(sp_loca_df$Lon)
lon_strips <- seq(lim_lon[1], lim_lon[2], length.out = 5)

sp_loca_df$lon_strip <- cut(sp_loca_df$Lon, lon_strips, labels = FALSE, include.lowest = TRUE)

head(sp_loca_df)
tail(sp_loca_df)
nrow(sp_loca_df) # 27384


#-----------------#
# Plot Covariance
#-----------------#

# Need a Math.cbrt function to enable the calculation of cbrt of negative number
Math.cbrt <- function(x) {
  sign(x) * abs(x) ^ (1/3)
}


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
               #zlim = c(-1, 1),
               xlab = "Latitude", ylab = "Latitude",
               col = tim.colors(500), cex = 800)
  }
}



EMP_cor_lat_plt(cor_PM25, sp_loca_df)



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






