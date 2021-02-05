#*****************************#
# Sp_T_Exploratory:
# Empirical Covariance Matrix
#*****************************#

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


#------------------#
# Detrend the data
#------------------#

lm_BC <- lm(BC ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
           Year + I(Year) ^ 2, data = Comps_lst[[1]])

lm_DU <- lm(DU ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
              Year + I(Year) ^ 2, data = Comps_lst[[2]])

lm_OM <- lm(OM ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
              Year + I(Year) ^ 2, data = Comps_lst[[3]])

lm_SS <- lm(SS ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
              Year + I(Year) ^ 2, data = Comps_lst[[4]])

lm_SU <- lm(SU ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
              Year + I(Year) ^ 2, data = Comps_lst[[5]])

lm_PM25 <- lm(PM25 ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
                Year + I(Year) ^ 2, data = Comps_lst[[6]])


lm_lst <- list(lm_BC, lm_DU, lm_OM, lm_SS, lm_SU, lm_PM25)


#-------------------#
# Extract Residuals
#-------------------#

res_lst <- list()
for (i in seq_along(lm_lst)) {
  res_lst[[i]] <- lm_lst[[i]]$residuals
  Comps_lst[[i]]$Residuals <- res_lst[[i]]
}

head(Comps_lst[[1]], 2)
#         Lon   Lat          BC Year ID Residuals
# 4503 -42.75 83.25 0.007900851 2016  1 0.3865688
# 4504 -42.00 83.25 0.006894099 2016  2 0.3848602


#-------------------------------------#
# Range of residuals (Component-wise)
#-------------------------------------#

range(res_lst[[1]])

Range_lst <- list()
Med_lst <- list()
for (i in seq_along(res_lst)) {
  Range_lst[[i]] <- range(res_lst[[i]])
  Med_lst[[i]] <- median(res_lst[[i]])
}


lapply(Range_lst, print)
# [1] -1.153176 53.767906
# [1] -18.95653 953.58374
# [1] -16.6651 202.1438
# [1] -1.449135  9.988610
# [1] -2.093822 94.155612
# [1] -35.81649 938.14513

lapply(Med_lst, print)
# [1] -0.1853067
# [1] -5.706429
# [1] -2.411226
# [1] -0.2750538
# [1] -0.328552
# [1] -5.849559



#===========================#
# Get the spatial locations
#===========================#

str(Comps_lst[[1]])
sp_loca_df <- Comps_lst[[1]] %>% 
  filter(Year == 2016) %>%
  select(Lon, Lat) %>% 
  arrange(Lon, Lat)

head(sp_loca_df) # each location 5 repetitions (yrs)
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



lm_lst <- list()
for (i in seq_along(Comps_lst)) {
  lm_lst[[i]] <- lm(Comps_lst[[i]][3] ~ Lon + I(Lon ^ 2) + Lat + I(Lat ^ 2) + 
                      Year + I(Year) ^ 2, data = Comps_lst[[i]])
}
#~~~~~~~~~~~~~~~~~~



#============================#
# Empirical Covariance Matrix
#============================#

head(Comps_lst[[1]])

#----------------------#
# Put df into space-wide
#----------------------#

## single ref
X_bc <- select(Comps_lst[[1]], Lon, Lat, Year, Residuals) %>%
  spread(key = Year, value = Residuals) %>%
  select(-Lon, -Lat) %>%
  t()

str(X_bc) # num [1:5, 1:27384]: 
# each row: Year
# each col: Residuals at each pair of (lon, lat)


## all components
X_lst <- list()
for (i in seq_along(Comps_lst)) {
  X_lst[[i]] <- select(Comps_lst[[i]], Lon, Lat, Year, Residuals) %>%
    spread(key = Year, value = Residuals) %>%
    select(-Lon, -Lat) %>%
    t()
}

str(X_lst)
# List of 6
# $ : num [1:5, 1:27384]



#============================#
# Covariance (Component-wise)
#============================#


#------------------#
# Lag-0 covariance
#------------------#

library(parallel)
Lag_0_cov_lst<- mclapply(X_lst, cov, use = "complete.obs", mc.cores = 6)

Lag_0_cov_lst <- list()
for(i in seq_along(X_lst)) {
  Lag_0_cov_lst[[i]] <- cov(X_lst[[i]], use = "complete.obs")
}



cov(X_lst, use = "complete.obs")
