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


lm_PM25_try <- lm(PM25 ~  Year , data = Comps_lst[[6]])


lm_lst <- list(lm_BC, lm_DU, lm_OM, lm_SS, lm_SU, lm_PM25)



all(lm_BC$residuals - residuals(lm_BC) == 0) # [1] TRUE


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


Comps_lst[[6]]$Residuals_try <- lm_PM25_try$residuals



#------------------------------------------#
# Range/Median of residuals (Component-wise)
#------------------------------------------#

range(res_lst[[1]])

Range_lst <- list()
Med_lst <- list()
Q_lst <- list()
for (i in seq_along(res_lst)) {
  Range_lst[[i]] <- range(res_lst[[i]])
  Med_lst[[i]] <- median(res_lst[[i]])
  Q_lst[[i]] <- quantile(res_lst[[i]])
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

print(Q_lst)

quantile(res_pm25)
#         0%        25%        50% 
# -20.522429 -15.726288 -10.378623 
# 75%       100% 
#  2.798213 952.304762



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

#~~~~~~~~~~~#
# single ref
#~~~~~~~~~~~#

X_bc <- select(Comps_lst[[1]], Lon, Lat, Year, Residuals) %>%
  spread(key = Year, value = Residuals) %>%
  select(-Lon, -Lat) %>%
  t()

str(X_bc) # num [1:5, 1:27384]: 
# each row: Year
# each col: Residuals at each pair of (lon, lat)


#~~~~~~~~~~#
# PM25 try
#~~~~~~~~~~#
head(Comps_lst[[6]], 2)

X_pm25 <- select(Comps_lst[[6]], Lon, Lat, Year, Residuals_try) %>%
  spread(key = Year, value = Residuals_try) %>%
  select(-Lon, -Lat) %>%
  t()

str(X_pm25) # num [1:5, 1:27384]


#~~~~~~~~~~~~~~~~~#
# filter only 2016 
#~~~~~~~~~~~~~~~~~#

X_bc_16 <- select(Comps_lst[[1]], Lon, Lat, Year, Residuals) %>%
  filter(Year == 2016) %>%
  spread(key = Year, value = Residuals) %>%
  select(-Lon, -Lat) %>%
  t()


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


## 2016
X_lst_16 <- list()
for (i in seq_along(Comps_lst)) {
  X_lst_16[[i]] <- select(Comps_lst[[i]], Lon, Lat, Year, Residuals) %>%
    filter(Year == 2016) %>%
    spread(key = Year, value = Residuals) %>%
    select(-Lon, -Lat) %>%
    t()
}
str(X_lst_16) # List of 6
# $ : num [1, 1:27384]



#============================#
# Covariance (Component-wise)
#============================#

#load("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_Sp_T_Exp_Covariance.RData")


#------------------#
# Lag-0 covariance
#------------------#

Lag_0_cov_BC <- cov(X_lst[[1]], use = "complete.obs")
Lag_0_cov_DU <- cov(X_lst[[2]], use = "complete.obs")
Lag_0_cov_OM <- cov(X_lst[[3]], use = "complete.obs")
Lag_0_cov_SS <- cov(X_lst[[4]], use = "complete.obs")
Lag_0_cov_SU <- cov(X_lst[[5]], use = "complete.obs")
Lag_0_cov_PM25 <- cov(X_lst[[6]], use = "complete.obs")

rm(Lag_0_cov_BC)



#~~~~~~~~#
# PM25 try
#~~~~~~~~#

Lag_0_cov_PM25_try <- cov(X_pm25, use = "complete.obs")

quantile(Lag_0_cov_PM25_try)


#~~~~~~~~~~~~~~~~~~~~#
# Understand results
#~~~~~~~~~~~~~~~~~~~~#

## BC
range(Lag_0_cov_BC) # [1] -60.96732 565.41161
#(c(-60.96732, 565.41161)) ^ (1/3)
#-(60) ^ (1/3)

median(Lag_0_cov_BC) # 4.841124e-05
str(Lag_0_cov_BC) # num [1:27384, 1:27384]

a <- which(Lag_0_cov_BC > 0)
length(a) # 391061578

a1 <- which(Lag_0_cov_BC > 1)
length(a1) # 525748

a2 <- which(Lag_0_cov_BC > 2)
length(a2) # 183123

a3 <- which(Lag_0_cov_BC > 3) # 95767
a50 <- which(Lag_0_cov_BC > 50) # 319


range(Lag_0_cov_DU) # [1] -2576.171 15125.574

range(Lag_0_cov_PM25) # [1] -6077.674 15113.980


#-----------------------------#
# Lag-0 Cross-Covariance (2016)
#------------------------------#

cov

## BC vs PM25
Lag_0_cov_BP_16 <- cov(X_lst_16[[1]], X_lst_16[[6]], use = "complete.obs")

range(Lag_0_cov_BP_16)
str(Lag_0_cov_BP_16)

#------------------#
# Lag-1 covariance
#------------------#

nrow(X_lst[[1]]) # 5

Lag_1_cov_BC <- cov(X_lst[[1]][-1, ], X_lst[[1]][-5, ], use = "complete.obs")
Lag_1_cov_DU <- cov(X_lst[[2]][-1, ], X_lst[[2]][-5, ], use = "complete.obs")
Lag_1_cov_OM <- cov(X_lst[[3]][-1, ], X_lst[[3]][-5, ], use = "complete.obs")



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


#-----------------#
# Plot Covariance
#-----------------#

# Need a Math.cbrt function to enable the calculation of cbrt of negative number
Math.cbrt <- function(x) {
  sign(x) * abs(x) ^ (1/3)
}


EMP_cov_lat_plt <- function(C, sp_loca_df) {
  require(fields)
  
  for (i in seq_along(unique(sp_loca_df$lon_strip))) {
    sp_strip <-  filter(sp_loca_df,  lon_strip  == i) %>% 
      arrange (Lat)
    
    idx <- sp_strip$n
    jitter <- seq(0, 1e-4, length = length(idx))
    
    image.plot(x = sp_strip$Lat + jitter, 
               y = sp_strip$Lat + jitter,
               z = Math.cbrt(C[idx, idx]),
               zlim = c(-1, 25),
               xlab = "Latitude", ylab = "Latitude",
               col = tim.colors(10), cex = 200)
  }
}


pt <- plot_cov_strips(C = Lag_0_cov_PM25, spat_df = sp_loca_df)

str(Lag_0_cov_BC)


pt <- EMP_cov_lat_plt(Lag_0_cov_PM25, sp_loca_df = sp_loca_df)


EMP_cov_lat_plt(Lag_0_cov_PM25_try, sp_loca_df)



########
image.plot(x = order(sp_loca_df$Lat), y = order(sp_loca_df$Lon),
           z = (Lag_0_cov_BC) ^ (1/3),
           xlab = "Latitude", ylab = "Longitutde",
           col = tim.colors(10))

head(sp_loca_df$Lat)

range(sp_loca_df$Lat)
range(sp_loca_df$Lon)
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






