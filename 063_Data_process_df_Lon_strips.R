#===========
# 9 May 2024
#===========

# Aim:
  # process the data further for parallelisation
  
# Method:
  # 1. divide the whole data into 8 strips by Lon_1, ..., Lon_8
  # 2. within each Lon strip, Lat decrease
  # 3. and Lat for Lon_3, Lon_5 will be further cut by Lat > 0, Lat < 0
    # Lon3_pos, Lon3_neg, Lon5_pos, Lon5_neg, Lon7_pos, Lon7_neg, Lon8_pos, Lon8_neg

# Reason：
  # 1. the HPC job longest running time is 10 days
  # 2. 200 locations with 5 CAMS variates, parallelize over 51 CPU cores require
      # 10 minutes 29 seconds
  # so need to limit the data chunk for HPC to be able to finish within 10 days


#=============
# Pre-settings
#=============
install.packages("dplyr")
library(dplyr)

load("~/Library/CloudStorage/OneDrive-UniversityofExeter/XC_PhD/Data/Processed/XC_WORK/RData/df_Cmpt_Res_log_16.rda")


head(df_Cmpt_Res_log_16)
# Lon Lat Year BC_Residuals_log DU_Residuals_log OM_Residuals_log
# SS_Residuals_log SU_Residuals_log PM25_Residuals_log


# save to a new copy
df_Res_log_16 <- df_Cmpt_Res_log_16

range(df_Res_log_16$Lon)
# [1] -179.25  180.00

range(df_Res_log_16$Lat)
# [1] -55.50  83.25

# Define the latitude range and step
start_lon <- -179.25
end_lon <- 180

# The number of division
num_div <- 8 

# division size
div_size_lon = (end_lon - start_lon) / num_div
# [1] 44.90625


breaks_lon = seq(start_lon, end_lon, length.out = num_div + 1)
length(breaks_lon) # [1] 9
breaks_lon
# [1] -179.25000 -134.34375 -89.43750  -44.53125
# [5]    0.37500   45.28125 90.18750  135.09375  180.00000
# 

diff(breaks_lon)
# [1] 44.90625 44.90625 44.90625
# [4] 44.90625 44.90625 44.90625
# [7] 44.90625 44.90625


df_Res_log_16 <- df_Res_log_16 %>%
  mutate(Lon_div_indx = cut(Lon, breaks = breaks_lon, labels = F))


#--------------------------------
# check NA for the boundary cut
#--------------------------------

which(is.na(df_Res_log_16$Lon_div_indx)) 

# [1] 1772 2782 3114 3474

df_Res_log_16[1770:1779,]
df_Res_log_16[2780:2784, ]
df_Res_log_16[3112:3116, ]
df_Res_log_16[3472:3475,]
#        Lon   Lat   Year
# 1772 -179.25 71.25 2016
# 2782 -179.25 68.25 2016
# 3114 -179.25 67.50 2016
# 3474 -179.25 66.75 2016

## all corresponds to Lon: -179.25
# so shall categorize to Lon_div_indx = 1

df_Res_log_16$Lon_div_indx[1772] <- 1
df_Res_log_16$Lon_div_indx[2782] <- 1
df_Res_log_16$Lon_div_indx[3114] <- 1
df_Res_log_16$Lon_div_indx[3474] <- 1

nrow(df_Res_log_16) # [1] 27384

df_div_Lon_lst <- split(df_Res_log_16, df_Res_log_16$Lon_div_indx)
str(df_div_Lon_lst)

#List of 8
# $ 1:'data.frame':	630 obs. of  10 variables:
# $ 2:'data.frame':	3163 obs. of  10 variables:
# $ 3:'data.frame':	4648 obs. of  10 variables:
# $ 4:'data.frame':	1933 obs. of  10 variables:
# $ 5:'data.frame':	5910 obs. of  10 variables:
# $ 6:'data.frame':	4245 obs. of  10 variables:
# $ 7:'data.frame':	4978 obs. of  10 variables:
# $ 8:'data.frame':	1877 obs. of  10 variables:


#=====================
# GPU running strategy
#=====================

# 1, 2 list corresponds to Lon_strip 1: 3793 locations;
# 3, 4 list corresponds to Lon_strip 2: 6581 locations;
# 5, corresponds to half of Lon_strip 3: 5910 locations;
# 6 corresponds to the other half of Lon_strip 3: 4245 locations;
# 7, 8 list corresponds to Lon_strip 4: 6855 locations.


#------------------------------------
# bind into Lon_strip_1, 2, 31, 32, 4
#------------------------------------

df_Lon_strip_1 <- do.call(rbind, list(df_div_Lon_lst[[1]], df_div_Lon_lst[[2]]))
str(df_Lon_strip_1)
#'data.frame':	3793 obs. of  10 variables:
unique(df_Lon_strip_1$Lon_div_indx)
# [1] 1 2


df_Lon_strip_2 <- do.call(rbind, list(df_div_Lon_lst[[3]], df_div_Lon_lst[[4]]))
str(df_Lon_strip_2)
# 'data.frame':	6581 obs. of  10 variables:
unique(df_Lon_strip_2$Lon_div_indx)
# [1] 3 4


df_Lon_strip_31 <- df_div_Lon_lst[[5]]
str(df_Lon_strip_31)
unique(df_Lon_strip_31$Lon_div_indx)
# [1] 5


df_Lon_strip_32 <- df_div_Lon_lst[[6]]
str(df_Lon_strip_32)
# 'data.frame':	4245 obs. of  10 variables
unique(df_Lon_strip_32$Lon_div_indx)
# [1] 6



df_Lon_strip_4 <- do.call(rbind, list(df_div_Lon_lst[[7]], df_div_Lon_lst[[8]]))
str(df_Lon_strip_4)
# 'data.frame':	6855 obs. of  10 variables:
unique(df_Lon_strip_4$Lon_div_indx)
# [1] 7 8


#---------------------------
# sort within each Lon strip
#---------------------------

# Lon: ascending (small-large)
# if Lat: descending (large-smal) not good for DSP, lat will have all neg values, not good for further B calculation
# Lat: ascending (small-large)
df_Lon_strip_1_sort <- arrange(df_Lon_strip_1, Lon, Lat)
str(df_Lon_strip_1_sort)

# 'data.frame':	3793 obs. of  10 variables:
#$ Lon               : num  -179 -179 -179 -179 -178 ...
#$ Lat               : num  66.8 67.5 68.2 71.2 66 ...

head(df_Lon_strip_1_sort, 30)




df_Lon_strip_2_sort <- arrange(df_Lon_strip_2, Lon, desc(Lat))
str(df_Lon_strip_2_sort)
# 'data.frame':	6581 obs. of  10 variables:
#$ Lon               : num  -89.2 -89.2 -89.2 -89.2 -89.2 ...
#$ Lat               : num  81.8 80.2 79.5 78.8 76.5 ...


df_Lon_strip_31_sort <- arrange(df_Lon_strip_31, Lon, desc(Lat))
str(df_Lon_strip_31_sort)
# 'data.frame':	5910 obs. of  10 variables:
#$ Lon               : num  0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 ...
#$ Lat               : num  52.5 51.8 51 49.5 48.8 ...


df_Lon_strip_32_sort <- arrange(df_Lon_strip_32, Lon, desc(Lat))
str(df_Lon_strip_32_sort)
# 'data.frame':	4245 obs. of  10 variables:
#$ Lon               : num  45.8 45.8 45.8 45.8 45.8 ...
#$ Lat               : num  68.2 66.8 66 65.2 64.5 ...

df_Lon_strip_4_sort <- arrange(df_Lon_strip_4, Lon, Lat)
str(df_Lon_strip_4_sort)
# 'data.frame':	6855 obs. of  10 variables:
#$ Lon               : num  90.8 90.8 90.8 90.8 90.8 ...
#$ Lat               : num  22.5 23.2 24 24.8 25.5 ...


#=======================================================
# Re-organize the BC,DU,etc.order to align with the DAG
#=======================================================

# 1: DU; 2: SU; 3: BC; 4: OM; 5:SS


#------------
# Lon_Strip_1
#------------
df_Lon_Strip_1_Sort <- data.frame(Lon = df_Lon_strip_1_sort$Lon,
           Lat = df_Lon_strip_1_sort$Lat,
           Year = df_Lon_strip_1_sort$Year,
           Z1 = df_Lon_strip_1_sort$DU_Residuals_log,
           Z2 = df_Lon_strip_1_sort$SU_Residuals_log,
           Z3 = df_Lon_strip_1_sort$BC_Residuals_log,
           Z4 = df_Lon_strip_1_sort$OM_Residuals_log,
           Z5 = df_Lon_strip_1_sort$SS_Residuals_log,
           PM25_res = df_Lon_strip_1_sort$PM25_Residuals_log,
           Lon_div_indx = df_Lon_strip_1_sort$Lon_div_indx)


str(df_Lon_Strip_1_Sort)
# 'data.frame':	3793 obs. of  10 variables:

df_Lon_Strip_1_Sort_new <- df_Lon_Strip_1_Sort

saveRDS(df_Lon_Strip_1_Sort_new, file = "df_Lon_Strip_1_Sort_new.rds")


range(df_Lon_Strip_1_Sort_new$Z1)
# [1] -3.308540  4.668558
quantile(df_Lon_Strip_1_Sort_new$Z1)
#     0%        25%        50%        75%        100% 
#  -3.3085399 -0.7825191 -0.2199278  0.5530506   4.6685585
 

range(df_Lon_Strip_1_Sort_new$Z2)
# [1] -1.227585  4.270820
quantile(df_Lon_Strip_1_Sort_new$Z2)
#          0%         25%         50%         75%        100% 
# -1.22758542 -0.24994444 -0.05963849  0.24968515  4.27082005 


range(df_Lon_Strip_1_Sort_new$Z3)
# [1] -2.302818  3.978425
quantile(df_Lon_Strip_1_Sort_new$Z3)
#        0%        25%        50%        75%       100% 
#-2.3028181 -0.5080320 -0.1693029  0.2157480  3.9784250 


range(df_Lon_Strip_1_Sort_new$Z4)
# [1] -2.318861  2.494690
quantile(df_Lon_Strip_1_Sort_new$Z4)
#         0%         25%         50%         75%        100% 
#-2.31886054 -0.38052471 -0.06471533  0.30779855  2.49468983 


range(df_Lon_Strip_1_Sort_new$Z5)
# [1] -2.192750  2.578429
quantile(df_Lon_Strip_1_Sort_new$Z5)
#          0%         25%         50%         75%        100% 
# -2.19275010 -0.42407572 -0.06511489  0.39396397  2.57842950 




#------------
# Lon_Strip_2
#------------
df_Lon_Strip_2_Sort <- data.frame(Lon = df_Lon_strip_2_sort$Lon,
                                  Lat = df_Lon_strip_2_sort$Lat,
                                  Year = df_Lon_strip_2_sort$Year,
                                  Z1 = df_Lon_strip_2_sort$DU_Residuals_log,
                                  Z2 = df_Lon_strip_2_sort$SU_Residuals_log,
                                  Z3 = df_Lon_strip_2_sort$BC_Residuals_log,
                                  Z4 = df_Lon_strip_2_sort$OM_Residuals_log,
                                  Z5 = df_Lon_strip_2_sort$SS_Residuals_log,
                                  PM25_res = df_Lon_strip_2_sort$PM25_Residuals_log,
                                  Lon_div_indx = df_Lon_strip_2_sort$Lon_div_indx)


str(df_Lon_Strip_2_Sort)
# 'data.frame':	6581 obs. of  10 variables:

saveRDS(df_Lon_Strip_2_Sort, file = "df_Lon_Strip_2_Sort.rds")


#------------
# Lon_Strip_31
#------------
df_Lon_Strip_31_Sort <- data.frame(Lon = df_Lon_strip_31_sort$Lon,
                                  Lat = df_Lon_strip_31_sort$Lat,
                                  Year = df_Lon_strip_31_sort$Year,
                                  Z1 = df_Lon_strip_31_sort$DU_Residuals_log,
                                  Z2 = df_Lon_strip_31_sort$SU_Residuals_log,
                                  Z3 = df_Lon_strip_31_sort$BC_Residuals_log,
                                  Z4 = df_Lon_strip_31_sort$OM_Residuals_log,
                                  Z5 = df_Lon_strip_31_sort$SS_Residuals_log,
                                  PM25_res = df_Lon_strip_31_sort$PM25_Residuals_log,
                                  Lon_div_indx = df_Lon_strip_31_sort$Lon_div_indx)


str(df_Lon_Strip_31_Sort)
# 'data.frame':	5910 obs. of  10 variables:
saveRDS(df_Lon_Strip_31_Sort, file = "df_Lon_Strip_31_Sort.rds")


#------------
# Lon_Strip_32
#------------

df_Lon_Strip_32_Sort <- data.frame(Lon = df_Lon_strip_32_sort$Lon,
                                   Lat = df_Lon_strip_32_sort$Lat,
                                   Year = df_Lon_strip_32_sort$Year,
                                   Z1 = df_Lon_strip_32_sort$DU_Residuals_log,
                                   Z2 = df_Lon_strip_32_sort$SU_Residuals_log,
                                   Z3 = df_Lon_strip_32_sort$BC_Residuals_log,
                                   Z4 = df_Lon_strip_32_sort$OM_Residuals_log,
                                   Z5 = df_Lon_strip_32_sort$SS_Residuals_log,
                                   PM25_res = df_Lon_strip_32_sort$PM25_Residuals_log,
                                   Lon_div_indx = df_Lon_strip_32_sort$Lon_div_indx)

str(df_Lon_Strip_32_Sort)
# 'data.frame':	4245 obs. of  10 variables:
saveRDS(df_Lon_Strip_32_Sort, file = "df_Lon_Strip_32_Sort.rds")


#------------
# Lon_Strip_4
#------------
df_Lon_Strip_4_Sort <- data.frame(Lon = df_Lon_strip_4_sort$Lon,
                                   Lat = df_Lon_strip_4_sort$Lat,
                                   Year = df_Lon_strip_4_sort$Year,
                                   Z1 = df_Lon_strip_4_sort$DU_Residuals_log,
                                   Z2 = df_Lon_strip_4_sort$SU_Residuals_log,
                                   Z3 = df_Lon_strip_4_sort$BC_Residuals_log,
                                   Z4 = df_Lon_strip_4_sort$OM_Residuals_log,
                                   Z5 = df_Lon_strip_4_sort$SS_Residuals_log,
                                   PM25_res = df_Lon_strip_4_sort$PM25_Residuals_log,
                                   Lon_div_indx = df_Lon_strip_4_sort$Lon_div_indx)


df_Lon_Strip_4_Sort_new <- df_Lon_Strip_4_Sort


str(df_Lon_Strip_4_Sort_new)
# 'data.frame':	6855 obs. of  10 variables:
saveRDS(df_Lon_Strip_4_Sort_new, file = "df_Lon_Strip_4_Sort_new.rds")

head(df_Lon_Strip_4_Sort_new, 30)

