#***************#
# Sp_T_Exp_EOF
#***************#

head(Comps_lst[[1]])
tail(Comps_lst[[1]])
str(Comps_lst[[1]])
str(Comps_lst)


#====================#
# Temporal-wide Data
#====================#

tt <- Comps_lst[[1]][, -c(5, 6)]
head(tt)

dim(Comps_lst[[1]]) # 136920 * 6
tp_wd <- Comps_lst[[1]][, -c(5, 6)] %>% 
  spread(key = Year, value = BC)

head(tp_wd)
tail(tp_wd)


tp_wd2 <- tp_wd[ ,-c(1, 2)]
head(tp_wd2)

Z <- t(tp_wd2)
str(Z) # num [1:5, 1:27384]


#=======================#
# Detrend & Rescale
#=======================#

# Spatial Mean over time
spat_mean <- apply(tp_wd, 1, mean) # a vector 1:27384
nT <- ncol(tp_wd[ ,-c(1, 2)]) # 5

Zspat_detrend <- Z - outer(rep(1, nT), spat_mean) # 5 * 27384
str(Zspat_detrend) # num [1:5, 1:27384]


# standardize
Zt <- 1/sqrt(nT - 1) * Zspat_detrend



#========#
# SVD Zt
#========#

E <- svd(Zt)
# V: constains EOFs in temp-wide format

V <- E$v
str(V) # num [1:27384, 1:5] 
colnames(V) <- paste0("EOF", 1:ncol(tp_wd2))
head(V)

EOFs <- cbind(tp_wd[, c(1, 2)], V) # possible problem: order of lon lat
head(EOFs)


#============#
# Visualize
#============#

ggplot(EOFs) + 
  geom_tile(aes(x = Lon, y = Lat, fill = EOF1)) + 
  fill_scale(name = "mug/m^3") +
  theme_bw() + 
  xlab("Longitude") + 
  ylab("Latitude")

ggplot(EOFs) + 
  geom_tile(aes(x = Lon, y = Lat, fill = EOF4)) + 
  fill_scale(name = "mug/m^3") +
  theme_bw() + 
  xlab("Longitude") + 
  ylab("Latitude")


