#*********#
# Detrend
#*********#

# prepare for the Lag-0, and lat-1 covariance or cross-cov

# ref: https://bookdown.org/ronsarafian/IntrotoDS/lme.html


df_all # 136920 = 27384 * 5
df_all_long # 821520 = 136920 * 6 components
head(df_all)

install.packages("stargazer")
library(stargazer)

install.packages("FRK")
library(FRK)

install.packages("sp")
library(sp)





#===================#
# Stepwise selection
#===================#

BC_step_lm <- list()
for(i in 0:5){
   BC_step_lm[[i + 1]] <- step(lm(BC ~ 1, data = df_all),
        scope = BC ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2) ,
        direction = "forward",
        steps = i)
}

str(BC_step_lm) # List of 5



stargazer(BC_step_lm[[1]], BC_step_lm[[2]], 
          BC_step_lm[[3]], BC_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)

          #covariate.labels = c("Lat", "Lon", "Year", "Lat:Lon", "Intercept"),
          #order = c("Intercept"))

stargazer(BC_step_lm[[5]], BC_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)




DU_step_lm <- list()
for(i in 0:5){
   DU_step_lm[[i + 1]] <- step(lm(DU ~ 1, data = df_all),
                               scope = DU ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                               direction = "forward",
                               steps = i)
}


stargazer(DU_step_lm[[1]], DU_step_lm[[2]], 
          DU_step_lm[[3]], DU_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)


stargazer(DU_step_lm[[5]], DU_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)




OM_step_lm <- list()
for(i in 0:5){
   OM_step_lm[[i + 1]] <- step(lm(OM ~ 1, data = df_all),
                               scope = OM ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                               direction = "forward",
                               steps = i)
}


stargazer(OM_step_lm[[1]], OM_step_lm[[2]], 
          OM_step_lm[[3]], OM_step_lm[[4]],
          title = "Stepwise Selection Results",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)


stargazer(OM_step_lm[[5]], OM_step_lm[[6]], 
          title = "Stepwise Selection Results",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)





SU_step_lm <- list()
for(i in 0:5){
   SU_step_lm[[i + 1]] <- step(lm(SU ~ 1, data = df_all),
                               scope = SU ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                               direction = "forward",
                               steps = i)
}


stargazer(SU_step_lm[[1]], SU_step_lm[[2]], 
          SU_step_lm[[3]], SU_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)


stargazer(SU_step_lm[[5]], SU_step_lm[[6]],
          title = "Stepwise Selection Results (5-6)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)





SS_step_lm <- list()
for(i in 0:5){
   SS_step_lm[[i + 1]] <- step(lm(SS ~ 1, data = df_all),
                               scope = SS ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                               direction = "forward",
                               steps = i)
}


stargazer(SS_step_lm[[1]], SS_step_lm[[2]], 
          SS_step_lm[[3]], SS_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)

stargazer(SS_step_lm[[5]], SS_step_lm[[6]],
          title = "Stepwise Selection Results (5-6)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)



PM25_step_lm <- list()
for(i in 0:5){
   PM25_step_lm[[i + 1]] <- step(lm(PM25 ~ 1, data = df_all),
                               scope = PM25 ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                               direction = "forward",
                               steps = i)
}

str(PM25_step_lm)

stargazer(PM25_step_lm[[1]], PM25_step_lm[[2]], 
          PM25_step_lm[[3]], PM25_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)

stargazer(PM25_step_lm[[5]], PM25_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          align = T,
          omit.stat = c("LL", "ser", "f"),
          no.space = T)



#=========================#
# 'Add lable col' to df_all
#=========================#

head(df_all)
str(df_all) # 136920 = 27384 * 5

df_all_long <- gather(df_all, key = Cmpts, value = Values, -Lon, -Lat, -Year, -ID )
str(df_all_long) # 821520 = 136920 * 6 components
head(df_all_long)
tail(df_all_long)



#=======================#
# Fit 6 models in one go
#=======================#

form_cmpts <- list()
lm_fit_cmpts <- list()
df_Cmpts <-list()
for(i in unique(df_all_long$Cmpts)) {
   
   df_Cmpts[[i]] <- filter(df_all_long, Cmpts == i)
   form_cmpts[[i]] <- as.formula(paste0(i, "~", "Lon + Lat + I(Lat ^ 2)"))
   lm_fit_cmpts[[i]] <- lm(form_cmpts[[i]], data = df_all)
   df_Cmpts[[i]]$Residual <- lm_fit_cmpts[[i]]$residuals
  
}

df_Cmpts_Res_long <- do.call("rbind", df_Cmpts) # 821520 * 7 = 136920 * 6 components * 7 var



#=====================#
# Save R obj as .RData
#=====================#

save(df_Cmpts_Res_long, file = "df_Cmpts_Res_long.RData")
load("df_Cmpts_Res_long.RData")


saveRDS(df_Cmpts_Res_long, file = "df_Cmpts_Res_long.RDS")
readRDS("df_Cmpts_Res_long.RDS")



head(df_Cmpts_Res_long, 3)
#        Lon   Lat Year ID Cmpts      Values  Residual
#BC.1 -42.75 83.25 2016  1    BC 0.007900851 0.3067757
#BC.2 -42.00 83.25 2016  2    BC 0.006894099 0.3046130
#BC.3 -39.75 83.25 2016  3    BC 0.006666577 0.3009178

tail(df_Cmpts_Res_long, 3)
#                Lon    Lat Year    ID Cmpts   Values Residual
#PM25.136918 -66.00 -54.75 2012 27382  PM25 5.111783 41.65979
#PM25.136919 -65.25 -54.75 2012 27383  PM25 5.291296 41.79623
#PM25.136920 -68.25 -55.50 2012 27384  PM25 5.179908 43.27677







##----------------##
str(df_Cmpts)
str(df_Cmpts[["BC"]])
str(form_cmpts)
str(lm_fit_cmpts[["BC"]])

df_Cmpts[["BC"]]
form_cmpts[["BC"]]
##----------------##



##########################################################

#========================#
# To add more covariates
#========================#

getwd()
load("GM_dat.RData")
str(GM_dat)

GM_dat_5yr <- GM_dat %>% filter(Year == c(2012, 2013, 2014, 2015, 2016))
# 5143

GM_dat_5yr_selt <- GM_dat_5yr %>% select(Year, Longitude, Latitude, MonitorType, POP, ELEVATION)
str(GM_dat_5yr_selt)

colnames(GM_dat_5yr_selt[2]) <- "Lon"
colnames(GM_dat_5yr_selt[3]) <- "Lat"


# head(GM_dat_5yr_selt)





#===========#
# Split data
#===========#

df_train_idx <- sample(unique(df_all$ID), length(unique(df_all$ID)) * 0.8)
# 21907

df_train <- df_all[ID %in% df_train_idx, ] # 109535
df_test <- df_all[!ID %in% df_train_idx, ] # 27385

par(mfrow = c(1, 1))
par(mai = rep(0.8, 4))

plot(unique(df_train[, c("Lon", "Lat")]), col = "LightSkyBlue",
     pch = 1, cex = 0.3)
points(unique(df_test[, c("Lon", "Lat")]), col = "darkblue",
       cex = 0.15)


legend("top", inset = c(200, -0.3), xpd = T, horiz = T,
       pch = 1, col = c("LightSkyBlue", "darkblue"),
       legend = c("Train", "Test"), cex = 0.8, bty = "n")


# ref: https://r-coder.com/add-legend-r/#Change_legend_size



#========================#
# Detrend Model seletion
#=======================#

#-------------------------#
# Linear fix-effect model
#-------------------------#

# model 1
mod_BC_1 <- lm(BC ~ Lon + Lat + I(Lat ^ 2) + Year + I(Year ^ 2), 
   data = df_train)

mod_BC_1_pred <- predict(mod_BC_1, newdata = df_test)
str(mod_BC_1_pred) # 27385

mean((df_test$BC - mod_BC_1_pred) ^ 2) # [1] 1.032977

test.response <- eval(mod_BC_1, envir = df_test)
str(test.response )

mean((test.response - mod_BC_1_pred)^2)

all(mod_BC_1[[2]] == mod_BC_1$residuals) # T
mod_BC_1$residuals

## model 2
mod_BC_2 <- lm(BC ~ Lon + Lat + I(Lat ^ 2) + cos(Year), 
   data = df_train)

mod_BC_2_pred <- predict(mod_BC_2, newdata = df_test)

mean((df_test$BC - mod_BC_2_pred) ^ 2) # [1] 1.033209



mod_BC_5 <- lm(BC ~ Lon + Lat + I(Lat ^ 2) + Year, 
               data = df_train)
mod_BC_5_pred <- predict(mod_BC_5, newdata = df_test)

mean((df_test$BC - mod_BC_5_pred) ^ 2) # [1] 1.033038


#---------------#
# Random model
#---------------#

# rd intercept model

install.packages("lme4")
library(lme4)


mod_BC_3 <- lmer(BC ~ Lon + Lat + Year + (1 | ID), 
               data = df_train)

mod_BC_3_pred <- predict(mod_BC_3, newdata = df_test, allow.new.levels = T)
mean((df_test$BC - mod_BC_3_pred) ^ 2) # [1] 1.089686


## rd intercept rd slope model
mod_BC_4 <- lmer(BC ~ Lon + Lat + Year + ( 1 + Year | ID), 
                 data = df_train)

## boundary (singular) fit: see ?isSingular



#====================#
# Plot of Lon_strips
#====================#

coords <- unique(df_all[, c("ID", "lon_strips", "Lon", "Lat")])
head(coords)

ggplot(data = coords, aes(x = Lon, y = Lat, color = as.factor(lon_strips))) + 
  geom_point()



#=====================#
# construct some basis
#=====================#

G <- auto_basis(data = df_all[, c("Lon", "Lat")] %>% SpatialPoints(),
                nres = 2, type = "Gaussian")

show_basis(G)

S <- eval_basis(basis = G, s = df_all[, c("Lon", "Lat")] %>% as.matrix()) %>%
   as.matrix()

colnames(S) <- paste0("B", 1:ncol(S))
head(S) # B1: B24

df_all_bas <- cbind(df_all, S)

rm(df_all_bas)
rm(S)


