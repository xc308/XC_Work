#*********#
# Detrend
#*********#

# prepare for the Lag-0, and lat-1 covariance or cross-cov

# ref: https://bookdown.org/ronsarafian/IntrotoDS/lme.html


df_all # 136920 = 27384 * 5
df_all_long # 821520 = 136920 * 6 components


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


