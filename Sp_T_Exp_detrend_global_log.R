#******************************#
# Global Detrend on log scale
#******************************#

library(dplyr)


df_all # 136920 = 27384 * 5
head(df_all)

df_all_lg <- df_all
head(df_all_lg)

df_all_lg <- df_all_lg %>% mutate(across(.cols = c(BC:PM25), .fns = log))
head(df_all_lg)

# across apply a function/s to multiple cols of df



install.packages("stargazer")
library(stargazer)

install.packages("FRK")
library(FRK)

install.packages("sp")
library(sp)

library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(STRbook)


#=============================#
# Fit 6 globle models in one go
#=============================#

# global model fitted using step-wise selection, see below section

df_all_long_log <- gather(df_all_lg, key = Cmpts, value = Values, -Lon, -Lat, -Year, -ID )
str(df_all_long_log)



form_cmpts <- list()
lm_fit_cmpts <- list()
df_Cmpts <-list()
for(i in unique(df_all_long_log$Cmpts)) {
  
  df_Cmpts[[i]] <- filter(df_all_long_log, Cmpts == i)
  form_cmpts[[i]] <- as.formula(paste0(i, "~", "Lon + Lat + I(Lat ^ 2)"))
  lm_fit_cmpts[[i]] <- lm(form_cmpts[[i]], data = df_all)
  df_Cmpts[[i]]$Residual <- lm_fit_cmpts[[i]]$residuals
  
}

df_Cmpts_Res_long_log <- do.call("rbind", df_Cmpts) # 821520 * 7 = 136920 * 6 components * 7 var


#=====================#
# Save R obj as .RData
#=====================#

save(df_Cmpts_Res_long_log, file = "df_Cmpts_Res_long_log.RData")
load("df_Cmpts_Res_long_log.RData")



#=========================================#
# Stepwise selection model on global level
#=========================================#

BC_step_lm <- list()
for(i in 0:5){
  BC_step_lm[[i + 1]] <- step(lm(BC ~ 1, data = df_all_lg),
                              scope = BC ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2) ,
                              direction = "forward",
                              steps = i)
}

str(BC_step_lm) # List of 5


stargazer(BC_step_lm[[1]], BC_step_lm[[2]], 
          BC_step_lm[[3]], BC_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)

stargazer(BC_step_lm[[5]], BC_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          align = T,
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)



DU_step_lm <- list()
for(i in 0:5){
  DU_step_lm[[i + 1]] <- step(lm(DU ~ 1, data = df_all_lg),
                              scope = DU ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                              direction = "forward",
                              steps = i)
}


stargazer(DU_step_lm[[1]], DU_step_lm[[2]], 
          DU_step_lm[[3]], DU_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          dep.var.labels=c("DU (on log scale)"),
          order = c("Constant"),
          align = T,
          omit.stat = c("LL", "f"),
          no.space = T)


stargazer(DU_step_lm[[5]], DU_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          dep.var.labels=c("DU (on log scale)"),
          order = c("Constant"),
          align = T,
          omit.stat = c("LL", "f"),
          no.space = T)



OM_step_lm <- list()
for(i in 0:5){
  OM_step_lm[[i + 1]] <- step(lm(OM ~ 1, data = df_all_lg),
                              scope = OM ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                              direction = "forward",
                              steps = i)
}


stargazer(OM_step_lm[[1]], OM_step_lm[[2]], 
          OM_step_lm[[3]], OM_step_lm[[4]],
          title = "Stepwise Selection Results",
          dep.var.labels=c("OM (on log scale)"),
          order = c("Constant"),
          align = T,
          omit.stat = c("LL", "f"),
          no.space = T)


stargazer(OM_step_lm[[5]], OM_step_lm[[6]], 
          title = "Stepwise Selection Results",
          align = T,
          dep.var.labels=c("OM (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)



SU_step_lm <- list()
for(i in 0:5){
  SU_step_lm[[i + 1]] <- step(lm(SU ~ 1, data = df_all_lg),
                              scope = SU ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                              direction = "forward",
                              steps = i)
}


stargazer(SU_step_lm[[1]], SU_step_lm[[2]], 
          SU_step_lm[[3]], SU_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          dep.var.labels=c("SU (on log scale)"),
          order = c("Constant"),
          align = T,
          omit.stat = c("LL", "f"),
          no.space = T)


stargazer(SU_step_lm[[5]], SU_step_lm[[6]],
          title = "Stepwise Selection Results (5-6)",
          align = T,
          dep.var.labels=c("SU (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL",  "f"),
          no.space = T)




SS_step_lm <- list()
for(i in 0:5){
  SS_step_lm[[i + 1]] <- step(lm(SS ~ 1, data = df_all_lg),
                              scope = SS ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                              direction = "forward",
                              steps = i)
}


stargazer(SS_step_lm[[1]], SS_step_lm[[2]], 
          SS_step_lm[[3]], SS_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          dep.var.labels=c("SS (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)

stargazer(SS_step_lm[[5]], SS_step_lm[[6]],
          title = "Stepwise Selection Results (5-6)",
          align = T,
          dep.var.labels=c("SS (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)



PM25_step_lm <- list()
for(i in 0:5){
  PM25_step_lm[[i + 1]] <- step(lm(PM25 ~ 1, data = df_all_lg),
                                scope = PM25 ~ (Lon + Lat + Year)^2 + I(Lon^2) + I(Lat^2),
                                direction = "forward",
                                steps = i)
}

str(PM25_step_lm)

stargazer(PM25_step_lm[[1]], PM25_step_lm[[2]], 
          PM25_step_lm[[3]], PM25_step_lm[[4]],
          title = "Stepwise Selection Results (1-4)",
          align = T,
          dep.var.labels=c("PM25 (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)

stargazer(PM25_step_lm[[5]], PM25_step_lm[[6]], 
          title = "Stepwise Selection Results (5-6)",
          align = T,
          dep.var.labels=c("PM25 (on log scale)"),
          order = c("Constant"),
          omit.stat = c("LL", "f"),
          no.space = T)


