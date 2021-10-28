#=======================#
# Univariate correlation
#=======================#

# Uni-var
# Lag0 only

Lag0_corr_BC <- cor(X[["BC"]])
save(Lag0_corr_BC, file = "Lag0_corr_BC")
# save(/CorrData/Lag0_corr_BC, file = "Lag0_corr_BC")
rm(Lag0_corr_BC)

Lag0_corr_DU <- cor(X[["DU"]])
save(Lag0_corr_DU, file = "Lag0_corr_DU")
rm(Lag0_corr_DU)

Lag0_corr_OM <- cor(X[["OM"]])
save(Lag0_corr_OM, file = "Lag0_corr_OM")
rm(Lag0_corr_OM)

Lag0_corr_SU <- cor(X[["SU"]])
save(Lag0_corr_SU, file = "Lag0_corr_SU")
rm(Lag0_corr_SU)

Lag0_corr_SS <- cor(X[["SS"]])
save(Lag0_corr_SS, file = "Lag0_corr_SS")
rm(Lag0_corr_SS)


Lag0_corr_PM25 <- cor(X[["PM25"]])
save(Lag0_corr_PM25, file = "Lag0_corr_PM25")
rm(Lag0_corr_PM25)





## BC VS DU
Lag0_corr_BC_DU <- cor(X[["BC"]], X[["DU"]])
save(Lag0_corr_BC_DU, file = "Lag0_corr_BC_DU")
rm(Lag0_corr_BC_DU)

## DU VS SU
Lag0_corr_DU_SU <- cor(X[["DU"]], X[["SU"]])
save(Lag0_corr_DU_SU, file = "Lag0_corr_DU_SU")
rm(Lag0_corr_DU_SU)

## SU VS PM25
Lag0_corr_SU_PM25 <- cor(X[["SU"]], X[["PM25"]])
save(Lag0_corr_SU_PM25, file = "Lag0_corr_SU_PM25")
rm(Lag0_corr_SU_PM25)

## PM VS OM
Lag0_corr_PM25_OM <- cor(X[["PM25"]], X[["OM"]])
save(Lag0_corr_PM25_OM, file = "Lag0_corr_PM25_OM")
rm(Lag0_corr_PM25_OM)

## OM VS SS
Lag0_corr_OM_SS <- cor(X[["OM"]], X[["SS"]])
save(Lag0_corr_OM_SS, file = "Lag0_corr_OM_SS")
rm(Lag0_corr_OM_SS)
