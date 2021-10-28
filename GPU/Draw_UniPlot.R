
#strip_num <- i

#strip_num <- 1
#strip_num <- 2
#strip_num <- 3

#strip_num_input <- 4



draw_fun <- function(strip_num_input){

load("Lag0_corr_BC")
png(paste0("Lag0_corr_BC-0",strip_num_input,".png"))
#png(paste0("Figures/Lag0_corr_BC-0", strip_num_input, ".png"))
# will create folder Figures automatically

EMP_cor_lat_plt2(Lag0_corr_BC, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_BC)

load("Lag0_corr_DU")
png(paste0("Lag0_corr_DU-0",strip_num_input,".png"))
EMP_cor_lat_plt2(Lag0_corr_DU, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_DU)

load("Lag0_corr_OM")
png(paste0("Lag0_corr_OM-0",strip_num_input,".png"))
EMP_cor_lat_plt2(Lag0_corr_OM, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_OM)

load("Lag0_corr_SS")
png(paste0("Lag0_corr_SS-0",strip_num_input,".png"))
EMP_cor_lat_plt2(Lag0_corr_SS, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_SS)

load("Lag0_corr_SU")
png(paste0("Lag0_corr_SU-0",strip_num_input,".png"))
EMP_cor_lat_plt2(Lag0_corr_SU, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_SU)

load("Lag0_corr_PM25")
png(paste0("Lag0_corr_PM25-0",strip_num_input,".png"))
EMP_cor_lat_plt2(Lag0_corr_PM25, sp_loc_df, strip_num = strip_num_input)
dev.off()
rm(Lag0_corr_PM25)
}

for(i in 1:4){
    draw_fun(strip_num_input = i)
}

