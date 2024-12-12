#=============
# 12 Dec. 2024
#=============
# Aim:
  # plot the SG, SG_inv generated in 032e

# Obj:
  # SG_SG_inv_TriWave_p10_n400
  # SG_SG_inv_TriWave_p10_n800

# Method:
  # use plot functions in 032d


image.path <- "./Results/"



#==========================
# plot functions with Main
#==========================

#-------------
# Sigma plot
#-------------

Plot_SG_TW_Main <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma, atop("Tri-wave (Orig); p = 6")))
        main = expression(atop(Sigma, atop("Tri-wave (SpN + Reg)", 
                                           atop("p = 10; Reg = 1e-9"))))
  )
}


## SpN+Reg Tri-Wave SG
jpeg(paste0(image.path, "SG_TW_SpNReg_p10.jpeg"), 
     width = 10, height = 9, units = "in", res = 300)
#par(mfrow = c(2, 2))
par(mfrow = c(1, 1))
Plot_SG_TW_Main(SG_SG_inv_TriWave_p10_n40$SIGMA, p = 10)
dev.off()



#----------------
# Sigma_inv plot
#----------------
Plot_SG_Inv_TW_Main <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma^{-1}~ (log), atop("Tri-wave(Orig); Thres = 1e-3")))
        main = expression(atop(Sigma^{-1}~ (log), atop("Tri-wave (SpN + Reg)", 
                                                       atop("p = 10; Reg = 1e-9; Thres = 1e-3"))))
  )
}


## SpN+Reg Tri-Wave SG_inv
jpeg(paste0(image.path, "SG_inv_TW_p10.jpeg"),
     width = 10, height = 9, units = "in", res = 300)
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_TriWave_p10_n40$SIGMA_inv)), p = 10)
dev.off()


#-----------------------------
# SG, SG_inv side by side plot
#-----------------------------
jpeg(paste0(image.path, "SG_SG_inv_p10.jpeg"),
     width = 10, height = 9, units = "in", res = 300)
par(mfrow = c(2, 2))
Plot_SG_TW_Main(SG_SG_inv_TriWave_p10_n40$SIGMA, p = 10)
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_TriWave_p10_n40$SIGMA_inv)), p = 10)
dev.off()


getwd()







