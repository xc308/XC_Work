#============
# 20 Feb.2024
#============

# Aim:
  # Use the simu results from 032c to plot the corresponding plts for
  # SIGMA, SIGMA_inv, with B_original and B_SpN+Reg, 
  # SIGMA_inv with thres

# Method:
  # SG_SG_inv results from 032c
    # SG_SG_inv_6_a01d05_TriWave_Orig_Thres
    # SG_SG_inv_6_a01d05_Wend_Orig_Thres
    # SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres
    # SG_SG_inv_6_a01d05_Wend_SpNReg_Thres


# Settings:
  # p = 6
  # ds = 0.05, D = [-1, 1]
  # A = 0.1, dlt = 0.5


#==========
# Settings
#==========

#------------------------------------
# Location, displacements, distance
#------------------------------------
#ds <- 0.1 # for with SpN plts
ds <- 0.05 # for w/o SpN plts, esp. SIGMA; also try for plts with SpN
# ds = 0.05 has better visualization effect for both with, w/o SpN

s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20]; num [1:40]


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]; num [1:40, 1:40]


# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]; num [1:1600]


#----------------
# data structure
#----------------

p = 6
hierarchy_data6 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 3, 5))
)


#========
# Plots
#========

# B_Origin: Tri-Wave
plt_Sig(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA, p = 6)
#plt_Sig(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA_inv, p = 6)
plt_Sig(log(abs(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA_inv)), p = 6)


# B_Origin: Wendland
plt_Sig(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA, p = 6)
#plt_Sig(log(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA_inv), p = 6)
plt_Sig(log(abs(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA_inv)), p = 6)


# B_SpN+Reg: Tri-Wave
plt_Sig(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA, p = 6)
#plt_Sig(log(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA_inv), p = 6)
plt_Sig(log(abs(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA_inv)), p = 6)


# B_SpN+Reg: Wendland
plt_Sig(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA, p = 6)
#plt_Sig(log(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv), p = 6)
plt_Sig(log(abs(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv)), p = 6)


#--------------------------
# plot functions with Main
#--------------------------

Plot_SG_TW_Main <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma, atop("Tri-wave (Orig); p = 6")))
        main = expression(atop(Sigma, atop("Tri-wave (SpN + Reg)", 
                                           atop("p = 6; Reg = 1e-9"))))
  )
}

## Origin Tri-Wave SG
Plot_SG_TW_Main(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA, p = 6)

## SpN+Reg Tri-Wave SG
Plot_SG_TW_Main(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA, p = 6)




Plot_SG_Inv_TW_Main <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma^{-1}~ (log), atop("Tri-wave(Orig); Thres = 1e-3")))
        main = expression(atop(Sigma^{-1}~ (log), atop("Tri-wave (SpN + Reg)", 
                                           atop("Reg = 1e-9; Thres = 1e-3"))))
  )
}

# Origin Tri-Wav SG_inv
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA_inv)), p = 6)

## SpN+Reg Tri-Wave SG_inv
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA_inv)), p = 6)




Plot_SG_Main_Wend <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma, atop("Wendland (Orig); p = 6")))
        main = expression(atop(Sigma, atop("Wendland (SpN + Reg)", 
                                           atop("p = 6; Reg = 1e-9"))))
        
  )
}

# Origin Wend SG
Plot_SG_Main_Wend(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA, p = 6)
## SpN+Reg Wend SG_inv
Plot_SG_Main_Wend(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA, p = 6)



Plot_SG_Inv_Main_Wend <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 5.5, 1), cex.main = 2)
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        #main = expression(atop(Sigma^{-1}~ (log), atop("Wendland (Orig); Thres = 1e-3")))
        main = expression(atop(Sigma^{-1}~ (log), atop("Wendland (SpN + Reg)", 
                                                       atop("Reg = 1e-9; Thres = 1e-3"))))
        
  )
}

# # Origin Wend SG_inv
Plot_SG_Inv_Main_Wend(log(abs(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA_inv)), p = 6)

## SpN+Reg Wend SG_inv
Plot_SG_Inv_Main_Wend(log(abs(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv)), p = 6)


#-----------------------
# Side by side with main
#-----------------------

# Original 
jpeg(paste0(image.path, "SG_SGinv_TW_WL_Orig_sbs.jpeg"), 
     width = 10, height = 9, units = "in", res = 300)
par(mfrow = c(2, 2))

Plot_SG_TW_Main(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA, p = 6)
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_6_a01d05_TriWave_Orig_Thres$SIGMA_inv)), p = 6)
Plot_SG_Main_Wend(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA, p = 6)
Plot_SG_Inv_Main_Wend(log(abs(SG_SG_inv_6_a01d05_Wend_Orig_Thres$SIGMA_inv)), p = 6)

dev.off()




## SpN + Reg 
jpeg(paste0(image.path, "SG_SGinv_TW_WL_SpNReg_sbs.jpeg"), 
     width = 10, height = 9, units = "in", res = 300)
par(mfrow = c(2, 2))

Plot_SG_TW_Main(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA, p = 6)
Plot_SG_Inv_TW_Main(log(abs(SG_SG_inv_6_a01d05_TriWave_SpNReg_Thres$SIGMA_inv)), p = 6)

Plot_SG_Main_Wend(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA, p = 6)
Plot_SG_Inv_Main_Wend(log(abs(SG_SG_inv_6_a01d05_Wend_SpNReg_Thres$SIGMA_inv)), p = 6)


dev.off()







