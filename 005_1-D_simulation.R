#==================================================
# 1-D simulation for Mixed Spatial Graphical Model 
#   (22 Aug 2023)
#==================================================

## Aim:
      # 1-D simulation for mixed spatial graphical model

rm(list = ls())


#==========
# Settings
#==========

image_path <- "./Results/"

install.packages("Matrix")
Yes
library(Matrix)

install.packages("dplyr")
library(dplyr)

#===================
# Simulation set up
#===================

#------------------
# Simulation domain
#------------------

# D = [-1, 1]
# spacing eta_i = 0.01, i = 1, 2, ..., 200

# n1: # of grid cells for process 1, Y_vec1 
# n2: # of grid cells for process 2, Y_vec2
# n = n1 + n2 = 400


#---------------------
# Process & Obs grids
#---------------------

#ds <- 0.01
ds <- 0.1 
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)

df <- data.frame(s = s, area = ds)

n1 <- n2 <- nrow(df)
n <- n1 + n2


#-----------
# Parameters
#-----------

# Matern for C11, and conditional variances C2_1, C3_2, C4_3, C5_4, C6_5, C7_6
sig2_11 <- sig2_21 <- sig2_32 <- sig2_43 <- sig2_54 <- 1 # marginal var for Matern
sig2_65 <- 1
sig2_76 <- 1
kappa11 <- kappa_21 <- kappa_32 <- kappa_43 <- kappa_54 <- 2
kappa_65 <- 2
kappa_76 <- 2



# Tri-Wave function
set.seed(50) 
delta <- sample(seq(-1, 1, by = 0.1), 1)
A = 1



#====================
# Matrix construction
#====================

#---------
# b() * ds
#---------

source("Fn_Wave_V4.R")

# displacements between pairs of points
H <- outer(df$s, df$s, FUN = "-")
H <- t(H)  

B <- A * wave_v4(h = H, delta = delta, A = 1) * ds
# [1:20, 1:20]


#-----------------
# D_vec: distance
#-----------------
D <- abs(H)
D_vec <- as.double(c(D)) #[1:400]


#-------------------------------
# C11 and conditional variances
#-------------------------------
source("Fn_Matern_32.R")
C11 <- Matern_32(Var = sig2_11, Kappa = kappa11, d_vec = D_vec) #[1:20, 1:20]
C2_1 <- Matern_32(Var = sig2_21, Kappa = kappa_21, d_vec = D_vec)
C3_2 <- Matern_32(Var = sig2_32, Kappa = kappa_32, d_vec = D_vec)
C4_3 <- Matern_32(Var = sig2_43, Kappa = kappa_43, d_vec = D_vec)
C5_4 <- Matern_32(Var = sig2_54, Kappa = kappa_54, d_vec = D_vec)
C6_5 <- Matern_32(Var = sig2_65, Kappa = kappa_65, d_vec = D_vec)
C7_6 <- Matern_32(Var = sig2_76, Kappa = kappa_76, d_vec = D_vec)


D <- bdiag(C11, C2_1, C3_2, C4_3, C5_4)
D <- bdiag(C11, C2_1, C3_2, C4_3, C5_4, C6_5)
D <- bdiag(C11, C2_1, C3_2, C4_3, C5_4, C6_5, C7_6)
str(D)


#===========
# Algorithm
#===========

p <- 5  # the number of variates
p <- 6
p <- 7

#-----------------------------
# Hierarchical data structure
#-----------------------------

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)

hierarchy_data2 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5))
)

cond <- hierarchy_data2$node_id == 6
hierarchy_data2[cond, ]$par_id
Check_par_node(Node = 6, data = hierarchy_data2) #[1] 1 5


hierarchy_data3 <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5, 6, 6, 7, 7),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4, c(1, 5), c(5, 3))
)





#--------
# Algo
#-------
source("Fn_Check_par_node.R")
Check_par_node(Node = 2, data = hierarchy_data) # [1] 1
Check_par_node(Node = 3, data = hierarchy_data) # [1] 2 1
Check_par_node(Node = 4, data = hierarchy_data) # [1] 2 3
Check_par_node(Node = 5, data = hierarchy_data)
Check_par_node(Node = 6, data = hierarchy_data2)
Check_par_node(Node = 7, data = hierarchy_data3) # [1] 5 3



n <- nrow(H) #20

make_SIGMA <- function(p, data) {
  
  ## bivariate
  B21 <- B
  C21 <- B21 %*% C11
  C12 <- t(C21)
  C22 <- B21 %*% t(C21) + C2_1
  
  SIGMA <- cbind(rbind(C11, C21), rbind(C12, C22))
  
  if (p == 2) {return (SIGMA)}
 
    
    for(r in seq(3, p, by = 1)) {
      
      PN = Check_par_node(Node = r, data = data)
      R <- C <- NULL
      
      for (c in seq(1, (r-1), by = 1)) {
        
        ## highly multivariate
        
        # C_rc = Sum_t B_rt C_tc
        BT <- NULL
        set.seed(0823)  # Brt remains the same for different c
        C_rc <- matrix(0, nrow = nrow(B), ncol = ncol(B))
        for (t in c(PN)) {
          
          delta <- sample(seq(-1, 1, by = 0.1), 1)
          B_rt <- A * wave_v4(h = H, delta = delta, A = 1) * ds
          
          BT <- rbind(BT, t(B_rt))
          
          C_rc <- C_rc + B_rt %*% SIGMA[((t - 1) * n + 1) : (t * n), ((c - 1) * n + 1) : (c * n)]
        }
        
        R <- cbind(R, C_rc)
        
        C_cr <- t(C_rc)
        C <- rbind(C, C_cr)
      }
      
      # C_rr <- Sum_t C_rt B_rt^T + D_rr
      t <- c(PN)
      Subset_cols <- function(t) {
        start_col <- (t - 1) * n + 1
        end_col <- t * n
        
        result <- R[, start_col:end_col]
      }
      
      result_lst <- lapply(t, FUN = Subset_cols)
      R_subset <- do.call(cbind, result_lst)
      
      C_rr <- R_subset %*% BT + D[((r-1)*n + 1) : (r*n), ((r-1)*n + 1) : (r*n)]
      
      Col <- rbind(C, C_rr)
      Row <- rbind(SIGMA, R)
      SIGMA <- cbind(Row, Col)
      
      
      if (r == p) return(as.matrix(SIGMA))
    }
    
}



#========
# Test
#========
SIGMA_1212 <- make_SIGMA(p = 2, data = hierarchy_data2)
SIGMA_1313 <- make_SIGMA(p = 3, data = hierarchy_data2)
SIGMA_1414 <- make_SIGMA(p = 4, data = hierarchy_data2)
SIGMA_1515 <- make_SIGMA(p = 5, data = hierarchy_data2)
SIGMA_1616 <- make_SIGMA(p = 6, data = hierarchy_data2)
SIGMA_1717 <- make_SIGMA(p = 7, data = hierarchy_data3)



#~~~~~~~~~~
# sym & pd
#~~~~~~~~~~

source("Fun_Tst_Sym_Pd.R")
Test_sym_pd(SIGMA_1414)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

Test_sym_pd(SIGMA_1515)
Test_sym_pd(SIGMA_1616)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

Test_sym_pd(SIGMA_1717)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


Test_sym_pd(SIGMA_Chain_7)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"



#=======================
# Visualisation of SIGMA
#=======================

#------
# Opt1
#------
plt_Sig <- function(Sigma, p) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  par(mar = c(3, 3, 3, 1))
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        main = paste("p = ", p))
}

plt_Sig(SIGMA_1212, p = 2)
plt_Sig(SIGMA_1313, p = 3)
plt_Sig(SIGMA_1414, p = 4)
plt_Sig(SIGMA_1515, p = 5)
plt_Sig(SIGMA_1616, p = 6)
plt_Sig(SIGMA_1717, p = 7)




#--------
# Save
#--------

png(paste0(image_path, "SIGMAS.png"), res = 300,
    width = 12, height = 10, units = "in")

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))

plt_Sig(SIGMA_1212, p = 2)
plt_Sig(SIGMA_1313, p = 3)
plt_Sig(SIGMA_1414, p = 4)
plt_Sig(SIGMA_1515, p = 5)

dev.off()


png(paste0(image_path, "SIGMA1212.png"), res = 300,
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_1212, p = 2)
dev.off()


png(paste0(image_path, "SIGMA1313.png"), res = 300, 
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_1313, p = 3)
dev.off()


png(paste0(image_path, "SIGMA1414.png"), res = 300, 
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_1414, p = 4)
dev.off()


png(paste0(image_path, "SIMGA1515.png"), res = 300, 
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_1515, p = 5)
dev.off()


png(paste0(image_path, "SIGMA1616.png"), res = 300, 
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_1616, p = 6)
dev.off()

#--------------
# piece together
#--------------


require("cowplot")
install.packages("cowplot")
library(cowplot)
plot_grid(S2, S3, nrow = 1)



#==================
# Chain structure
#==================

## Data
hierarchy_data4 <- data.frame(
  node_id = c(1, 2, 3, 4, 5, 6, 7),
  par_id = c(NA, 1, 2, 3, 4, 5, 6)
)


## Construct
SIGMA_Chain_7 <- make_SIGMA(p = 7, data = hierarchy_data4)


## Test
Test_sym_pd(SIGMA_Chain_7)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


## plot and save
png(paste0(image_path, "SIGMA7_Chain.png"), res = 300, 
    width = 8, height = 7, units = "in")
plt_Sig(SIGMA_Chain_7, p = 7)
dev.off()


## Inverse
SIGMA_7_C_inv <- chol2inv(chol(SIGMA_Chain_7))
SIGMA_7_C_inv_mat <- as.matrix(SIGMA_7_C_inv)
plt_Sig(SIGMA_7_C_inv_mat, p = 7)

png(paste0(image_path, "SIGMA7_Chain_Inv_log.png"),
    res = 300, width = 8, height = 7, units = "in")
plt_Sig(log(abs(SIGMA_7_C_inv_mat)), p = 7)
dev.off()














