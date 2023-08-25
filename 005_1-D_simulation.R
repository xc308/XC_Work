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

image_path <- "./Results/Fig/"

install.packages("Matrix")
Yes
library(Matrix)


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

# Matern for C11, and conditional variances C2_1, C3_2, C4_3, C_54
sig2_11 <- sig2_21 <- sig2_32 <- sig2_43 <- sig2_54 <- 1 # marginal var for Matern
kappa11 <- kappa_21 <- kappa_32 <- kappa_43 <- kappa_54 <- 2

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


D <- bdiag(C11, C2_1, C3_2, C4_3, C5_4)
str(D)


#===========
# Algorithm
#===========

p <- 5  # the number of variates


#-----------------------------
# Hierarchical data structure
#-----------------------------

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 3, 4, 4, 5),
  par_id = c(NA, 1, c(2, 1), c(2, 3), 4)
)


#--------
# Algo
#-------
source("Fn_Check_par_node.R")
Check_par_node(Node = 2, data = hierarchy_data) # [1] 1
Check_par_node(Node = 3, data = hierarchy_data) # [1] 2 1
Check_par_node(Node = 4, data = hierarchy_data) # [1] 2 3
Check_par_node(Node = 5, data = hierarchy_data)


n <- nrow(H) #20

make_SIGMA <- function(p) {
  
  ## bivariate
  B21 <- B
  C21 <- B21 %*% C11
  C12 <- t(C21)
  C22 <- B21 %*% t(C21) + C2_1
  
  SIGMA <- cbind(rbind(C11, C21), rbind(C12, C22))
  
  if (p == 2) {return (SIGMA)}
 
    
    for(r in seq(3, p, by = 1)) {
      
      PN = Check_par_node(Node = r, data = hierarchy_data)
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
SIGMA_1212 <- make_SIGMA(p = 2)
SIGMA_1313 <- make_SIGMA(p = 3)
SIGMA_1414 <- make_SIGMA(p = 4)
SIGMA_1515 <- make_SIGMA(p = 5)


source("Fun_Tst_Sym_Pd.R")
Test_sym_pd(SIGMA_1414)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"

Test_sym_pd(SIGMA_1515)


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



#------
# Opt2
#------

plt_Sigma <- function(Sigma, p) {
  #p = 4
  indx <- seq(1, p, by = 1)
  
  Sigma_df <- expand.grid(s1 = df$s, comp1 = c(paste0("Y", indx)), 
                          s2 = df$s, comp2 = c(paste0("Y", indx))) %>%
    mutate(cov = c(Sigma))
  
  Sigma_plt <- ggplot() + 
    geom_tile(data = Sigma_df, aes(x = s1, s2, fill = cov)) + 
    facet_grid(comp1 ~ comp2) +
    scale_fill_gradient(low = "#FFFFCC", high = "#000080") +
    ylab("s") + xlab("u") + 
    scale_y_reverse()
  
  print(Sigma_plt)
  
}


















