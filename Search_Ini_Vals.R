#=============
# 23 May 2023
#=============

# Aim:
  # Test initial values for neg_logL 
  # to avoid Inf for 064a

# Method:
  # local evaluation of neg_logL
  # grid search
  # preventative coding in neg_logL function .Machine$double.xmax

# Function:
  # 060_2D_Inf_neg_logL_CAR_GPU

# Data:
  # df_Lon_Strip_1_new in 063



library(Matrix)

#==============
# GPU settings
#==============
#------
# torch
#------
# Set the library path to the desired directory
.libPaths("/bask/projects/v/vjgo8416-xchen")

# Load the torch library
library(torch)


#torch_set_num_interop_threads(2)
#torch_set_num_threads(2)

cat("inter threads:", "\n")
torch_get_num_interop_threads()

cat("intra threads:", "\n")
torch_get_num_threads()


#-------------------------------
# check BLAS and OPENBLAS info: only 36 works
#-------------------------------
#install.packages("RhpcBLASctl")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(RhpcBLASctl)


#cat("Check Current BLAS Library", "\n")
#sessionInfo()

cat("Check the current number of BLAS threads", "\n")
blas_get_num_procs()

#blas_set_num_threads(48)

#cat("Updated BLAS threads:", "\n")
#blas_get_num_procs()



#-----------
# GPUmatrix
#------------

#install.packages("GPUmatrix", lib="/bask/projects/v/vjgo8416-xchen")
.libPaths("/bask/projects/v/vjgo8416-xchen")
library(GPUmatrix)

system("nvidia-smi")



#==============
# Pre-settings
#==============

#-----
# df
#-----
df_Lon_Strp_1_Srt <- readRDS("df_Lon_Strip_1_Sort_new.rds")


#---------
# data str
#---------
hierarchy_data_CAMS <- data.frame(
  node_id = c(1, 2, 3, 4,  5, 5),
  par_id = c(NA, 1, 2, 3, c(4, 1))
)

p = 5
data_str <- hierarchy_data_CAMS


#-----
# pars
#-----
source("Fn_para_mat_construct.R")
all_pars_lst_CAR_2D_CMS <- All_paras_CAR_2D(p = 5, data = hierarchy_data_CAMS)
all_pars_lst <- all_pars_lst_CAR_2D_CMS


#----------
# 2D coords
#----------
crds <- cbind(df_Lon_Strp_1_Srt$Lon, df_Lon_Strp_1_Srt$Lat)
#head(crds)


#----------------------------------------
# Construct displacement matrix (DSP_mat)
#----------------------------------------
# Aim:
# for TST12 function to construct shft dist matrix

source("Fn_make_DSP_mat.R")
DSP <- make_DSP_mat(crds = crds)
#str(DSP[, , 1]) # num [1:3793, 1:3793]
#str(DSP[, , 2]) # num [1:3793, 1:3793]

#DSP[, , 1][1, 1:20]  #for dlt_lon ini
# [1] 0.00 0.00 0.00 0.00 0.75
#[6] 0.75 0.75 0.75 0.75 1.50
#[11] 1.50 1.50 1.50 2.25 2.25
#[16] 2.25 3.00 3.00 3.00 3.75

#DSP[, , 2][1, 1:20] #for dlt_lat ini
# [1]  0.00  0.75  1.50  4.50
#[5] -0.75  0.00  0.75  1.50
#[9]  4.50 -0.75  0.00  0.75
#[13]  4.50 -0.75  0.00  0.75
#[17] -0.75  0.00  0.75 -1.50



#---------------------------
# Construct distance matrix 
#---------------------------
# Aim:
# for H_adj and phi for UniCAR

DIST <- as.matrix(dist(crds, diag = T, upper = T))

Nb_radius <- 1.5 # degree 

H_adj <- matrix(as.numeric(abs(DIST) < Nb_radius), nrow(DIST), nrow(DIST))
diag(H_adj) <- 0

#dim(H_adj) # 3793 3793
#length(which(H_adj != 0)) # 27772

#spec <- eigen(H_adj, symmetric = T, only.values = T)$val
#max(abs(spec)) # [1] 7.976106

#phi <- 1/max(abs(spec)) # [1] 0.1253745
#phi <- trunc(phi * 100)/100 # [1] 0.12
phi <- 0.12

#=========
# ini_vals
#=========
#Vals <- c(0.7, 0.8, 0.8, 0.65, 0.65,    #A 
 #         0.2, 0.3, 0.3, 0.4, 0.4, 0.5, # dlt_lon
#          rep(0.2, 5),                  # dlt_lat
 #         rep(0.5, 5)                   # sig2
  #        ) # w/o tau2s


#all_ini_Vals <- c(Vals, rep(1, p)) # with tau2s



#=============================================
# Grid search pars for non-Inf neg_logL value
#=============================================
source("Fn_neg_logL_CAR_2D_GPU.R")

# possible pars choices
A_vals <- c(1, 2, 10)
dlt_lon_vals <- c(5, 10, 20)
dlt_lat_vals <- c(5, 10, 20)
sig2_vals <- 0.5
tau2_vals <- 0.5


# Initialize minimum negative log-likelihood to a large value
min_neg_ll <- Inf
best_params <- c() 
for(A in A_vals) {
  cat("A", A)
  for(dlt_lon in dlt_lon_vals) {
    cat("dlt_lon", dlt_lon)
    for(dlt_lat in dlt_lat_vals) {
      cat("dlt_lat", dlt_lat)
      for(sig2 in sig2_vals){
        cat("sig2", sig2)
        for(tau2 in tau2_vals){
            cat("tau2", tau2)
            all_ini_Vals <- c(rep(A, p), rep(dlt_lon, p), rep(dlt_lat, p),
                              rep(sig2, p), rep(tau2, p))
            
            
            neg_ll <- neg_logL_CAR_2D_GPU(theta = all_ini_Vals, p = 5, 
                                data_str = hierarchy_data_CAMS,
                                all_pars_lst = all_pars_lst_CAR_2D_CMS,
                                dsp_lon_mat = DSP[, , 1],
                                dsp_lat_mat = DSP[, , 2],
                                b = "Tri-Wave",
                                phi = phi, H_adj = H_adj, 
                                df = df_Lon_Strp_1_Srt)
            
            neg_ll <- as.numeric(neg_ll, device = "cpu")
            
            # Check if current neg-log likelihood is less than the minimum found so far
            if (neg_ll < min_neg_ll) {
              min_neg_ll <- neg_ll
              best_params <- all_ini_Vals
              cat("Yes", "\n")
            } else {
              cat("No", "\n")
            }
            
          }
        }
      }
    }
}

cat("best pars:", "\n")  
best_params

cat("neg_ll:", "\n")  
min_neg_ll  
  
  











 













