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
  # 055

# Data:
  # df_Lon_Strip_1_new in 063


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

spec <- eigen(H_adj, symmetric = T, only.values = T)$val
#max(abs(spec)) # [1] 7.976106

phi <- 1/max(abs(spec)) # [1] 0.1253745
phi <- trunc(phi * 100)/100 # [1] 0.12


#=========
# ini_vals
#=========
Vals <- c(0.7, 0.8, 0.8, 0.65, 0.65,    #A 
          0.2, 0.3, 0.3, 0.4, 0.4, 0.5, # dlt_lon
          rep(0.2, 5),                  # dlt_lat
          rep(0.5, 5)                   # sig2
          ) # w/o tau2s


all_ini_Vals <- c(Vals, rep(1, p)) # with tau2s



#==================
# neg_logL function
#==================
source("Fn_neg_logL_CAR_2D.R")

neg_logL_CAR_2D(theta = all_ini_Vals, p = 5, 
                data_str = hierarchy_data_CAMS,
                all_pars_lst = all_pars_lst_CAR_2D_CMS,
                dsp_lon_mat = DSP[, , 1],
                dsp_lat_mat = DSP[, , 2],
                b = "Tri-Wave",
                phi = phi, H_adj = H_adj, 
                df = df_Lon_Strp_1_Srt)
 













