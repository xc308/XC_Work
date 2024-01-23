#=============
# 19 Jan. 2024
#=============

# Aim:
  # inference for 1D
  # 1. neg_logL function return neg log likelihood function
  # 2. optim such a neg_logL function to obtain the optimized parameters
  # 3. use the optimized parameters to cokrig true processes jointly


#==========
# settings
#==========

source("Fn_para_mat_construct.R")
#all_paras_mat <- All_paras(p = 7, data = hierarchy_data4)

all_pars_lst <- All_paras(p = 5, data = hierarchy_data)


# source("035_1D_simu_Matern_Chain_6F.R") 
# TST11_SpNormPert_SG_SGInv



ds <- 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20]

s <- seq(-10 + ds/2, 10 - ds/2, by = ds)


# displacements between pairs of points
# a vector quantity has magnitude and direction
H <- outer(s, s, FUN = "-")
H <- t(H)  
str(H) # num [1:20, 1:20]; num [1:200, 1:200]

# distance
# a scalar quantity
D_vec <- as.double(c(abs(H))) 
str(D_vec) # num [1:400]; num [1:40000]



#==========
# neg_logL
#==========
# Arguments:  
  # theta: a vector of all parameters
  # p: the number of variates
  # data_str: hierarchy data structure
  # all_pars_lst: a list of different parameter matrices
  # Z: observations for all processes


neg_logL <- function(theta, p, data_str, all_pars_lst, Z, obs_ind){
  
  #? Z1 ? 
  
  
  # connect each component of theta to all_pars_lst 
  # to incoporate each theta component into the neg log L function
 
  theta_indx <- 1
  for (lst in 1:length(all_pars_lst[[lst]])){
    for (i in 1:nrow(all_pars_lst[[lst]])){
      for (j in 1:ncol(all_pars_lst[[lst]])){
        if (is.na(all_pars_lst[[lst]][i, j])){
          all_pars_lst[[lst]][i, j] <- theta[theta_indx]
          theta_indx <- theta_indx + 1
        } 
      }
    }
  }
  
  
  # construct SIGMA_Y, SIGMA_Y_inv for process Y
  SG_SG_inv_Y <- TST11_SpNormPert_SG_SGInv(p = p, data = data_str, 
                            A_mat = all_pars_lst[[1]],
                            dlt_mat = all_pars_lst[[2]],
                            sig2_mat = all_pars_lst[[3]],
                            kappa_mat = all_pars_lst[[4]], 
                            d_vec = D_vec, h = H)
  
   
  SG_Y <- SG_SG_inv_Y$SIGMA
  SG_Y_inv <- SG_SG_inv_Y$SIGMA_inv
  
  # calculate SG_Ng
  ## 1st calcuate the # of parameters accumulated so far, 
    # so can connect theta components on top of current index
    # with measurement error tau2
  
  SUM <- 0 
  for (i in 1:length(all_pars_lst)){
    s <- sum(is.na(all_pars_lst[[i]]))
    SUM <- SUM + s
  }
  
  # tau2 diag matrix
  #tau2_mat <- diag(theta[SUM+1], theta[SUM+2], ..., theta[SUM+p] )
  THETA <- c()
  for(i in 1:p){
    THETA <- c(THETA, theta[SUM + i])
  }
  tau2_mat <- diag(THETA)
  
  #n1 <- length(Z1) # total # of locations of univariate process
  I_sp_mat <- I_sparse(size = n1, value = 1)
  SG_Ng <- kronecker(tau2_mat, I_sp_mat)
  SG_Ng_inv <- solve(SG_Ng)
  
  
  # SIGMA, SIGMA_inv for observation Z
  SG_Z = SG_Y + SG_Ng
  
  
  # SG_Z_inv = SG_Ng_inv - SG_Ng_inv(SG_Y_inv +SG_Ng_inv)^{-1}SG_Ng_inv
  SG_Y_Ng <- SG_Y_inv + SG_Ng_inv
  SG_Y_Ng_inv <- chol2inv(chol(SG_Y_Ng))
  
  SG_Z_inv <- SG_Ng_inv - SG_Ng_inv %*% SG_Y_Ng_inv %*% SG_Ng_inv
  
  
  # log_det(SG_Z)
  source("Fn_log_det.R")
  chol_SG_Z <- chol(SG_Z)
  log_SG_Z_det <- log_det(chol_SG_Z)
  
  
  # neg_logL
  L <- length(Z) # different from n1 = length(Z1)
  neg_logL <- - (- (L/2) * log(2*pi) - 1/2 * log_SG_Z_det - 
                   1/2 * t(Z) %*% SG_Z_inv %*% Z) 
  
  
  # return scalar
  return(neg_logL)
  
}



#========
# Optim
#========

# Aim:
  # optimize the above neg_logL function
  
# Method:
  # 1. optim function
  # 2. in which par aguement indicate the initial values for all parameters
      # is a vector, so need to collect matrices of parameters into a vector
      # according to their sequence. 

  # 3. initial value assign according to types such as A, dlt, sig2, kappa


#---------------------------------
# Initial values matrix to vector
#---------------------------------

ini <- c(1, 0.1, 1, 2) # A, dlt, sig2, kappa
Vals <- c()
for (i in 1:length(all_pars_lst)){
  value <- rep(ini[i], sum(is.na(all_pars_lst[[i]])))
  Vals <- c(Vals, value)
}

# Vals
#[1] 1.0 1.0 1.0 0.1 0.1 0.1
#[7] 1.0 1.0 1.0 2.0 2.0 2.0

## lower bound for each parameters, 
  # NA: no lower bound
  # kappa: must > 0
lower_bound <- c(rep(NA, sum(is.na(all_pars_lst[[1]]))),
                 rep(0.05, sum(is.na(all_pars_lst[[2]]))),
                 rep(0.001, sum(is.na(all_pars_lst[[3]]))),
                 rep(0.001, sum(is.na(all_pars_lst[[4]]))),
                 rep(0.001, p)
                 )
#  [1]    NA    NA    NA    NA
#  [5]    NA    NA 0.050 0.050
#  [9] 0.050 0.050 0.050 0.050
#  [13] 0.001 0.001 0.001 0.001
#  [17] 0.001    NA    NA    NA
#  [21]    NA    NA 0.001 0.001
#  [25] 0.001 0.001 0.001


optm_pars <- optim(par = c(Vals, rep(1, p)), # ini guess
      fn = neg_logL,
      method = "L-BFGS-B",
      lower = lower_bound,
      control = list(trace = 1, 
                     pgtol = 1e-5,
                     maxit = 3000))$par




#==========
# Cokrig
#==========

# Aim:
  # use the above optimzed parameters krig the desired 
  # multiple true processes simultaneously

# Method:
  # 1. construct SIGMA, SIGMA_inv for process 
      # using optimized parameters
  # 2. but TST function require parameter matrices
      # while optimized parameters are in vector form


#---------------------------
# optm_pars vector 2 matrix
#---------------------------
vec_2_mat <- function(vector, all_pars_lst) {
  
  vec_indx <- 1
  
  for(lst in 1:length(all_pars_lst)){
    for (i in 1:nrow(all_pars_lst[[lst]])){
      for (j in 1:ncol(all_pars_lst[[lst]])){
        if (is.na(all_pars_lst[[lst]][i, j])){
          all_pars_lst[[lst]][i, j] <- vector[vec_indx]
          vec_indx <- vec_indx + 1
        }  
      }
    }
  }
  return(all_pars_lst)
}


optm_pars_lst <- vec_2_mat(vector = optm_pars, all_pars_lst)


#-------------------------------------------------
# construct SIGMA, SIGMA_inv using optm_pars_lst
#-------------------------------------------------
optm_SG_SG_inv <- TST11_Matern_Chain(p = p, data = hierarchy_data, 
                   A_mat = optm_pars_lst[[1]],
                   dlt_mat = optm_pars_lst[[2]],
                   sig2_mat = optm_pars_lst[[3]],
                   kappa_mat = optm_pars_lst[[4]],
                   d_vec = D_vec, h = H)


SG_Y <- optm_SG_SG_inv$SIGMA
SG_Y_inv <- optm_SG_SG_inv$SIGMA_inv


#------------------
# SG_Ng, SG_Ng_inv 
#------------------
# optimized measurement error are the final p parameters
#Totl_Vals <- c(Vals, rep(1, p))
#Totl_Vals[-(1:length(Vals))]
tau2_optm <- optm_pars[-(1:length(Vals))]

tau2_optm_mat <- diag(tau2_optm)

n1 <- length(Z1)
I_mat <- I_sparse(size = n1, value = 1)
SG_Ng <- kronecker(tau2_optm_mat, I_mat)
SG_Ng_inv <- solve(SG_Ng)

## (SG_Y_inv + SG_Ng_inv)^{-1}
SG_Y_Ng <- SG_Y_inv + SG_Ng_inv
SG_Y_Ng_inv <- chol2inv(chol(SG_Y_Ng))

## SG_Z_inv
SG_Z_inv <- SG_Ng_inv - SG_Ng_inv %*% SG_Y_Ng_inv %*% SG_Ng_inv


#-------------------------------
# Cokrig true processes jointly
#-------------------------------

mu <- SG_Y %*% SG_Z_inv %*% Z
var <- diag(SG_Y - SG_Y %*% SG_Z_inv %*% SG_Y)


#----------------------------------
# save each true process separately
#----------------------------------













