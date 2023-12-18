#=============
# 18 Dec. 2023
#=============

# Aim:
  # 1. want to know what str does wave function returns
  # 2. want to know what str does B after spectral normalization returns

source("Fn_Waves.R")
source("Fn_Wendland_32.R")

wv <- wave_v5(h = H, delta = 0.5, A = 1)
str(wv)
# num [1:20, 1:20]
wv_d01 <- wave_v5(h = H, delta = 0.1, A = 1)

eigen(wv_d01, only.values = T)
# all zero
# singular

max(wv) # 1
min(wv) # -1

spectral_norm(wv)
# NaN or -Inf


#------------
# conclusion
#------------
# 1. since spectral normalization requires dividing
  # largest eigen value of B,
  # however, B must be low-rank hence singular
  # leading to 0 eigen values
  # therefore, scaling factor 1/max_singular must be
  # NaN or -Inf

# leading to the B_rt after spectral normalization 
  # becomes NaN or -Inf


#----------
# Remedies
#----------
# use all non-neg-valued B function, 
  # and then using row-stochastic 
  # i.e., each element of B divided by each row-sum

source("Fn_row_stoc.R")


#--------------------
# Check each row == 0
#--------------------
# 1. since to apply row-stochast to matrix
  # none of the rows of such a matrix 
  # must be all 0, 
  # as long as each row of the elements are non-neg
  # and not a row having all zero elements
  # can the matrix be apply row-stochast


# 2. The above wv value of B function always have
  # the last low being all zeor due to shift
  # consider this might also be the reason 
  # the B matrix being low rank 


# 3. Try to add the last value of the row with all-zero
  # elements 

check_set <- function(A_mat) {
  
  for(i in 1:nrow(A_mat)){
    if(all(A_mat[i, ] == 0)){
      #cat(i, "all zero", "\n")
      A_mat[i, ncol(A_mat)] <- A_mat[i-1, ncol(A_mat)]
    }
  }
  return(A_mat)
}

wv_set <- check_set(wv)

eigen(wv_set, only.values = T)$val
# the eigen value with the largest magnitude now is 0.28
# so can do spectral normalization

spectral_norm(wv_set) # normalizable


#--------------------------------------------
# Investigate which row has all zero elements
#--------------------------------------------
for (d in seq(0.1, 1, by = 0.1)){
  for (a in seq(0.1, 1, by = 0.1)){
    
    wv5 <- wave_v5(h = H, delta = d, A = a)
    
    cat("dlt:", d, "A:", a, "\n")
    check_set(wv5)
  }
}

# all the 20th row of wv5 has all-zero elements


#--------
# want to know if minor adjust the last element of 20th row
# will have normal largest singular value
# such that spectral normalization is viable
#--------
for (d in seq(0.1, 1, by = 0.1)){
  for (a in seq(0.1, 1, by = 0.1)){
    
    wv5 <- wave_v5(h = H, delta = d, A = a)
    
    cat("dlt:", d, "A:", a, "\n")
    wv5_set <- check_set(wv5)
    
    eig_val_wv5_set <- eigen(wv5_set, only.values = T)$val
    max_singula <- max(Mod(eig_val_wv5_set))
    cat("max singular:", max_singula, "\n")
    
    min_singula <- min(Mod(eig_val_wv5_set))
    cat("min singular:", min_singula, "\n")
  }
}

## Result:
  # 1. all max singular value are non-zero
  # 2. all min singular value are 0

  # 3. min singular value be 0 means Inf condition number
  # 4. computation involving such matrix would be very
    # numerical instable



#-------
# want to know the max and min singular value of B
# without normalization
# The problem might lie in the smallest eigen value being 0
#-------

for (d in seq(0.1, 1, by = 0.1)){
  for (a in seq(0.1, 1, by = 0.1)){
    
    wv5 <- wave_v5(h = H, delta = d, A = a)
    
    cat("dlt:", d, "A:", a, "\n")
    # set 
    wv5_set <- check_set(wv5)
    
    eig_val_wv5_set <- eigen(wv5_set, only.values = T)$val
    max_singula <- max(Mod(eig_val_wv5_set))
    cat("max singular:", max_singula, "\n")
    
    min_singula <- min(Mod(eig_val_wv5_set))
    cat("min singular:", min_singula, "\n")
    
    all_sig_val <- Mod(eig_val_wv5_set)
    cat("all singular:",all_sig_val, "\n")
    
    svd_d <- svd(eig_val_wv5_set)$d
    cat("svd all d:", svd_d, "\n")
    
    kappa_wv5_set <- kappa(eig_val_wv5_set)
    cat("C.N:", kappa_wv5_set, "\n")
  }
}

# Conclusion:
  # 1. the matrix B using wave function has only 1
  # non-zero singular value
  # the condition number of such a matrix is Inf

  # 2. SO, the problem is not lie in the B
  # have too large norm but
  # in the B's smallest eigen value is exact 0
  # such that the condition number of B is Inf

  # such that computation invloving such matrix B
  # would be very numerical unstable


# Remedies:
  # 1. spectral normalization
  # 2. svd reconstrut


#------------------------------------------------------------
# Function to merge both check&set and spectral nomalization
#------------------------------------------------------------

check_set_SpNorm <- function(A_mat) {
  # check each row of the matrix
  # if a row has all zero elements
  # set the last element of such a row 
  # with the last element of the last row
  # to ensure B is not low-rank or singular
  
  for(i in 1:nrow(A_mat)){
    if(all(A_mat[i, ] == 0)){
      A_mat[i, ncol(A_mat)] <- A_mat[i-1, ncol(A_mat)]
    }
  }
  #return(A_mat)
  
  # spectral normalization
  singular_val <- svd(A_mat)$d
  max_sg_val <- max(singular_val)
  scale_factor <- 1/max_sg_val
  
  # element-wise scaling
  normed_A <- scale_factor * A_mat
  return(normed_A)
}


wv_spnm_d05_a1 <- check_set_SpNorm(wv)
max(wv_spnm_d05_a1) # 0.2426644
min(abs(wv_spnm_d05_a1)) #0


# original
# dlt: 0.5 A: 1 
#max singular: 0.28 
#min singular: 0 
#all singular: 0.28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
#svd all d: 0.28 

svd(wv_spnm_d05_a1)$d
# [1] 1.000000e+00 9.346484e-01
#[3] 8.478413e-01 7.009933e-01
#[5] 4.738735e-01 2.671377e-01
#[7] 2.029344e-01 1.680954e-01
#[9] 1.414017e-01 1.280659e-01
#[11] 1.203791e-01 1.047980e-01
#[13] 8.790467e-02 8.368375e-02
#[15] 5.489871e-02 5.423650e-02
#[17] 4.592267e-02 4.558934e-02
#[19] 3.642992e-09 0.000000e+00

kappa(wv_spnm_d05_a1)
# [1] Inf
# still have zero and nearly zero svd$d

## remedies: Low rank approximation


#---------------------------
# SVD low-rank approximation
#---------------------------
svd_res <- svd(wv_spnm_d05_a1)
thres <- 1e-4

k <- which(svd_res$d > thres)
# [1]  1  2  3  4  5  6  7  8  9 10 11 12
# [13] 13 14 15 16 17 18

U_k <- svd_res$u[, k]
D_k <- diag(svd_res$d[k])
V_k <- svd_res$v[, k]

LRApp <- U_k %*% D_k %*% t(V_k)
str(LRApp)# num [1:20, 1:20]
max(LRApp) #1
min(abs(LRApp)) #0

svd(LRApp)$d
kappa(LRApp)

## results:
  # still have 0 or very close to 0 svd$d
  # each element of becomes very small
  # meaning of B might be not good


# set diag of LRApp to be 1
diag(LRApp) <- 1
svd(LRApp)$d
# [1] 1.8798341 1.5210528
#[3] 1.2047165 1.1179271
#[5] 1.1120803 1.1057883
#[7] 1.0945585 1.0834335
#[9] 1.0688879 1.0406467
#[11] 1.0012302 0.9995648
#[13] 0.9976593 0.9931084
#[15] 0.9753975 0.9654475
#[17] 0.9157302 0.8750985
#[19] 0.4759078 0.4502985
kappa(LRApp)
# [1] 5.387829


#------------------------------
# low rank approximation on wv
#------------------------------

wv_input <- check_set(wv)

svd_res <- svd(wv_input)

thres <- 1e-4

k <- which(svd_res$d > thres)

U_k <- svd_res$u[, k]
D_k <- diag(svd_res$d[k])
V_k <- svd_res$v[, k]

LRApp <- U_k %*% D_k %*% t(V_k)
# each element of reconstrcuted matrix B
# is very small
max(LRApp) # 1
min(abs(LRApp)) #0

kappa(LRApp) # [1] Inf




#-----------------
# Further remedies: regularize
#-----------------
wv
dim(wv_spnm_d05_a1)


diag(wv_spnm_d05_a1) <- 1


svd(wv_spnm_d05_a1)$d
# [1] 1.8736123 1.5148458 1.2010459 1.1092155
#[5] 1.1076608 1.1032075 1.0980463 1.0814592
#[9] 1.0598665 1.0492390 1.0349956 1.0290149
#[13] 1.0161879 1.0050470 0.9670579 0.9560416
#[17] 0.9046041 0.8683802 0.4654519 0.4457206

kappa(wv_spnm_d05_a1)
# [1] 5.420932


#===============
# Final strategy
#===============

# 1.check row with all zero elements
  # input with same position element of last row + 
# 2.spectral normalization +
# 3.diag element set to 1

# reason:
  # b_lk(s, s) or b_lk(h = 0) = 1

  # in addition, by structure of function b
  # b_lk(/h - dlt/) = 1

# Impact:
  # lower the condition number of B
  # via 1. spectral normalization cap B's 
            # norm to 1, 
      # 2. diag(B) set to 1, get rid off
            # the 0 and close to zero svd$d 


#--------------------------------------
# Want to know direct set B[20, 20] = 1
# then spectral norm
#--------------------------------------

wv[20, 20] <- 1
spN_wv <- SpNorm(wv)

svd(spN_wv)$d
# [1] 1.000000e+00 9.349626e-01 8.540419e-01
#[4] 7.093053e-01 4.842905e-01 2.838643e-01
#[7] 2.078682e-01 1.736915e-01 1.489516e-01
#[10] 1.306379e-01 1.255745e-01 1.187071e-01
#[13] 8.813919e-02 8.724716e-02 5.492000e-02
#[16] 5.446144e-02 4.592114e-02 4.572582e-02
#[19] 6.090703e-09 0.000000e+00

## conclusion:
  # first input 1 to wv[20, 20]
  # then spectral norm
  # still get 0 eigen values


#------------------------
# Final strategy function
#------------------------
check_set_SpNorm_pert <- function(A_mat) {
  # check each row of the matrix
  # if a row has all zero elements
  # set the last element of such a row 
  # with the last element of the last row
  # to ensure B is not low-rank or singular
  
  for(i in 1:nrow(A_mat)){
    if(all(A_mat[i, ] == 0)){
      A_mat[i, ncol(A_mat)] <- A_mat[i-1, ncol(A_mat)]
    }
  }
  
  
  # spectral normalization
  singular_val <- svd(A_mat)$d
  max_sg_val <- max(singular_val)
  scale_factor <- 1/max_sg_val
  
  # element-wise scaling
  normed_A <- scale_factor * A_mat
  
  # diag(norm_A) set to 1
  diag(normed_A) <- 1
  
  return(normed_A)
}


wv_remedy <- check_set_SpNorm_pert(wv)
kappa(wv_remedy)

wv_01_remedy <- check_set_SpNorm_pert(wv_d01)
kappa(wv_01_remedy)


wv_0101 <- wave_v5(h = H, delta = 0.1, A = 0.1)
wv_0101_remedy <- check_set_SpNorm_pert(wv_0101)
kappa(wv_0101_remedy)


#======================
# Test all combinations
#======================
# Aim:
  # each combination the condition number of B


for (d in seq(0.1, 1, by = 0.1)){
  for (a in seq(0.1, 1, by = 0.1)){
    
    wv5 <- wave_v5(h = H, delta = d, A = a)
    # set 
    wv5_remedy <- check_set_SpNorm_pert(wv5)
    wv5_spnorm <- check_set_SpNorm(wv5)
    
    cat("dlt:", d, "A:", a, "\n")
    cat("condition number:", kappa(wv5_remedy), "\n")
    
    cat("max singular val:", max(svd(wv5_remedy)$d), "\n")
    cat("min singular val:", min(svd(wv5_remedy)$d), "\n")
    
    
    cat("No pert", "\n")
    cat("condition number:", kappa(wv5_spnorm), "\n") # Inf
    cat("max singular val:", max(svd(wv5_spnorm)$d), "\n") # 1
    cat("min singular val:", min(svd(wv5_spnorm)$d), "\n") # 0
  }
}


# dlt: 0.1 A: 0.1 
#condition number: 9.111521 
#dlt: 0.1 A: 0.2 
#condition number: 9.111521 
#dlt: 0.1 A: 0.3 
#condition number: 9.111521 
#dlt: 0.1 A: 0.4 
#condition number: 9.111521 
#dlt: 0.1 A: 0.5 
#condition number: 9.111521 
#dlt: 0.1 A: 0.6 
#condition number: 9.111521 
#dlt: 0.1 A: 0.7 
#condition number: 9.111521 
#dlt: 0.1 A: 0.8 
#condition number: 9.111521 
#dlt: 0.1 A: 0.9 
#condition number: 9.111521 
#dlt: 0.1 A: 1 
#condition number: 9.111521 
#dlt: 0.2 A: 0.1 
#condition number: 5.458719 
#dlt: 0.2 A: 0.2 
#condition number: 5.458719 
#dlt: 0.2 A: 0.3 
#condition number: 5.458719 
#dlt: 0.2 A: 0.4 
#condition number: 5.458719 
#dlt: 0.2 A: 0.5 
#condition number: 5.458719 
#dlt: 0.2 A: 0.6 
#condition number: 5.458719 
#dlt: 0.2 A: 0.7 
#condition number: 5.458719 
#dlt: 0.2 A: 0.8 
#condition number: 5.458719 
#dlt: 0.2 A: 0.9 
#condition number: 5.458719 
#dlt: 0.2 A: 1 
#condition number: 5.458719 
#dlt: 0.3 A: 0.1 
#condition number: 6.826081 
#dlt: 0.3 A: 0.2 
#condition number: 6.826081 
#dlt: 0.3 A: 0.3 
#condition number: 6.826081 
#dlt: 0.3 A: 0.4 
#condition number: 6.826081 
#dlt: 0.3 A: 0.5 
#condition number: 6.826081 
#dlt: 0.3 A: 0.6 
#condition number: 6.826081 
#dlt: 0.3 A: 0.7 
#condition number: 6.826081 
#dlt: 0.3 A: 0.8 
#condition number: 6.826081 
#dlt: 0.3 A: 0.9 
#condition number: 6.826081 
#dlt: 0.3 A: 1 
#condition number: 6.826081 
#dlt: 0.4 A: 0.1 
#condition number: 5.499038 
#dlt: 0.4 A: 0.2 
#condition number: 5.499038 
#dlt: 0.4 A: 0.3 
#condition number: 5.499038 
#dlt: 0.4 A: 0.4 
#condition number: 5.499038 
#dlt: 0.4 A: 0.5 
#condition number: 5.499038 
#dlt: 0.4 A: 0.6 
#condition number: 5.499038 
#dlt: 0.4 A: 0.7 
#condition number: 5.499038 
#dlt: 0.4 A: 0.8 
#condition number: 5.499038 
#dlt: 0.4 A: 0.9 
#condition number: 5.499038 
#dlt: 0.4 A: 1 
#condition number: 5.499038 
#dlt: 0.5 A: 0.1 
#condition number: 5.420932 
#dlt: 0.5 A: 0.2 
#condition number: 5.420932 
#dlt: 0.5 A: 0.3 
#condition number: 5.420932 
#dlt: 0.5 A: 0.4 
#condition number: 5.420932 
#dlt: 0.5 A: 0.5 
#condition number: 5.420932 
#dlt: 0.5 A: 0.6 
#condition number: 5.420932 
#dlt: 0.5 A: 0.7 
#condition number: 5.420932 
#dlt: 0.5 A: 0.8 
#condition number: 5.420932 
#dlt: 0.5 A: 0.9 
#condition number: 5.420932 
#dlt: 0.5 A: 1 
#condition number: 5.420932 
#dlt: 0.6 A: 0.1 
#condition number: 5.420437 
#dlt: 0.6 A: 0.2 
#condition number: 5.420437 
#dlt: 0.6 A: 0.3 
#condition number: 5.420437 
#dlt: 0.6 A: 0.4 
#condition number: 5.420437 
#dlt: 0.6 A: 0.5 
#condition number: 5.420437 
#dlt: 0.6 A: 0.6 
#condition number: 5.420437 
#dlt: 0.6 A: 0.7 
#condition number: 5.420437 
#dlt: 0.6 A: 0.8 
#condition number: 5.420437 
#dlt: 0.6 A: 0.9 
#condition number: 5.420437 
#dlt: 0.6 A: 1 
#condition number: 5.420437 
#dlt: 0.7 A: 0.1 
#condition number: 5.406821 
#dlt: 0.7 A: 0.2 
#condition number: 5.406821 
#dlt: 0.7 A: 0.3 
#condition number: 5.406821 
#dlt: 0.7 A: 0.4 
#condition number: 5.406821 
#dlt: 0.7 A: 0.5 
#condition number: 5.406821 
#dlt: 0.7 A: 0.6 
#condition number: 5.406821 
#dlt: 0.7 A: 0.7 
#condition number: 5.406821 
#dlt: 0.7 A: 0.8 
#condition number: 5.406821 
#dlt: 0.7 A: 0.9 
#condition number: 5.406821 
#dlt: 0.7 A: 1 
#condition number: 5.406821 
#dlt: 0.8 A: 0.1 
#condition number: 5.250274 
#dlt: 0.8 A: 0.2 
#condition number: 5.250274 
#dlt: 0.8 A: 0.3 
#condition number: 5.250274 
#dlt: 0.8 A: 0.4 
#condition number: 5.250274 
#dlt: 0.8 A: 0.5 
#condition number: 5.250274 
#dlt: 0.8 A: 0.6 
#condition number: 5.250274 
#dlt: 0.8 A: 0.7 
#condition number: 5.250274 
#dlt: 0.8 A: 0.8 
#condition number: 5.250274 
#dlt: 0.8 A: 0.9 
#condition number: 5.250274 
#dlt: 0.8 A: 1 
#condition number: 5.250274 
#dlt: 0.9 A: 0.1 
#condition number: 4.767257 
#dlt: 0.9 A: 0.2 
#condition number: 4.767257 
#dlt: 0.9 A: 0.3 
#condition number: 4.767257 
#dlt: 0.9 A: 0.4 
#condition number: 4.767257 
#dlt: 0.9 A: 0.5 
#condition number: 4.767257 
#dlt: 0.9 A: 0.6 
#condition number: 4.767257 
#dlt: 0.9 A: 0.7 
#condition number: 4.767257 
#dlt: 0.9 A: 0.8 
#condition number: 4.767257 
#dlt: 0.9 A: 0.9 
#condition number: 4.767257 
#dlt: 0.9 A: 1 
#condition number: 4.767257 
#dlt: 1 A: 0.1 
#condition number: 4.189142 
#dlt: 1 A: 0.2 
#condition number: 4.189142 
#dlt: 1 A: 0.3 
#condition number: 4.189142 
#dlt: 1 A: 0.4 
#condition number: 4.189142 
#dlt: 1 A: 0.5 
#condition number: 4.189142 
#dlt: 1 A: 0.6 
#condition number: 4.189142 
#dlt: 1 A: 0.7 
#condition number: 4.189142 
#dlt: 1 A: 0.8 
#condition number: 4.189142 
#dlt: 1 A: 0.9 
#condition number: 4.189142 
#dlt: 1 A: 1 
#condition number: 4.189142 


#===========
# Think: what if we only set the diag(B) = 1
#===========
# diag(B) = 1 + 
# spectral normalization

diag(wv) <- 1 
svd(wv)$d
# [1] 4.87334711 4.05444952
#[3] 3.15994665 2.35838593
#[5] 1.87478950 1.48305480
#[7] 1.45831815 1.45738450
#[9] 1.45501534 1.45342429
#[11] 1.43128852 1.25259172
#[13] 1.15245117 1.13645990
#[15] 1.09385360 1.05613285
#[17] 0.89515767 0.82886789
#[19] 0.02394332 0.01127690

kappa(wv)
# [1] 676.5478

spN_wv <- SpNorm(wv)
svd(spN_wv)$d
kappa(spN_wv)
# [1] 676.5478

## conclusion:
  # this stragey does not change kappa






