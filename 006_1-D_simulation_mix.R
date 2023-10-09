#=============
# 8 Oct. 2023
#=============

# Aim: 
  # C11 modelled using homogenous CAR
  # to show the mixed graphical model framework 
    # C11 is using CAR
    # while the partial cross-correlation (conditioned on the known)
    # b_lk is modelled using geostats method, 
  # such that both processes can be modelled within the same framework


# Methods:
  # two different Ajacency matrix: 
    # 1: spatial distance are ordered, using outer(s, s, "-")
    # 2: spatial distance hence Neighbour are not ordered, rnorm


# Sigma:
  # = (I - C)^{-1} M
  # C = phi * AJ
        # phi: spatial dependece parameter
        # AJ: Adjacency matrix corresponding to Neighbour
  # M = tau_sq * In
        # tau_sq: the V[Xl(si) | Xl(sj)], 
                  # different from V[Xl(si) | Xk(sj)]


#============
# 9 Oct. 2023
#============

# Aim:
# investigate the mixed feature in the mixed spatial graphical model
# when some of the univariate covriance is modelled using CAR
# when some of the them were modelled using Matern.

# Conditions:
# they are are grided data, 
# continuous lattice and discrete lattice
# each location then represents the centroid of the each grid


#---------
# settings
#---------

library(Matrix)




#--------------
# Data 1: ordered spatial seperation lag using outer(s, s, "-")
#--------------
ds <- 0.1
s <- seq(-1 + ds/2, 1 + ds/2, ds)


#----------------
# Ajacency matrix
#----------------

D <- abs(t(outer(s, s, "-")))

quantile(D)
#  0%  25%  50%  75% 100% 
# 0.0  0.3  0.6  1.0  1.9 

N_radius <- 0.5
N_index <- as.matrix(D <= N_radius)

Aj <- matrix(as.numeric(N_index), nrow(D), ncol(D))
diag(Aj) <- 0

str(Aj)


#-----------
# Parameters
#-----------

# C11 and C3_2: (I - phi * A)^{-1} * tau^2
set.seed(09-10)
#phi_11 <- round(rnorm(1), 1)
phi_11 <- -0.1
#phi_32 <- round(rnorm(1), 1)
phi_32 <- -0.3 
tau11_sq <- 1
tau32_sq <- 0.5


# C2_1 and C4_3: Matern
sig2_21 <- sig2_43 <- 1
kappa_21 <- kappa_43 <- 2



#-----------------------------
# C11  and Conditional variance
#-----------------------------

#~~~~~~~~~~~~~~~~~
# C11 and C3_2 CAR
#~~~~~~~~~~~~~~~~~
I_mat <- Diagonal(length(s))
mid <- solve(I_mat - phi_11 * Aj)
C11 <- tau11_sq * mid
#mid@p
# "p" attribute is an integer vector that stores the starting position of each column in the data arrays.

C3_2 <- tau32_sq * solve(I_mat - phi_32 * Aj)

Test_sym_pd(C11)
Test_sym_pd(C3_2)

#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"


#~~~~~~~
# Investigate the condition for CAR sigma to be p.d.
#~~~~~~~
# Conjecture:
  # phi has range: 1/lambda_1, 1/lambda_n
  # ordered eigen value of Aj matrix: lambda_1 < 0, lambda_n > 0


eigen_Aj <- eigen(Aj, only.values = T) 
range(eigen_Aj$values) # [1] -1.957639  9.625191

1 / range(eigen_Aj$values) 
# [1] -0.5108193  0.1038940
# the range for the choice of phi



#~~~~~~~~~~~~~~~~~~~~~
# C2_1 and C4_2 Matern
#~~~~~~~~~~~~~~~~~~~~~

source("Fn_Matern_32.R")
C2_1 <- Matern_32(Var = sig2_21, Kappa = kappa_21, d_vec = D_vec)
C4_3 <- Matern_32(Var = sig2_43, Kappa = kappa_43, d_vec = D_vec)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D : block diag containing all above
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- bdiag(C11, C2_1, C3_2, C4_3)

str(D)
# Formal class 'dgCMatrix' [package "Matrix"] with 6 slots



#===========
# Algorithm
#===========

p <- 4

#-----------------------------
# Hierarchical data structure
#-----------------------------

hierarchy_data_mix <- data.frame(
  node_id = c(1, 2, 3, 3, 4),
  par_id = c(NA, 1, c(2, 1), 1)
)


#--------
# Check parent node
#-------
source("Fn_Check_par_node.R")
Check_par_node(Node = 2, data = hierarchy_data_mix) # [1] 1
Check_par_node(Node = 3, data = hierarchy_data_mix) # [1] 2 1


#---------
# b() * ds
#---------

source("Fn_Wave_V4.R")

# displacements between pairs of points
H <- outer(df$s, df$s, FUN = "-")
H <- t(H)  

A <- 1
B <- A * wave_v4(h = H, delta = delta, A = 1) * ds
# [1:20, 1:20]


#----------------
# Construct SIGMA
#----------------

SIGMA_1212_Mix <- make_SIGMA(p = 2, data = hierarchy_data_mix)
SIGMA_1212_Mix <- as(SIGMA_1212_Mix, "matrix")
Test_sym_pd(SIGMA_1212_Mix)


SIGMA_1313_Mix <- make_SIGMA(p = 3, data = hierarchy_data_mix)
SIGMA_1313_Mix <- as(SIGMA_1313_Mix, "matrix")
Test_sym_pd(SIGMA_1313_Mix)

#[1] "Symmetric: Yes"
#[1] "p.d.: Yes"

SIGMA_1414_Mix <- make_SIGMA(p = 4, data = hierarchy_data_mix)
SIGMA_1414_Mix <- as(SIGMA_1414_Mix, "matrix")
Test_sym_pd(SIGMA_1414_Mix)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"



#----------
# Visualize Sigma with CAR and Matern Mixed
#----------
plt_Sig(Sigma = SIGMA_1212_Mix, p = 2)
plt_Sig(Sigma = SIGMA_1313_Mix, p = 3)
plt_Sig(Sigma = SIGMA_1414_Mix, p = 4)



