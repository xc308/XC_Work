## Matern function when nu = 3/2
# var = variance
# kappa = land scale
# d_vec = vectorzied distance 


Matern_32 <- function(Var, Kappa, d_vec) {
  M <- Var * ( 1 + Kappa * d_vec) * exp(-(Kappa * d_vec))
  matrix(M, nrow = sqrt(length(d_vec)))
}
